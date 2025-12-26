#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------- #
#  Unified Electro–Optical Transient Deconvolution
#  ----------------------------------------------
#  • Separates surface-limited and diffusion-limited charge storage from
#    chrono-amperometric data **and** links them to time-resolved UV-Vis OD.
#  • Three-step kinetic extraction (ln I fit → pre-fit → global fit).
#  • Generates three publication-ready figures:
#        1)  I(t) + fitted components + parameter table
#        2)  Accumulated charge components Q(t)
#        3)  Optical ΔOD(t) components
#
#  Author     :  Luhan WEI (Solid State Ionics Lab @ Westlake University) · 30-May-2025

#  License    :  MIT  – see text at end of file
#  Cite this  :  Luhan Wei, Yang Hu, Yiwei Huang, Ying Lu, Zihan Xu, Nian Zhang and Qiyang Lu*, 
#              Deconvoluting Surface and Bulk Charge Storage Processes in Redox-Active Oxides by 
#              Integrating Electrochemical and Optical Insights, 
#              Journal of the American Chemical Society, 146 (2024), 24167–24176
#              https://pubs.acs.org/doi/10.1021/jacs.4c09261
#
#  Python Version : 3.9.7
# --------------------------------------------------------------------------- #
"""
Quick Start
===========

1.  Install dependencies::

        pip install numpy pandas matplotlib scipy

2.  Point *both* paths below to your csv files:

        CSV_PATH     – current-vs-time data  (cols: time/s , I/mA)
        CSV_OD_PATH  – UV–Vis ΔOD time trace (cols: time/s , OD)

3.  Run::

        Unified Electro–Optical Transient Deconvolution.py
"""
# --------------------------------------------------------------------------- #
# Imports
# --------------------------------------------------------------------------- #
from __future__ import annotations

from math import pi
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import special
from scipy.optimize import curve_fit, fsolve

# --------------------------------------------------------------------------- #
# User-configurable constants
# --------------------------------------------------------------------------- #
CSV_PATH: str = (
    r"path\to\Anodic_Current_vs_Time.csv"              # ← current trace
)
CSV_OD_PATH: str = (
    r"path\to\Anodic_Optical_Density_vs_Time.csv"      # ← UV–Vis trace
)


L_CM: float          = 8e-5     # Film thickness  [cm]   (650 nm)
T_SHORT_MAX: float   = 60.0     # Upper time bound for Step-2/3 & plots
T_LIN_REGION: tuple  = (20.0, 80.0)  # ln I linear window on full curve

# Aesthetics — tweak if desired
FIG_SIZE = (9, 6)
FONT_SZ  = 16
MARKER_SZ = 12
LW_MAIN   = 3


# --------------------------------------------------------------------------- #
# Helper classes & util functions
# --------------------------------------------------------------------------- #
class RootFinder:
    """
    Lightweight root finder for ``x·tan(x) – Bi = 0`` used to obtain Λ₁.

    Parameters
    ----------
    start, stop : float
        Search interval.
    step        : float
        Increment for initial guesses.
    xtol        : float
        Solver tolerance.
    """

    def __init__(self,
                 start: float,
                 stop: float,
                 step: float = 0.01,
                 xtol: float = 1e-9) -> None:
        self.start, self.stop, self.step, self.xtol = start, stop, step, xtol
        self.roots: np.ndarray = np.array([], dtype=float)

    def _add(self, root: float) -> None:
        """Add root if unique and inside range."""
        if (self.start <= root <= self.stop
                and not np.any(np.abs(self.roots - root) < self.xtol)):
            self.roots = np.append(self.roots, root)

    def _solve(self,
               func: Callable,
               guess: float,
               *args) -> float | None:
        """Single call to ``fsolve`` with a starting guess."""
        root, *_ , flag, _ = fsolve(func,
                                    guess,
                                    args=args,
                                    full_output=True,
                                    xtol=self.xtol)
        return root[0] if flag == 1 else None

    def find(self,
             func: Callable,
             *args) -> np.ndarray:
        """Return all roots found in ``[start, stop]``."""
        current = self.start
        for g in np.arange(self.start, self.stop + self.step, self.step):
            if g < current:
                continue
            root = self._solve(func, g, *args)
            if root is not None:
                current = root
                self._add(root)
        return self.roots


def linear_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """y = a·x + b — for ln I regression."""
    return a * x + b


def fill_nan_neighbors(arr: np.ndarray) -> np.ndarray:
    """
    Replace NaNs by the mean of nearest non-NaN neighbors *in-place*.

    If one side is missing, fall back to the available neighbor.
    """
    out = arr.copy()
    for idx in np.where(np.isnan(out))[0]:
        left, right = idx - 1, idx + 1
        while left >= 0 and np.isnan(out[left]):
            left -= 1
        while right < out.size and np.isnan(out[right]):
            right += 1
        if 0 <= left < out.size and right < out.size:
            out[idx] = 0.5 * (out[left] + out[right])
        elif 0 <= left < out.size:
            out[idx] = out[left]
        elif right < out.size:
            out[idx] = out[right]
    return out


def cumulative_trapz(t: np.ndarray,
                     I: np.ndarray) -> np.ndarray:
    """
    Cumulative trapezoidal ∫I·dt.  Return Q(t) with Q(0) = 0 mC.
    """
    dq = 0.5 * (I[1:] + I[:-1]) * (t[1:] - t[:-1])
    return np.insert(np.cumsum(dq), 0, 0.0)


def Charge_Integration(x_data, y_data):
    
    dQ_lt = []
    
    x = np.array(x_data)
    
    y = np.array(y_data)
        
    for j in range(1, len(x)):
        
        Up_edge = y[j-1]
        
        Low_edge = y[j]
        
        H = x[j] - x[j-1]
        
        dQ = ((Up_edge + Low_edge) * H)/2
        
        dQ_lt.append(dQ)
     
    
    Q_Exp = []
    
    for j in range(0, len(dQ_lt)):
        
        sub_dQ = np.array(dQ_lt[:j])
        
        Q_Exp.append(sum(sub_dQ))
    
    return np.array(Q_Exp)

# --------------------------------------------------------------------------- #
# Model definitions (original formulas, untouched)
# --------------------------------------------------------------------------- #
def edlc_pseudo(tau_d: float,
                t: np.ndarray,
                i0: float,
                tau_edl: float,
                a: float,
                lamda: float,
                C: float) -> np.ndarray:
    """Short-time analytical expression (pre-fit)."""
    y_surface = i0 * np.exp(-t / tau_edl)
    y_bulk = (lamda * (a / tau_d)
              * np.exp(lamda**2 * t / tau_d)
              * special.erfc(lamda * np.sqrt(t / tau_d))
              + C)
    return y_surface + y_bulk


def edlc_common_k(Q_tot: float,
                  t: np.ndarray,
                  i0: float,
                  k: float,
                  l_surf: float,
                  D: float,
                  C: float) -> np.ndarray:
    """Global model with common surface rate constant k."""
    y_surface = i0 * np.exp(-t * (k / l_surf))
    y_bulk = (Q_tot * (k / L_CM)
              * np.exp(k ** 2 * t / D)
              * special.erfc(k * np.sqrt(t / D))
              + C)
    return y_surface + y_bulk


# --------------------------------------------------------------------------- #
# Main routine
# --------------------------------------------------------------------------- #
def main() -> None:
    # --- 1 · Load & sanitize ------------------------------------------------
    df = (pd.read_csv(CSV_PATH)
            .loc[lambda d: (d["time/s"] > 0) & (d["I/mA"] > 0)])
    
    t_full = df["time/s"].values
    I_full = df["I/mA"].values
    
    # --- 2 · Step-1 : ln I linear fit (20-80 s on *full* trace) ------------
    lin_mask = (T_LIN_REGION[0] <= t_full) & (t_full <= T_LIN_REGION[1])
    slope = curve_fit(linear_model,
                      t_full[lin_mask],
                      np.log(fill_nan_neighbors(I_full[lin_mask])))[0][0]
    tau_d = 1.0 / (abs(slope) * 4 / pi**2)
    D_linear = L_CM**2 / tau_d
    
    # --- 3 · Restrict to ≤ 60 s for Steps 2 & 3 ----------------------------
    short_mask = t_full <= T_SHORT_MAX
    t = t_full[short_mask]
    I = I_full[short_mask]
    Q_tot = cumulative_trapz(t, I).max()
    
    # --- 4 · Step-2 : pre-fit ---------------------------------------------
    bounds = ([0, 0, 0, 0, -0.5],
              [1, np.inf, np.inf, np.inf, np.inf])
    
    p_pre, _ = curve_fit(
        lambda tt, i0, tau_edl, a, lam, C:
            edlc_pseudo(tau_d, tt, i0, tau_edl, a, lam, C),
        t, I, bounds=bounds, maxfev=100_000)
    
    i0_init, tau_edl_init, a_init, lam_init, C_init = p_pre
    k_init = lam_init * L_CM / tau_d
    L_surf_init = k_init * tau_edl_init
    
    # --- 5 · Step-3 : global fit ------------------------------------------
    p_opt, p_cov = curve_fit(
        lambda tt, i0, k, l_surf, D, C:
            edlc_common_k(Q_tot, tt, i0, k, l_surf, D, C),
        t, I,
        p0=[i0_init, k_init, L_surf_init, D_linear, C_init],
        maxfev=100_000)
    
    i0, k, L_surf_cm, D_fit, C_res = p_opt
    
    # --- 6 · Derived quantities ------------------------------------------
    L_surf_nm = L_surf_cm * 1e7
    sub_L_cm  = (np.sqrt(t.max()) / np.sqrt(tau_d)) * L_CM
    D_mod     = sub_L_cm**2 / (L_CM**2 / D_fit)
    Bi        = k * sub_L_cm / D_mod
    
    # Λ₁ (kept for completeness)
    Λ1 = RootFinder(0.01, 5_000).find(lambda x, l=Bi: x * np.tan(x) - l)[0]
    D_linear_mod = -(slope * L_CM**2) / Λ1**2  # stored, not plotted
    
    # --- 7 · Build fitted curves / charges -------------------------------
    t_dense = np.linspace(t.min(), T_SHORT_MAX, len(t))
    
    I_surface = i0 * np.exp(-t_dense * (k / L_surf_cm))
    I_bulk    = edlc_common_k(Q_tot, t_dense, 0, k, L_surf_cm, D_fit, C_res)
    I_total   = I_surface + I_bulk
    
    Q_exp = Charge_Integration(t, I)
    Q_surf = Charge_Integration(t_dense, I_surface)
    Q_bulk = Charge_Integration(t_dense, I_bulk)
    Q_total = Charge_Integration(t_dense, I_total)
    
    Q_fit_time = t_dense[1:]
    
    # 8 · Load UV–Vis data  (0 – 60 s)  & fit ΔOD(t)
    # ------------------------------------------------------------------ #
    df_uv = (pd.read_csv(CSV_OD_PATH)
                .loc[lambda d: (d["time/s"] >= 0) & (d["time/s"] <= T_SHORT_MAX)])
    
    uv_t  = df_uv["time/s"].values
    uv_od = df_uv.iloc[:, 1].values            # 2nd column = OD(t)
    
    uv_exp = uv_od - min(uv_od)           # legacy normalisation
    
    f_surf = np.interp(uv_t, Q_fit_time, Q_surf)
    f_bulk = np.interp(uv_t, Q_fit_time, Q_bulk)
    
    def F(x, a, b):
        
        return a * f_surf +  b * f_bulk
    
    popt, pcov = curve_fit(F, uv_t, uv_exp)
    
    a = popt[0]
    b = popt[1]
    
    UV_Surf = a * f_surf
    UV_Bulk = b * f_bulk
    UV_Total = a * f_surf + b * f_bulk
    
    
    # --- 8 · Figure 1 : current components -------------------------------
    plt.figure(figsize=FIG_SIZE)
    plt.semilogx(t, I, "o", ms=MARKER_SZ, mfc="white", mew=2,
                 color="#4b5563", label="Experiment")
    plt.semilogx(t_dense, I_total,  "-", lw=LW_MAIN, color="#ec5156",
                 label=r"$I_{\mathrm{Total}}$")
    plt.semilogx(t_dense, I_surface, "-", lw=LW_MAIN, color="#1f77b4",
                 label=r"$I_{\mathrm{Surface}}$")
    plt.semilogx(t_dense, I_bulk,   "--", lw=LW_MAIN, color="#eebb65",
                 label=r"$I_{\mathrm{Bulk}}$")
    
    plt.xlabel("Time / s", fontsize=FONT_SZ)
    plt.ylabel("Current / mA", fontsize=FONT_SZ)
    plt.legend(frameon=False, fontsize=FONT_SZ)
    
    # Embedded table ------------------------------------------------------
    table_vals = [[f"{k:.2e}",      "cm/s"],
                  [f"{D_mod:.2e}",  "cm$^{2}$/s"],
                  [f"{L_surf_nm:.1f}", "nm"],
                  [f"{Bi:.2f}",     "–"]]
    
    tbl = plt.table(
        cellText=table_vals,
        rowLabels=["$\itk$", "$\itD$", "$\itL$$_{Surf}$", "$\it{Bi}$"],
        colLabels=["Value", "Unit"],
        loc="upper right",
        bbox=[0.45, 0.25, 0.50, 0.35])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(FONT_SZ)
    
    plt.tight_layout()
    plt.show()
    
    # --- 9 · Figure 2 : accumulated charge -------------------------------
    plt.figure(figsize=FIG_SIZE)
    plt.plot(Q_fit_time, Q_exp, "o", ms=MARKER_SZ, mfc="white", mew=2,
             color="#669ca8", label="Experiment")
    plt.plot(Q_fit_time, Q_total, "-", lw=LW_MAIN, color="#ec5156",
             label=r"$Q_{\mathrm{Total}}$")
    plt.plot(Q_fit_time, Q_surf,  "-", lw=LW_MAIN, color="#1f77b4",
             label=r"$Q_{\mathrm{Surface}}$")
    plt.plot(Q_fit_time, Q_bulk,  "--", lw=LW_MAIN, color="#eebb65",
             label=r"$Q_{\mathrm{Bulk}}$")
    
    plt.xlabel("Time / s", fontsize=FONT_SZ)
    plt.ylabel("Accumulated charge / mC", fontsize=FONT_SZ)
    plt.legend(frameon=False, fontsize=FONT_SZ)
    plt.tight_layout()
    plt.show()
    
    # ------------------------------------------------------------------ #
    # 11 · Figure 3 : UV–Vis ΔOD components
    # ------------------------------------------------------------------ #
    plt.figure(figsize=(9.5, 6))
    plt.plot(uv_t, uv_exp, "o", markersize=10,
              markeredgecolor="#67a3b3", markeredgewidth=2,
              markerfacecolor="#b3d7e0", alpha=0.8,
              label=r"ΔOD$_{Exp}$")
    plt.plot(uv_t, UV_Surf,  "-", lw=3.5, color="#1588d7",
              label=r"ΔOD$_{Surf}$")
    plt.plot(uv_t, UV_Bulk,  "-", lw=3.5, color="#ec9f0f",
              label=r"ΔOD$_{Bulk}$")
    plt.plot(uv_t, UV_Total, "-", lw=4.5, color="#7f7f7f",
              label=r"ΔOD$_{Total}$")
    
    plt.xlabel("Time / s", fontsize=FONT_SZ)
    plt.ylabel("$\Delta$OD (a.u.)", fontsize=FONT_SZ)
    plt.legend(frameon=False, fontsize=18)
    ax3 = plt.gca()
    for spine in ax3.spines.values():
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.show()



# --------------------------------------------------------------------------- #
# Command-line entry-point
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()

# --------------------------------------------------------------------------- #
# MIT License

# Copyright (c) 2025 Archer_LH (Luhan WEI)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------- #

