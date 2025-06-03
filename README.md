# Current–Transient & Optical Deconvolution Toolkit

A single-file Python demo that reproduces the electrochemical/optical analysis workflow reported in Wei et al. “Deconvoluting Surface and Bulk Charge-Storage Processes in Redox-Active Oxides by Integrating Electrochemical and Optical Insights” J. Am. Chem. Soc. 146 (2024) 24167-24176 (https://pubs.acs.org/doi/10.1021/jacs.4c09261). 

The script fits chrono-amperometric (CA) transients and operando ΔOD traces with one unified kinetic model, delivering the key parameters. 

(Note. This script is built on Python 3.9.7 Version)

---

## Why do we need this?

Redox-active transition-metal oxides (TMOs) show two concurrent charge-storage modes: a fast, near-surface pseudocapacitive layer and a slower, diffusion-limited bulk intercalation zone 
. Traditional electrochemical models seldom capture both length-scales. Wei et al. introduced an explicit surface-layer thickness L<sub>Surf</sub> and showed that the same parameter set fits CA and ΔOD simultaneously 
, unlocking truly quantitative surface/bulk deconvolution.

| mode | length‑scale | behaviour |
|------|--------------|-----------|
| surface pseudocapacitance | ≤ tens nm | fast, capacitive |
| bulk intercalation        | film thickness | slow, diffusive |

Wei *et al.* introduced an explicit surface layer of thickness *L*<sub>Surf</sub> and showed that a unified model explains *both* electrochemistry and operando spectroscopy.

---

## Repository layout

```
.
├─ src/
│   └─ current_transient_deconvolution.py   ← main demo script
├─ data/
│   ├─ Echem_Current/                       ← I(t)  @ 0.80 V vs Ag/AgCl
│   └─ Optical_Delta-OD/                    ← ΔOD(t) @ 0.80 V at the wavelength of 500 nm
└─ LICENSE
```

The bundled datasets reproduce Figure 5a of the paper for birnessite δ‑MnO2 (650 nm thick).

---

## Installation

```bash
pyenv install 3.9.7
pip install numpy pandas matplotlib scipy
```

---

## Running the demo

```bash
python src/Unified Electro–Optical Transient Deconvolution Code.py
```

Three figures appear:

1. **Current components** with parameter table (*k*, *D*, *L*<sub>Surf</sub>, *Bi*)
2. **Accumulated charge** curves Δ*Q*<sub>Surf/Bulk/Total</sub>
3. **ΔOD components** reconstructed from electrochemical fit

---

## Method summary

| Step | Time window | Purpose | Output |
|------|-------------|---------|--------|
| **1** | 20–80 s (Long-time estimation) | ln *I* fit | τ<sub>d</sub>, *D*<sub>linear</sub> |
| **2** | 0–60 s | pre‑fit surface+diffusion expression | initial guesses |
| **3** | 0–60 s | global fit, common *k* | *k*, *D*, *L*<sub>Surf</sub>, *Bi* |

ΔOD is then described by  

> ΔOD(*t*) = ε<sub>Surf</sub> Δ*Q*<sub>Surf</sub>(*t*) + ε<sub>Bulk</sub> Δ*Q*<sub>Bulk</sub>(*t*)

---

## Expected output (demo data)

| *k* (cm s‑¹) | *D* (cm² s‑¹) | *L*<sub>Surf</sub> (nm) | *Bi* |
|--------------|---------------|-------------------------|------|
| ~1E-5 | ~1E-11 | 5–25 | ~40 |

Matches Table 1 of Wei *et al.*.

---

## Using your own data

* Replace paths at the top of `current_transient_deconvolution.py`.
* Ensure electrochemical file has `time/s, I/mA`; UV‑Vis file `time/s, OD`.
* Adjust `L_CM`, `T_SHORT_MAX`, `T_LIN_REGION` if necessary.

---

## Citation

```text
@software{deconvolution_toolkit_2025,
  author  = {Rolen and ChatGPT-o3},
  title   = {Current–Transient & Optical Deconvolution Toolkit},
  year    = {2025},
  url     = {https://github.com/your‑org/ca‑opt‑deconvolution},
  note    = {MIT License}
}
```

and  

```text
@article{wei2024deconvoluting,
  title   = {Deconvoluting Surface and Bulk Charge Storage Processes in Redox-Active Oxides by Integrating Electrochemical and Optical Insights},
  author  = {Wei, Luhan and Hu, Yang and Huang, Yiwei *et al.*},
  journal = {J. Am. Chem. Soc.},
  volume  = {146},
  pages   = {24167--24176},
  year    = {2024},
  doi     = {10.1021/jacs.4c09261}
}
```

---

## License

Released under the MIT License – see `LICENSE` for details.


## Reference
[1] Li, J.; Xiao, X.; Yang, F.; Verbrugge, M. W.; Cheng, Y.-T. Potentiostatic Intermittent Titration Technique for Electrodes Governed by Diffusion and Interfacial Reaction. J. Phys. Chem. C 2012, 116, 1472– 1478,  DOI: 10.1021/jp207919q
[2] Montella, C. Discussion of the Potential Step Method for the Determination of the Diffusion Coefficients of Guest Species in Host Materials. Part I. Influence of Charge Transfer Kinetics and Ohmic Potential Drop. J. Electroanal. Chem. 2002, 518, 61– 83,  DOI: 10.1016/S0022-0728(01)00691-X
[3] Zhang, D.; Wang, R.; Wang, X.; Gogotsi, Y. In situ Monitoring Redox Processes in Energy Storage Using UV-Vis Spectroscopy. Nat. Energy 2023, 8, 567– 576,  DOI: 10.1038/s41560-023-01240-9
[4] Li, Y.; Chueh, W. C. Electrochemical and Chemical Insertion for Energy Transformation and Switching. Annu. Rev. Mater. Res. 2018, 48, 137– 165,  DOI: 10.1146/annurev-matsci-070317-124525



