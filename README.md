# Current–Transient & Optical Deconvolution Toolkit

A single-file Python toolkit for reproducing the electrochemical/optical analysis workflow reported in **Wei *et al.*, "Deconvoluting Surface and Bulk Charge-Storage Processes in Redox-Active Oxides by Integrating Electrochemical and Optical Insights," *J. Am. Chem. Soc.* (2024)** (https://doi.org/10.1021/jacs.4c09261).

This script fits chrono-amperometric (CA) transients and operando ΔOD traces using a unified kinetic model. It extracts key kinetic parameters, including the surface reaction rate constant ($k_{Surf}$) and the bulk chemical diffusion coefficient ($D_{Chem}$), by explicitly accounting for a surface redox layer ($L_{Surf}$).

---

## Overview

Redox-active transition metal oxides (TMOs) play central roles in electrochemical energy storage and conversion devices, where charge storage arises from coupled surface pseudocapacitive-like and bulk battery-type redox processes. Quantitative deconvolution of these processes is essential for understanding their kinetics and optimizing device performance, yet it remains challenging due to their distinct time and length scales.

In this work, we develop an electrochemical model that explicitly incorporates a surface redox layer thickness, denoted as $L_{Surf}$, to account for the different spatial extents of surface and bulk processes. By integrating this length-scale parameter with the surface reaction rate constant $k_{Surf}$ and the chemical diffusion coefficient $D_{Chem}$, the model enables a unified description of redox kinetics across multiple timescales. Using birnessite $\delta$-MnO2 as a representative redox-active TMO, we combine potential step chronoamperometry with operando ultraviolet-visible (UV-Vis) spectroscopy to capture transient electrochemical and optical responses associated with surface and bulk redox reactions.

We show that both the electrochemical currents and optical absorbance transients at each applied potential can be quantitatively reconstructed using a single set of kinetic parameters ($k_{Surf}$, $D_{Chem}$) together with $L_{Surf}$. This framework enables the explicit separation of surface and bulk contributions and reveals their distinct kinetic behaviors over multiple timescales. The extracted potential-dependent $D_{Chem}$, $k_{Surf}$, and $L_{Surf}$ provide a quantitative picture of how redox kinetics evolve during charge storage. By capturing both temporal and spatial aspects of redox processes, this methodology offers new insights for the rational design of advanced electrochemical energy storage and conversion systems.

| Mode | Length‑scale | Behavior |
| :--- | :--- | :--- |
| **Surface Pseudocapacitor-Like** | $\le$ tens of nm ($L_{Surf}$) | Fast, Capacitive |
| **Bulk Battery-Like** | Film thickness | Slow, Diffusive |

---

## Repository Layout

```text
.
├─ src/
│   └─ Unified Electro–Optical Transient Deconvolution Code.py    ← Main demo script
├─ data/
│   ├─ Echem_Current/                        ← I(t) @ 0.80 V vs Ag/AgCl
│   └─ Optical_Delta-OD/                     ← ΔOD(t) @ 0.80 V (500 nm wavelength)
└─ LICENSE

```

*The bundled datasets reproduce Figure 5a of the associated paper in https://doi.org/10.1021/jacs.4c09261.*

---

## Installation

The toolkit is validated on **Python 3.9** (e.g., 3.9.7). It is recommended to create a specific environment to avoid version conflicts.

```bash
# Example using pyenv
pyenv install 3.9.7

# Install dependencies
pip install numpy pandas matplotlib scipy

```

*(Note: Other Python versions have not been strictly validated and may cause numerical or plotting inconsistencies.)*

---

## Running the Demo

1. Open `Unified Electro–Optical Transient Deconvolution Code.py`.
2. Locate the **User-configurable constants** section.
3. Ensure `CSV_PATH` points to the demo current data:
`Data/Echem_Current/Anodic_Current_vs_Time.csv`
4. Ensure `CSV_OD_PATH` points to the demo optical data:
`Data/Optical_Delta-OD/Anodic_Optical_Density_vs_Time.csv`
5. Run the script:
```bash
python Unified Electro–Optical Transient Deconvolution Code.py

```



**Output:**
Three figures will be generated:

1. **Current Components:** Fitted current profiles including surface and bulk contributions with the fitted parameter table ($k_{Surf}$, $D_{Chem}$, $L_{Surf}$ and Biot number (Bi = $k_{Surf}$ $L$ / $D_{Chem}$)).
2. **Accumulated Charge:** Evolution of the fitted accumulated charge profiles during surface and bulk processes as well as the total charge storage component.
3. **ΔOD Components:** Optical density reconstructed from the linear combination fit of the accumulated charge profiles for surface and bulk storage processes.


## Method summary
| Step  | Time window                      | Purpose                              | Output                         | 
| ----- | -------------------------------- | ------------------------------------ | ------------------------------ |
| **1** | long-time region (e.g., 20–80 s) | linear fit of ln                     | τ<sub>d</sub>, D<sub>linear</sub> (initial D guess)                      |  
| **2** | short-time region (e.g., 0–60 s) | pre-fit surface+diffusion expression | improved initial guesses       | 
| **3** | short-time region (e.g., 0–60 s) | global fit with shared parameters    | **k, D, L<sub>Surf</sub>, Bi** |

ΔOD is then described by:

ΔOD(t) = ε<sub>Surf</sub> · ΔQ<sub>Surf</sub>(t) + ε<sub>Bulk</sub> · ΔQ<sub>Bulk</sub>(t)

where ΔQ components are time integrals derived from the fitted current components.

## Expected output (demo data)
| k (cm s-¹) | D (cm² s-¹) | L<sub>Surf</sub> (nm) | Bi  |
| ---------- | ----------- | --------------------- | --- |
| ~1E-5      | ~1E-11      | 5–25                  | ~40 |

---

## Using Your Own Data

To analyze your own experimental data, modify the **User-configurable constants** section at the top of `Unified Electro–Optical Transient Deconvolution Code.py` following these guidelines:

### 1. Data Preparation & Paths

* **File Format:**
* Electrochemical file: `time (s), I (mA)`
* UV-Vis file: `time (s), OD`


* **Data Requirements:**
* The input CSV must correspond to a **single applied potential step**.
* **Crucial:** The data must be trimmed to start exactly from the **current spike** (the onset of decay).
* *For Anodic steps:* Start from the highest positive peak current.
* *For Cathodic steps:* Start from the most negative peak current.




* **Update Paths:**
* Change `CSV_PATH` to your current transient file path.
* Change `CSV_OD_PATH` to your optical density file path.



### 2. Current-Only Analysis

If you **do not** have optical data, you must:

1. Comment out the `CSV_OD_PATH` variable.
2. Comment out the entire code block labeled `# 8 · Load UV–Vis data (0 – 60 s) & fit ΔOD(t)` (usually near the end of the script).
This ensures the script processes only the electrochemical data without raising errors.

### 3. Physical Parameters (`L_CM`)

* **Film Thickness:** Set `L_CM` to your sample's diffusion length in **cm** (e.g., for a 650 nm film, set `L_CM = 6.5e-5`).
* **Geometry Warning:** This model assumes a **planar thin-film geometry**. It is **not suitable** for describing ion intercalation kinetics in **spherical particles** (e.g., slurry electrodes).

### 4. Time Window Configuration

Adjust `T_SHORT_MAX` and `T_LIN_REGION` based on your experiment duration to optimize fitting.

* **Guideline (based on a 120s experiment):**
* `T_LIN_REGION = (20.0, 80.0)`: The long-time window used for the initial  vs.  linear fit to estimate the diffusion coefficient ().
* `T_SHORT_MAX = 60.0`: The short-time window (0–60s) used for the global non-linear fit to extract kinetic parameters (, ) and separate surface/bulk layers.


* **Scaling:** These values are not rigid. If your experiment is significantly longer or shorter than 120s, scale these windows proportionally. You may also fine-tune them to exclude noisy tails or initial instrument ringing.

---

## Method Summary

| Step | Time Window | Purpose | Output |
| --- | --- | --- | --- |
| **1** | Long-time (e.g., 20–80s) |  vs  linear fit | (Initial Guess) |
| **2** | Short-time (e.g., 0–60s) | Pre-fit surface + diffusion expression | Refined Initial Guesses |
| **3** | Short-time (e.g., 0–60s) | Global fit (common ) | **, , , ** |

For Optical Deconvolution:



*Where  components are integrated from the fitted current components.*

---

## Data Validity & Prerequisites
Before employing this code for kinetic fitting, users must possess a robust understanding of the electrochemical system under study. It is imperative to ensure—with high confidence—that the raw electrochemical current and optical absorbance data are **predominantly governed by electrochemical ion intercalation/deintercalation processes**.

Users must strictly exclude interference from parasitic side reactions to ensure the reliability of the quantitative results. Potential sources of error include, **but are not limited to**:

* Gas evolution (hydrogen or oxygen evolution) at high overpotentials.

* Anomalous current/optical signatures arising from electrode degradation.

* Sudden changes associated with phase transitions.

* Spontaneous chemical reactions (e.g., proton release).

**Recommendation**: Preliminary experimental validation is strongly advised. We recommend utilizing techniques such as in situ crystallographic or electronic structure characterization to phenomenologically analyze and rule out these interferences before performing the kinetic deconvolution.

---

## Citation

If you use this toolkit in your research, please cite the following software and paper:

```bibtex
@software{SurfBulk-Deconv-Toolkit_2025,
  author = {Luhan WEI},
  title = {Current–Transient & Optical Deconvolution Toolkit},
  year = {2025},
  url = {[https://github.com/LuhanWEI/SurfBulk-Deconv](https://github.com/LuhanWEI/SurfBulk-Deconv)},
  note = {MIT License}
}

@article{wei2024deconvoluting,
  title = {Deconvoluting Surface and Bulk Charge Storage Processes in Redox-Active Oxides by Integrating Electrochemical and Optical Insights},
  author = {Wei, Luhan and Hu, Yang and Huang, Yiwei and others},
  journal = {J. Am. Chem. Soc.},
  volume = {146},
  pages = {24167--24176},
  year = {2024},
  doi = {10.1021/jacs.4c09261}
}

```

---

## License

Released under the MIT License – see `LICENSE` for details.

---

## References
[1] Li, J.; Xiao, X.; Yang, F.; Verbrugge, M. W.; Cheng, Y.-T. Potentiostatic Intermittent Titration Technique for Electrodes Governed by Diffusion and Interfacial Reaction. J. Phys. Chem. C 2012, 116, 1472–1478. (https://doi.org/10.1021/jp207919q)

[2] Montella, C. Discussion of the Potential Step Method for the Determination of the Diffusion Coefficients of Guest Species in Host Materials. Part I. Influence of Charge Transfer Kinetics and Ohmic Potential Drop. J. Electroanal. Chem. 2002, 518, 61–83. (https://www.sciencedirect.com/science/article/abs/pii/S002207280100691X)

[3] Zhang, D.; Wang, R.; Wang, X.; Gogotsi, Y. In situ Monitoring Redox Processes in Energy Storage Using UV-Vis Spectroscopy. Nat. Energy 2023, 8, 567–576. (https://www.nature.com/articles/s41560-023-01240-9)

[4] Li, Y.; Chueh, W. C. Electrochemical and Chemical Insertion for Energy Transformation and Switching. Annu. Rev. Mater. Res. 2018, 48, 137–165. (https://doi.org/10.1146/annurev-matsci-070317-124525)
