# hygro-dual-damage-UMAT

Coupled hygro-viscoelastic dual-damage UMAT for 3D-printed CF/Onyx composites in Abaqus/Standard.

![Graphical Abstract](figures/graphical_abstract.jpg)

 **Four-step research narrative:** biphasic experimental degradation → Carter–Kibler dual-population transport → dual-damage UMAT calibrated 0–90 days → bootstrap identifiability criterion


## Overview
This repository contains a FORTRAN UMAT implementation for modelling moisture-driven degradation in fibre-reinforced polymer composites. The model combines coupled moisture transport, reversible plasticisation, irreversible hydrolytic damage, and viscoelastic response in a finite-element framework.

The project also includes Python tools for sensitivity analysis and parametric bootstrapping, plus publication-quality Matplotlib figures used in calibration and uncertainty quantification.

## Main features
- FORTRAN UMAT implementation for Abaqus/Standard.
- Hygro-mechanical degradation modelling for composite materials.
- Python surrogate for sensitivity analysis and uncertainty quantification.
- Bootstrap-based parameter identifiability study.
- Publication-quality plots and processed data.

## Repository structure
```text
UMAT/                FORTRAN source and Abaqus verification input files
Python_surrogate/    Python scripts for surrogate modelling and UQ
figures/             Matplotlib figures used in the study
data/                Experimental and processed CSV data
docs/                Theory notes and manuscript drafts
examples/            Simple Abaqus example and usage script
```

## Files
- `UMAT/H2_coupled.for` — main UMAT source file.
- `Python_surrogate/surrogate.py` — surrogate model and post-processing.
- `Python_surrogate/bootstrap_uq.py` — parametric bootstrap workflow.
- `Python_surrogate/sensitivity_oat.py` — one-at-a-time sensitivity analysis.
- `figures/` — calibration plots, damage evolution, bootstrap intervals, and sensitivity results.
- `data/` — raw and cleaned experimental datasets.

## How to use
1. Compile the UMAT with Abaqus/Standard.
2. Link the subroutine in your input file.
3. Run the verification example in `examples/`.
4. Use the Python scripts to reproduce sensitivity and bootstrap analyses.

## Example workflow
```text
Abaqus input file -> UMAT simulation -> Python surrogate -> sensitivity analysis -> bootstrap UQ -> figures
```

## Figures
Add your main plots here when ready:

- damage evolution
- calibration fit
- bootstrap prediction interval
- sensitivity bar chart

## Citation
If you use this code, please cite the corresponding thesis or manuscript.

## License
Add your preferred license here before public release.
