# hygro-dual-damage-UMAT

Coupled hygro-viscoelastic dual-damage UMAT for 3D-printed CF/Onyx composites in Abaqus/Standard.

![Graphical Abstract](figures/graphical_abstract.png)

![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-teal.svg)
![Language: FORTRAN](https://img.shields.io/badge/Language-FORTRAN-blue)
![Solver: Abaqus/Standard](https://img.shields.io/badge/Solver-Abaqus%2FStandard-orange)
![Thesis: LUTPub](https://img.shields.io/badge/Thesis-URN%3ANBN%3Afi--fe2026060260742-green)

**Four-step research narrative:**
biphasic experimental degradation → Carter–Kibler dual-population transport →
dual-damage UMAT calibrated 0–90 days → bootstrap identifiability criterion

## Overview

This repository contains a FORTRAN UMAT implementation for modelling
moisture-driven degradation in fibre-reinforced polymer composites.
The model combines coupled moisture transport, reversible plasticisation,
irreversible hydrolytic damage, and viscoelastic response in a
finite-element framework.

The project also includes Python tools for sensitivity analysis and
parametric bootstrapping, plus publication-quality Matplotlib figures
used in calibration and uncertainty quantification.

> **Thesis (open access):** Irfan Irfan, *Coupled Hygro-Viscoelastic
> Dual-Damage UMAT for 3D-Printed CF/Onyx Composites*, LUT University,
> 2026. [URN:NBN:fi-fe2026060260742](https://urn.fi/URN:NBN:fi-fe2026060260742)

---

## Main Features

- FORTRAN UMAT implementation for Abaqus/Standard.
- Hygro-mechanical degradation modelling for composite materials.
- Python surrogate for sensitivity analysis and uncertainty quantification.
- Bootstrap-based parameter identifiability study.
- Publication-quality plots and processed data.

---

## Repository Structure

```text
UMAT/                FORTRAN source and Abaqus verification input files
Python_surrogate/    Python scripts for surrogate modelling and UQ
figures/             Matplotlib figures used in the study
data/                Experimental and processed CSV data
docs/                Theory notes and manuscript drafts
examples/            Simple Abaqus example and usage script
```

---

## UMAT Implementation — Key Algorithms

The UMAT executes three coupled procedures at every integration point
using a backward-Euler implicit scheme. Full derivations are in
Appendix A of the thesis (URN:NBN:fi-fe2026060260742).

### A.1 — Carter–Kibler Moisture Diffusion

The coupled mobile/bound moisture ODE system is solved by a 2×2 linear
system inverted with Cramer's rule. A fallback to explicit Euler is
triggered if the determinant falls below 10⁻²⁰.

```fortran
IF (dtimedays .GT. 0.D0) THEN
   CKA(1,1) = 1.D0 + dtimedays*(1.D0/taudiffdays + betatrap)
   CKA(1,2) = dtimedays*(1.D0/taudiffdays) - kappa
   CKA(2,1) = -dtimedays*betatrap
   CKA(2,2) = 1.D0 + dtimedays*kappa
   CKB(1)   = Cm + dtimedays*(Csat/taudiffdays)
   CKB(2)   = Cb
   DET = CKA(1,1)*CKA(2,2) - CKA(1,2)*CKA(2,1)
   IF (ABS(DET) .LT. 1.D-20) THEN
      Cm_new = Cm + dtimedays*((Csat-Cm-Cb)/taudiffdays
     &          - betatrap*Cm + kappa*Cb)
      Cb_new = Cb + dtimedays*(betatrap*Cm - kappa*Cb)
   ELSE
      Cm_new = (CKB(1)*CKA(2,2) - CKA(1,2)*CKB(2)) / DET
      Cb_new = (CKA(1,1)*CKB(2) - CKA(2,1)*CKB(1)) / DET
   END IF
   Cm = MAX(0.D0, MIN(Csat, Cm_new))
   Cb = MAX(0.D0, MIN(Csat-Cm, Cb_new))
END IF
Ctotal = MAX(0.D0, MIN(Csat, Cm + Cb))
```

> **State variable layout:** `STATEV(7)` → Cₘ, `STATEV(8)` → C_b

---

### A.2 — Dual-Damage Evolution Laws

Reversible plasticisation `P_i` (driven by mobile moisture Cₘ) and
irreversible hydrolytic damage `H_i` (driven by bound moisture C_b)
are updated separately. The hydrolytic ratchet enforces `H_i` to be
monotonically non-decreasing.

```fortran
! Reversible plasticisation driven by mobile moisture Cm
DO I=1,3
   IF (Cm .GT. 1.D-10) THEN
      P_target = 1.D0 - DEXP(-kP(I)*(Cm**mP))
   ELSE
      P_target = 0.D0
   END IF
   IF (dtime_days .GT. 0.D0 .AND. tauP_days .GT. 0.D0) THEN
      expm = DEXP(-dtime_days/tauP_days)
      P(I) = P_target + (P(I) - P_target)*expm
   ELSE
      P(I) = P_target
   END IF
   P(I) = MAX(0.D0, MIN(0.99D0, P(I)))
END DO

! Irreversible hydrolysis driven by bound moisture Cb
DO I=1,3
   IF (Cb .GT. 1.D-10 .AND. dtime_days .GT. 0.D0) THEN
      expm    = DEXP(-kH(I)*(Cb**mH)*dtime_days)
      H_target = 1.D0 - (1.D0 - H(I))*expm
      H(I) = MAX(H(I), H_target)   ! irreversible ratchet
   END IF
   H(I) = MAX(0.D0, MIN(0.99D0, H(I)))
END DO
```

> **State variable layout:** `STATEV(1:3)` → P_i, `STATEV(4:6)` → H_i

---

### A.3 — Damage-Coupled Generalised Maxwell Update

Damage accelerates viscoelastic relaxation by reducing the effective
relaxation time via τ_J = τ₀ · exp(−s_P·P̄ − s_H·H̄), and scales
branch moduli through ξ = 1 − ½(P̄ + H̄).

```fortran
Pbar = (P(1)+P(2)+P(3))/3.D0
Hbar = (H(1)+H(2)+H(3))/3.D0
damage_factor = MAX(0.1D0, MIN(1.D0, 1.D0-0.5D0*(Pbar+Hbar)))

DO K=1,NBR
   tau0s = PROPS(idx_tau0 + (K-1))
   tauJ  = MAX(1.D-12, tau0s*sec2day
  &          * DEXP(-SPfac*Pbar - SHfac*Hbar))
   alpha = MIN(1.0D12, dtime_days/(tauJ + 1.D-18))
   IF (alpha .LT. 1.D-6) THEN
      gamma = alpha - alpha*alpha
   ELSE
      gamma = alpha/(1.D0 + alpha)
   END IF
   DO J=1,6
      iv     = i0 + (K-1)*6 + J
      ev_new = (STATEV(iv) + alpha*epsm(J))/(1.D0 + alpha)
      IF (ABS(ev_new).GT.1.D6 .OR. ev_new.NE.ev_new)
  &      ev_new = 0.5D0*STATEV(iv)
      STATEV(iv) = ev_new
      delta(J)   = epsm(J) - ev_new
   END DO
   DO I=1,6
      DO J=1,6
         gJk = PROPS(idx_g + (K-1)*6 + (J-1))
         sig_vis(I)  = sig_vis(I)
  &                  + gJk*damage_factor*Cij(I,J)*delta(J)
         DDSDDE(I,J) = DDSDDE(I,J)
  &                  + gJk*damage_factor*gamma*Cij(I,J)
      END DO
   END DO
END DO
```

> **State variable layout:** `STATEV(9 : 9+6·N_br)` → viscous overstress strains

---

## Python Surrogate — Bootstrap UQ Workflow

From `Python_surrogate/bootstrap_uq.py`. The surrogate reproduces
UMAT prediction intervals without an Abaqus runtime.

```python
import numpy as np
from scipy.optimize import minimize

def ck_surrogate(t, params):
    """Carter-Kibler dual-population moisture saturation model."""
    C_inf, kA, kB = params
    C_bound  = C_inf * kB / (kA + kB)
    C_mobile = C_inf * kA / (kA + kB)
    return (C_mobile * (1 - np.exp(-kA * t))
          + C_bound  * (1 - np.exp(-kB * t)))

# Parametric bootstrap (n=1000) for parameter identifiability
rng = np.random.default_rng(42)
boot_params = np.zeros((1000, 3))
for i in range(1000):
    idx = rng.integers(0, len(t_data), len(t_data))
    res = minimize(residual, p0,
                   args=(t_data[idx], E_data[idx]),
                   method='Nelder-Mead')
    boot_params[i] = res.x

ci_lower, ci_upper = np.percentile(boot_params, [2.5, 97.5], axis=0)
```

**Key finding:** A model achieving <1.2 % RMSE over 90 days still
carries ±22 pp one-year prediction uncertainty — reducible only by
data beyond the bound-moisture trapping timescale (τ_trap ≈ 166 days).

---

## Results Gallery

| Calibration Fit | Bootstrap Prediction Interval |
|:-:|:-:|
| ![Calibration](figures/calibration_fit.png) | ![Bootstrap CI](figures/bootstrap_ci.png) |

| Dual Damage Evolution | OAT Sensitivity |
|:-:|:-:|
| ![Damage](figures/damage_evolution.png) | ![Sensitivity](figures/sensitivity_oat.png) |

---

## Files

- `UMAT/H2_coupled.for` — main UMAT source file.
- `Python_surrogate/surrogate.py` — surrogate model and post-processing.
- `Python_surrogate/bootstrap_uq.py` — parametric bootstrap workflow.
- `Python_surrogate/sensitivity_oat.py` — one-at-a-time sensitivity analysis.
- `figures/` — calibration plots, damage evolution, bootstrap intervals, sensitivity results.
- `data/` — raw and cleaned experimental datasets.

---

## How to Use

1. Compile the UMAT with Abaqus/Standard.
2. Link the subroutine in your input file.
3. Run the verification example in `examples/`.
4. Use the Python scripts to reproduce sensitivity and bootstrap analyses.

```text
Abaqus input file → UMAT simulation → Python surrogate → sensitivity analysis → bootstrap UQ → figures
```

---

## Citation

If you use this code, please cite the thesis:

> Irfan Irfan. *Coupled Hygro-Viscoelastic Dual-Damage UMAT for
> 3D-Printed CF/Onyx Composites*. Master's Thesis, LUT University,
> Finland, 2026.
> URN:NBN:fi-fe2026060260742
> https://urn.fi/URN:NBN:fi-fe2026060260742

A `CITATION.cff` is included in this repository for automated citation
by GitHub, Zenodo, and Zotero.

---

## License

This work is licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to share and adapt this material for any purpose,
including commercially, as long as you give appropriate credit
(cite the thesis above) and indicate any changes made.

© 2026 Irfan Irfan, LUT University.
