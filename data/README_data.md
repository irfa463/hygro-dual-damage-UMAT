# Experimental Data

This folder contains tensile characterisation data for 3D-printed CF/Onyx
specimens conditioned in distilled water at 26 °C for up to 90 days.
Data supports model calibration and validation in the associated journal paper.

## Conditioning Protocol

| Parameter | Value |
|-----------|-------|
| Material | 3D-printed CF/Onyx (Markforged) |
| Immersion medium | Distilled water |
| Temperature | 26 °C |
| Conditioning periods | 0, 15, 30, 60, 90 days |
| Standard | (add your standard, e.g. ASTM D3039) |

## Files

| File | Conditioning (days) | Description |
|------|-------------------|-------------|
| `tensile_T0.csv` | 0 (dry baseline) | As-printed, no immersion |
| `tensile_T15.csv` | 15 | Post-immersion tensile data |
| `tensile_T30.csv` | 30 | Post-immersion tensile data |
| `tensile_T60.csv` | 60 | Post-immersion tensile data |
| `tensile_T90.csv` | 90 | Post-immersion tensile data |

## Columns (all files)

| Column | Description | Unit |
|--------|-------------|------|
| `specimen_id` | Specimen label | — |
| `gauge_length_mm` | Gauge length | mm |
| `width_mm` | Specimen width | mm |
| `thickness_mm` | Specimen thickness | mm |
| `cross_section_mm2` | Cross-sectional area | mm² |
| `max_load_N` | Peak load at failure | N |
| `UTS_MPa` | Ultimate tensile strength | MPa |
| `strain_at_break_pct` | Strain at failure | % |
| `elastic_modulus_GPa` | Young's modulus | GPa |

## Notes

- All specimens printed with Markforged Mark Two printer
- Fibre layup: (add your layup, e.g. [0/90]s or isotropic Onyx)
- Specimens dried and weighed before immersion (T0 baseline)
- After conditioning, surface-dried with tissue before testing
- Testing conducted at LUT University, 2025
