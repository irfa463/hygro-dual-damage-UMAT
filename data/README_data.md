# Experimental Data

This folder contains tensile-characterisation data for 3D-printed continuous
carbon-fibre-reinforced Onyx (CF/Onyx) specimens conditioned in distilled water
at \(26\,^\circ\mathrm{C}\) for up to 90 days. The data support the constitutive
modelling, uncertainty quantification, and validation discussed in the associated
journal manuscript.

## Conditioning Protocol

| Parameter | Value |
|---|---|
| Material | 3D-printed continuous CF/Onyx composite |
| Printer | Markforged Mark Two |
| Matrix | Onyx nylon-based polymer |
| Fibre architecture | Continuous carbon fibre aligned with the tensile direction; Onyx layers deposited at \(\pm45^\circ\) |
| Immersion medium | Distilled water |
| Temperature | \(26\,^\circ\mathrm{C}\) |
| Conditioning periods | 0, 15, 30, 60, and 90 days |
| Tensile standard | ASTM D3039 |
| Test rate | \(2\,\mathrm{mm\,min^{-1}}\) |
| Nominal specimen dimensions | \(250\times25\times2.2\,\mathrm{mm}\) |

## Primary Tensile Dataset

The primary ageing dataset contains tensile-characterisation results for specimens
tested after 0, 15, 30, 60, and 90 days of distilled-water immersion.

| File | Conditioning | Description |
|---|---:|---|
| `tensile_T0.csv` | 0 days | As-printed dry baseline |
| `tensile_T15.csv` | 15 days | Water-aged tensile data |
| `tensile_T30.csv` | 30 days | Water-aged tensile data |
| `tensile_T60.csv` | 60 days | Water-aged tensile data |
| `tensile_T90.csv` | 90 days | Water-aged tensile data |

## Short-Term Supplementary Campaign

The supplementary workbook contains the short-term gravimetric and tensile campaign
performed at 15 days to provide an independent check of early-time water uptake and
to assess the recoverability of the mechanical response after drying.

| File | Description |
|---|---|
| `short_term_experimental_compaign.xlsx` | Raw Instron/DIC tensile data and gravimetric mass measurements for unaged, 15-day immersed, and 15-day immersed–dried specimens |

The workbook contains the following worksheets:

| Worksheet | Description |
|---|---|
| `load-displacament Zero days` | Load, displacement, Instron, and DIC-related data for the unaged specimens |
| `stress-strain instron OBS` | Stress–strain data obtained from the Instron measurements |
| `stress-strain DIC GOM DATA` | Stress–strain data obtained using DIC/GOM strain measurements and Instron load data |
| `MASS` | Gravimetric measurements before immersion, during the 15-day immersion period, after drying, and on the tensile-testing day |

The short-term campaign includes three nominal conditions:

1. Unaged reference specimens.
2. Specimens immersed in distilled water for 15 days and tested without oven drying.
3. Specimens immersed for 15 days and subsequently oven-dried before tensile testing.

The gravimetric measurements indicate an apparent mass increase of approximately
\(5.7\text{--}5.9\,\mathrm{wt\%}\) after 15 days of immersion. Because the post-drying
masses were not consistently comparable with the initial masses, the short-term
campaign is not used to calculate a quantitative residual bound-moisture fraction or
to identify the mobile/bound moisture partition.

The short-term tensile results are provided as supplementary observations. They were
not pooled with the primary 0–90-day ageing dataset because the supplementary
campaign did not reproduce the strength reduction measured in the primary 15-day
campaign. The data are therefore intended to support discussion of early-time
moisture uptake, moisture-state control, and experimental limitations rather than
to provide an independent calibration of the irreversible degradation parameters.

## Data Columns

The primary CSV files contain the following columns:

| Column | Description | Unit |
|---|---|---|
| `specimen_id` | Specimen label | — |
| `gauge_length_mm` | Gauge length | mm |
| `width_mm` | Specimen width | mm |
| `thickness_mm` | Specimen thickness | mm |
| `cross_section_mm2` | Cross-sectional area | mm² |
| `max_load_N` | Peak load before failure | N |
| `UTS_MPa` | Ultimate tensile strength | MPa |
| `strain_at_break_pct` | Strain at failure | % |
| `elastic_modulus_GPa` | Longitudinal elastic modulus | GPa |

The Excel workbook contains raw and processed stress–strain, load–displacement,
strain, and mass measurements. Column names and worksheet structures may differ
from the standardized CSV files.

## Data-Processing Notes

- The primary tensile properties were calculated from the maximum valid stress on
  the pre-failure loading branch.
- Isolated records associated with post-failure motion, strain resets, or
  nonphysical terminal values should not be interpreted as additional tensile
  failure points.
- The short-term workbook is retained in its original form to preserve the raw
  experimental record.
- Surface moisture was removed before tensile testing according to the available
  laboratory procedure.
- The time between water removal, weighing, drying, and tensile testing was not
  identical for all short-term specimens; this should be considered when interpreting
  their moisture state.
- The short-term mass measurements provide an apparent total uptake at 15 days but
  do not independently identify mobile and bound moisture populations.
- All specimens were manufactured using the same nominal Markforged Mark Two
  printing route unless otherwise stated in the associated manuscript.

## Reproducibility

The CSV files provide the processed specimen-level tensile properties used for
modelling and statistical analysis. The Excel workbook provides the underlying
short-term Instron, DIC/GOM, stress–strain, and gravimetric records.

The data should be interpreted together with the manuscript, which describes the
conditioning protocol, constitutive model, uncertainty quantification, data
screening, and limitations of the short-term supplementary campaign.

