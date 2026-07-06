# Gut Microbiome-Diet-Type 2 Diabetes Computational Model

Computational model simulating the effect of diet on glycemic control in Type 2
Diabetes (T2D) via a gut microbiome / short-chain fatty acid (SCFA) mechanistic
pathway, calibrated using literature-derived constants and an original survey.

## Overview

This repository contains the simulation code, model parameters, and aggregate
outputs supporting the paper "Modeling the Effects of Dietary Patterns on Gut Microbiome-Derived Short-Chain Fatty Acids and Type 2 Diabetes Mellitus" (currently under
review). The model simulates four dietary patterns (Western, Mediterranean,
Plant-based, Ketogenic) across a bootstrapped virtual patient population, and
outputs are validated against published meta-analyses of dietary interventions
in T2D.

## Repository Contents

| File | Description |
|---|---|
| `simulate.py` | Main simulation: 4-layer model (meal-level glucose/insulin, gut fermentation, chronic adaptation, HbA1c kinetics), Monte Carlo over 4 diets |
| `plot_results.py` | Generates glucose/HbA1c comparison charts from simulation output |
| `compare_to_literature.py` | One-sample z-test comparing simulated vs. published HbA1c effect sizes |
| `parameters.csv` | Model rate constants and coefficients (literature-derived and survey-calibrated, see paper Methods) |
| `diet_targets.csv` | Per-diet SCFA/weight steady-state targets |
| `diet_meals.csv` | Per-diet meal carbohydrate/fiber composition |
| `survey_covariates_SYNTHETIC_EXAMPLE.csv` | Fabricated example covariates file, for verifying the code runs (NOT real data) |
| `survey_instrument.md` | Survey questions administered (see Data Availability) |
| `example_outputs/` | Sample simulation outputs (trajectories, validation summaries, sample intraday traces, charts) from a full verified run |

## Requirements

Python 3.x, standard library only -- no external packages required to run
`simulate.py` or `compare_to_literature.py`. `plot_results.py` requires
`matplotlib`.

## Usage

**Note:** `simulate.py` requires a `survey_covariates.csv` file in the same
directory, which is not included in this repository (see Data Availability
below). A synthetic example with fabricated data,
`survey_covariates_SYNTHETIC_EXAMPLE.csv`, is provided so the code can be
verified to run end-to-end; rename it to `survey_covariates.csv` to test. Do
not interpret output generated from the synthetic file as representing real
results -- it exists only to confirm the code runs correctly.

Place `simulate.py`, `parameters.csv`, `diet_targets.csv`, `diet_meals.csv`,
and a valid `survey_covariates.csv` in the same directory, then:

```bash
python simulate.py
```

This produces, for both a 120-day and 365-day horizon:
- `mean_trajectories_by_condition_120d.csv` / `_365d.csv`
- `validation_summary_120d.csv` / `_365d.csv`
- `sample_patient_intraday_120d.csv` / `_365d.csv`

Runtime is approximately 2-3 minutes.

To generate comparison charts:
```bash
python plot_results.py mean_trajectories_by_condition_120d.csv T2D
```

To compare simulated results against published literature:
```bash
python compare_to_literature.py mean_trajectories_by_condition_120d.csv T2D
```

## Model Validation

Condition-group labels (T2D, Hypertension, Atherosclerosis, None, Other) are
used only for post hoc grouping of simulation output and are never supplied to
the model's equations. See the paper's Methods and Results sections for full
validation methodology.

## Data Availability

Individual-level survey covariates (`survey_covariates.csv`) are not included
in this repository due to participant confidentiality commitments made under
IRB-approved informed consent. Aggregated, de-identified summary statistics are
available upon request. The survey instrument (questions only, no responses)
is included as `survey_instrument.md`.



## Citation

Manuscript currently under review. Citation will be added upon publication.
