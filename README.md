# Online BIS Estimation: Dimensionality Analysis

Source code accompanying the journal paper.

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox

## Quick Start

```matlab
cd('Source code')
run('run_analysis')
```

The script processes all patients in `patientDataFinal_auto.mat` and writes all outputs (figures, tables, and `.mat` results) to a `results/` directory.

## Data

`patientDataFinal_auto.mat` contains de-identified BIS, propofol infusion rate, remifentanil infusion rate, and demographic data from the VitalDB open dataset.

## Repository Structure

| Folder | Description |
|---|---|
| `core/` | BIS prediction math — Hill equation, Vanluchene/Greco models, C50 scaling |
| `ekf/` | Extended Kalman Filter estimators — 1D (k-scale), 2D (kP/kR), 3D (log-linear), 4D (RLS) |
| `physiology/` | Pharmacokinetic state propagation, effect-site dynamics, E0/BISmin estimation |
| `processor/` | Online orchestration — sample-by-sample processing, artifact gating, initialization |
| `data/` | Data loading from .mat file, demographics, results struct construction |
| `analysis/` | Statistical tests, Monte Carlo equivalence, reporting utilities |
| `figures/` | Journal figure and table generators (Figures 1–7, Tables 1–4) |
| `diagnostic_plots/` | Exploratory and debugging plots for development |
| `results/` | Generated at runtime — figures, tables, and `.mat` results |

## Algorithm Overview

Five model variants of increasing dimensionality are compared:

- **0D** — Population (no personalization)
- **1D** — Single k-scale EKF
- **2D** — Separate propofol/remifentanil scaling (kP, kR) via EKF with FIM-based identifiability gating
- **3D** — Log-linear EKF
- **4D** — Full Vanluchene/Greco parameter estimation

## Citation

*Citation details will be added upon publication.*
