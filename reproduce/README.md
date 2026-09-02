# reproduce/

Scripts that regenerate the cohort-level results reported in the paper. They are
drivers, not part of the estimator: nothing in `core/`, `ekf/`, `physiology/` or
`processor/` depends on this directory.

## Inputs

| input | in repo | note |
|---|---|---|
| `cfg_canonical.mat` | yes | configuration of the published run |
| `cohort_caseids.csv` | yes | the 209 analysed cases; column 1 is the patient index used throughout the code, column 2 the VitalDB CaseID |
| `patientDataFinal_auto.mat` | yes | 300-case VitalDB extract, ~67 MB, tracked at the repository root. Rebuildable from VitalDB using the CaseIDs in `cohort_caseids.csv`. |

`repro_paths.m` resolves all of these and reports by name anything that is
missing, rather than failing later inside a load.

Results are written to `reproduce/output/`.

## Scripts

| script | what it produces | approximate runtime |
|---|---|---|
| `run_cohort_montecarlo.m` | cohort-wide Monte Carlo equivalent-set fraction and L2 radius, replacing the single-patient numbers | minutes |
| `run_q_sweep_4d.m` | 4D MAE and parameter drift against process-noise scale `q` | ~6 x 15 min |
| `run_allmodels_sweep.m` | all five models against `q`, in-sample | ~n_q x 15 min |
| `run_holdout_sweep.m` | held-out split: adapt on the first half of each case, freeze, predict the second | ~n_q x 15 min |
| `run_bismin_sweep.m` | held-out accuracy against the fixed lower asymptote `BISmin` | ~8 x 15 min |
| `run_kgamma_experiment.m` | the 2D `(k, gamma)` model: in-sample and held-out accuracy, pairwise contrasts, parameter and bound diagnostics | ~n_q x 15 min |
| `run_fim_directions.m` | per-patient cumulative-FIM eigenvalues, leading eigenvectors and effective rank | minutes |
| `run_phase_split.m` | induction and maintenance accuracy, split at 10 minutes | ~15 min |
| `c50_ratio_test/run_c50_variants.m` | held-out effect of freezing the potency ratio `C50R/C50P` at its population value | ~4 x 15 min |

Each of the sweep drivers wraps a cohort function of the same name prefix
(`run_cohort_van`, `run_cohort_allmodels`, `run_cohort_holdout`), which loop the
209 cases and return per-patient metrics. The cohort functions take `cfg` and the
process-noise scale as arguments and hold no state.

`analysis/monte_carlo_cohort.m` stays in `analysis/` — it is an analysis routine,
not a driver.

## Conventions

- Patient index `i` is the row of `patientDataFinal_auto.mat`; `cohort_caseids.csv`
  maps it to the VitalDB CaseID.
- Cases shorter than 900 samples are excluded, as in `run_analysis.m`.
- A case that raises an error is reported with its index and message; it is not
  silently dropped.
- Every estimator sets `Q = cfg.q * P0` from its own initial covariance, so the
  single dimensionless scalar `cfg.q` means the same thing for all of them. The
  sweeps set `cfg.q` and nothing else; `param_scales.m` holds only the initial
  standard deviations, used to normalize parameter drift and roughness.
