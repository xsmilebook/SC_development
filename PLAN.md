# Plan

We will implement the reviewer-required ComBat-GAM pipeline with a dry-run on a small subset first, then full `sbatch` runs on `q_fat_c`, with clear logging and documentation updates per AGENTS. This plan covers HCP-D, ABCD (three covariate sets; cognition baseline-only), and Chinese Cohort (batch=`study`) using the `combat_gam/` code.

## Scope
- In: ComBat-GAM/Nonlinear-ComBat-GAM implementation, SLURM submission (q_fat_c), small-sample test before full runs, log directory setup, documentation updates.
- Out: Changes to preprocessing, SC construction, or non-ComBat analyses.

## Action items
[ ] Define a project log folder (e.g., `outputs/logs/combat_gam/`) and a consistent, searchable filename pattern (e.g., `combat_<dataset>_<variant>_<date>_<jobid>.log`).
[ ] Create a minimal subset test plan per dataset (e.g., fixed N or specific sites) and document the selection criteria.
[ ] Implement HCP-D ComBat-GAM using `combat_gam/neuroHarmonize` with GAMM settings aligned to `gamfunction/gammsmooth.R` and batch=`site`.
[ ] Implement Chinese Cohort ComBat-GAM with batch=`study`, same GAMM settings.
[ ] Implement ABCD Nonlinear-ComBat-GAM for three covariate sets (A/B/C), with cognition restricted to baseline only.
[ ] Add SLURM `sbatch` scripts for q_fat_c: one for small test, one for full run; ensure logs go to the new log folder.
[ ] Validate outputs after test and full runs (dimensions, covariate presence, missingness, site balance) and record anomalies.
[ ] Update `docs/workflow.md` and a results note to document the new ComBat runs and outputs; then update `PROGRESS.md` and session log.

## Open questions
- None.
