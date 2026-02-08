# 6th_pfactor
This folder contains scripts to evaluate the associations of a general psychopathological factor (*p*-factor) and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different *p*-factor levels. 

## S1st_pfactor_effect_continuous_ABCD.Rmd
This script aims to assess the relationship between *p*-factor and SC strength in the ABCD dataset. GAMMs were employed to control for potential confounders such as sex, head motion, and age smoothing. Subsequently, the GAMM analysis was conducted to investigate the significance of the association between *p*-factor and SC strength. In cases where the association proved significant for at least one edge, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. 

## S2nd_cbcl_totalraw_effect_continuous_ABCD.R
Sensitivity analysis using CBCL total raw score (`cbcl_scr_syn_totprob_r`) instead of *p*-factor in the ABCD dataset. In addition to saving model outputs (`outputs/results/cbcl_totprob/*.rds`), this script generates project-local figures:

- `outputs/figures/cbcl_totprob/Association/`: S-A rank scatter plots and 12×12 t-value matrices.
- `outputs/figures/cbcl_totprob/Interaction/`: high (90th) vs low (10th) CBCL trajectories across S-A deciles (1–10).

## run_abcd_lmm_agewp_agebp_baselineage_pfactor_interaction_decileavg.R
ABCD baseline-age decomposition LMM interaction workflow for *p*-factor (no t-value matrix output).  
- `age_bp`: baseline age (from baseline event; fallback to minimum age if baseline tag is missing).  
- `age_wp`: current age minus baseline age.  
- *p*-factor uses baseline-only `GENERAL_base` (aligned with cognition baseline strategy).  
- Two interaction models are fitted for each edge:
  - `y ~ age_wp * cov + age_bp + sex + mean_fd + (1 | subID)`
  - `y ~ age_wp + age_bp * cov + sex + mean_fd + (1 | subID)`
- Produces decile plots (1-10) for both `age_wp` and `age_bp` under low/high *p*-factor in two pipelines:
  - edge-first: fit per edge, then average predictions within decile;
  - decile-avg-SC-first: aggregate SC ratio within each decile first, then fit interaction models.
- Decile plotting uses the actual `age_wp` (for `age_wp*cov`) and `age_bp` (for `age_bp*cov`) values on x-axis; axis labels remain `Age` and `SC strength (ratio)` to match previous figures.
- In `decile-avg-SC-first`, interaction significance is tested for both `age_wp:cov` and `age_bp:cov` by LRT against
  `y ~ age_wp + cov + age_bp + sex + mean_fd + (1 | subID)`, followed by FDR correction (`p_agewp_cov_fdr`, `p_agebp_cov_fdr`); significance is not drawn on figures.

