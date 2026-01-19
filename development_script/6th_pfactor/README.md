# 6th_pfactor
This folder contains scripts to evaluate the associations of a general psychopathological factor (*p*-factor) and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different *p*-factor levels. 

## S1st_pfactor_effect_continuous_ABCD.Rmd
This script aims to assess the relationship between *p*-factor and SC strength in the ABCD dataset. GAMMs were employed to control for potential confounders such as sex, head motion, and age smoothing. Subsequently, the GAMM analysis was conducted to investigate the significance of the association between *p*-factor and SC strength. In cases where the association proved significant for at least one edge, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. 

## S2nd_cbcl_totalraw_effect_continuous_ABCD.R
Sensitivity analysis using CBCL total raw score (`cbcl_scr_syn_totprob_r`) instead of *p*-factor in the ABCD dataset. In addition to saving model outputs (`outputs/results/cbcl_totprob/*.rds`), this script generates project-local figures:

- `outputs/figures/cbcl_totprob/Association/`: S-A rank scatter plots and 12×12 t-value matrices.
- `outputs/figures/cbcl_totprob/Interaction/`: high (90th) vs low (10th) CBCL trajectories across S-A deciles (1–10).

