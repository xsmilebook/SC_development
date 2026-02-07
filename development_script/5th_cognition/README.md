# 5th_cognition
This folder contains scripts to evaluate the associations between fluid cognition and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different cognition levels. 

Additional ABCD variants are provided as standalone Rscript entry points (e.g., `run_abcd_cognition_fluid_uncorrected_S{1,2,3}.R`, `run_abcd_cognition_comp_agecorrected_S{1,2,3}.R`, and `run_abcd_cognition_fluid_fc_S{1,2,3}.R`).


## S1st_testCog_correlation_wholesample_*.Rmd
This script aims to assess the relationship between fluid cognition and SC strength in the HCP-D or ABCD dataset. GAMs were employed to control for potential confounders such as sex, head motion, and age smoothing. Subsequently, the GAM analysis was conducted to investigate the significance of the association between fluid cognition and SC strength. In cases where the strength of at least one edge has a significant association with cognition, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. 


## S2nd_compositescorePlot_scatterplot.Rmd
This script generated the schema of fluid cognition components and the scatter plots between SC strength and fluid cognition for 3 exemplified connections. Given that the associations were only significant in the ABCD dataset, the scatter plots were only depicted using the ABCD dataset.

## S3rd_SCdev_vary_by_cognition.Rmd
This script depicts developmental trajectories varying by baseline cognition levels in the ABCD dataset. We fitted an age-by-cognition interaction GAM model for each connection. Using the acquired models, we estimated SC strength by assigning cognitive performance as low and high levels, respectively. To define these levels, we used the 10th percentile of baseline cognitive performance for the low level and the 90th percentile for the high level. We then averaged trajectories for low and high cognition levels independently within deciles of the S-A connectional axis for visualization purposes.

## run_abcd_lmm_agewp_agebp_baselineage_cognition_interaction_decileavg.R
ABCD baseline-age decomposition LMM interaction workflow for cognition (no t-value matrix output).  
- `age_bp`: baseline age (from baseline event; fallback to minimum age if baseline tag is missing).  
- `age_wp`: current age minus baseline age.  
- Cognition uses baseline-only `nihtbx_fluidcomp_uncorrected_base`.  
- Two interaction models are fitted for each edge:
  - `y ~ age_wp * cov + age_bp + sex + mean_fd + (1 | subID)`
  - `y ~ age_wp + age_bp * cov + sex + mean_fd + (1 | subID)`
- Produces decile plots (1-10) for both `age_wp` and `age_bp` under low/high cognition in two pipelines:
  - edge-first: fit per edge, then average predictions within decile;
  - decile-avg-SC-first: aggregate SC ratio within each decile first, then fit interaction models.
- Decile plotting uses the actual `age_wp` (for `age_wp*cov`) and `age_bp` (for `age_bp*cov`) values on x-axis; axis labels remain `Age` and `SC strength (ratio)` to match previous figures.
- In `decile-avg-SC-first`, interaction significance is tested for both `age_wp:cov` and `age_bp:cov` by LRT against
  `y ~ age_wp + cov + age_bp + sex + mean_fd + (1 | subID)`, followed by FDR correction (`p_agewp_cov_fdr`, `p_agebp_cov_fdr`) and figure annotation.


