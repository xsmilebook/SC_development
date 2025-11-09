# 1st_dataclean
This folder contains codes to `extract the SC strength of large-scale SC matrices`, `conduct ComBat to harmonize SC strength from multiple acquisition sites`, `visualize the average SC strength matrix`, and `generate the age distribution plots for the final samples`. Separate scripts are provided for the ABCD and HCP-D datasets for the first three steps. For sensitivity analyses, users can specify the resolutions of large-scale SC matrices, and the additional covariate at the beginning of the scripts for the second and third steps.

## S1st_mergedata_*.R
**This script generated the index of spurious edges in the fine-grained SC network based on the thresholds of the 75th percentile (P75th) coefficient of variation (CV).** First, a dataframe is created where each column represents the streamline counts scaled by the volume of nodes in an edge. The normalized streamline counts are extracted from the fine-grained SC network built on the Schaefer-400 atlas. Edges connecting regions in the limbic system are excluded, resulting in 70,786 edges for the 376\*376 network. Second, CVs are computed for each edge. We apply consistency-based thresholds on the structural connectivity matrices to minimize the influence of false-positive connections. We select the thresholds of the P75th based on a previous study[(Baum, et. al., 2020)](https://www.pnas.org/doi/pdf/10.1073/pnas.1912034117), namely removing streamlines from the top quartile or three quartiles of inconsistent connections. The primary results are based on the structural connectivity matrices thresholded at the P75th CV.

## S2nd_mergedata_SA_ds_sumSC_*.R
This script extracts the SC strength of large-scale SC matrices. A dataframe will be generated in which each column represents the streamline counts scaled by node volumes of a connection in a large-scale structural connectivity matrix (e.g., 78 columns for 12\*12 matrix). Before running, ds.resolution should be determined. ds.resolution represents the resolution of the large-scale structural connectivity matrix. In this study, we computed the primary results based on the ds.resolution of 12. In addition, we also replicated our primary results based on the ds.resolution of 7 and 17, which are common resolutions used in previous studies on large-scale networks. The streamlines from the top quartile (P75th CV threshold) of inconsistent connections in the 376\*376 matrix were removed. The remaining streamline counts were summed up and divided by the average node volumes they connected.

## S3rd_combat_controlsite_*.R
We conducted ComBat to harmonize the structural connectivity strength from multiple acquisition sites, referring to a code improved by Richard Beare [code](https://github.com/PennLINC/Larsen_IronDevelopment/blob/master/combat.R). ComBat was performed separately for observations included in developmental models (including covariates for age, sex, and head motion), cognitive models (including covariates for fluid cognition, age, sex, and head motion), and p-factor models (including covariates for p-factor, age, sex, and head motion).

## S4th_plotSCdata_SA12_separateage_HCPD.R
This script plots the average structural connectivity matrices of 12\*12 at specific ages.

## S5th_demodescrip_plot.R
This script generates the age distribution plots for the final samples. In addition, we used this script to describe demographic and behavioral informations.

# merge_demography_info_and_screen
This folder contains codes for organizing demographic, cognitive, and psychopathologic variables and screening data for the HCP-D and ABCD datasets.
