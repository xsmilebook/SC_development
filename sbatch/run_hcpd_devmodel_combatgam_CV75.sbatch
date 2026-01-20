#!/bin/bash
#SBATCH -J hcpd_devmodel_cg_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH --mem=120G
#SBATCH -t 48:00:00
#SBATCH -o outputs/logs/hcpd_devmodel_combatgam_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

CVTHR="${CVTHR:-75}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"
EUCLID_CSV="${EUCLID_CSV:-${project_root}/wd/interdataFolder_HCPD/average_EuclideanDistance_12.csv}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment_r41

# Avoid user-level R libraries overriding conda packages (ABI mismatch).
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS=""
export R_LIBS_SITE="${CONDA_PREFIX}/lib/R/library"
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Avoid thread oversubscription (mclapply already uses cores).
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[INFO] project_root=${project_root}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] EUCLID_CSV=${EUCLID_CSV}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

echo "[INFO] Checking R package dependencies (system libraries masked)"
Rscript --vanilla -e 'pkgs<-c("mgcv","nlme","tidyverse","gratia","mvnfast","psych","R.matlab","patchwork","RColorBrewer","scales","reshape","ecostats"); cat(".libPaths():\n"); print(.libPaths()); miss<-pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(miss)>0) stop(paste0("Missing packages: ", paste(miss, collapse=", "))); for(p in pkgs) cat(p,":",as.character(packageVersion(p)),"\n")'

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  --input_rds="${INPUT_RDS}"

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}"

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}"

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  --input_rds="${INPUT_RDS}" \
  --euclid_csv="${EUCLID_CSV}"
