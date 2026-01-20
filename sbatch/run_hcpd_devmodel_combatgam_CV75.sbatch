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
N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"
MAKE_MATRIX_GRAPHS="${MAKE_MATRIX_GRAPHS:-0}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment_r41

# Avoid system/user-level R libraries overriding our controlled libraries.
# IMPORTANT: some compute nodes have older GLIBC than the login node; packages
# compiled on the login node (even inside conda env) can fail with
# `GLIBC_2.xx not found`. To make the run robust, we install/override required
# R packages into a project-local library (compiled on the compute node).
project_r_lib="${project_root}/outputs/r_libs/scdevelopment_r41"
mkdir -p "${project_r_lib}"

export R_LIBS_USER="${project_r_lib}"
export R_LIBS=""
export R_LIBS_SITE="${CONDA_PREFIX}/lib/R/library"
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null
export R_PROFILE=/dev/null
export R_ENVIRON=/dev/null

# Avoid thread oversubscription (mclapply already uses cores).
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

echo "[INFO] project_root=${project_root}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] EUCLID_CSV=${EUCLID_CSV}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"
echo "[INFO] CONDA_PREFIX=${CONDA_PREFIX}"
echo "[INFO] project_r_lib=${project_r_lib}"
echo "[INFO] Rscript=$(command -v Rscript || true)"
Rscript --version || true

echo "[INFO] Ensuring R package dependencies (system/user libraries masked; project-local installs allowed)"
Rscript --vanilla -e 'lib_user<-Sys.getenv("R_LIBS_USER"); lib_site<-Sys.getenv("R_LIBS_SITE"); dir.create(lib_user, showWarnings=FALSE, recursive=TRUE); .libPaths(unique(c(lib_user, lib_site))); cat(".libPaths():\n"); print(.libPaths()); options(repos=c(CRAN="https://cloud.r-project.org")); pkgs<-c("mgcv","nlme","tidyverse","gratia","mvnfast","psych","R.matlab","patchwork","RColorBrewer","scales","reshape","ecostats"); for(p in pkgs){cat("== ",p," ==\n", sep=""); ok<-tryCatch({suppressPackageStartupMessages(library(p, character.only=TRUE)); TRUE}, error=function(e){cat("LOAD_ERROR: ",conditionMessage(e),"\n", sep=""); FALSE}); if(!ok){cat("Installing into ",lib_user,": ",p,"\n", sep=""); install.packages(p, lib=lib_user, dependencies=TRUE); ok2<-tryCatch({suppressPackageStartupMessages(library(p, character.only=TRUE)); TRUE}, error=function(e){cat("RELOAD_ERROR: ",conditionMessage(e),"\n", sep=""); FALSE}); if(!ok2) stop("Failed to load after install: ", p)}; cat("VERSION: ", as.character(packageVersion(p)),"\n", sep="") }'

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  --input_rds="${INPUT_RDS}" \
  ${N_EDGES:+--n_edges=${N_EDGES}}

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --skip_posterior="${SKIP_POSTERIOR}"

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}"

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  --input_rds="${INPUT_RDS}" \
  --euclid_csv="${EUCLID_CSV}" \
  --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}"
