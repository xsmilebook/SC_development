#!/bin/bash
#SBATCH -J cbcl_assoc_full
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH --mem=120G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/cbcl_assoc_full_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment_r41

# Avoid user-level R libraries overriding conda packages (ABI mismatch).
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

Rscript --vanilla development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD.R
