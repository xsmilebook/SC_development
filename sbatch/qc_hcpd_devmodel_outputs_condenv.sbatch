#!/bin/bash
#SBATCH -J qc_hcpd_gamout
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 8
#SBATCH -o outputs/logs/qc_hcpd_gamout_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
CVTHR="${CVTHR:-75}"
FORCE="${FORCE:-0}"

# Conda env (adjust if needed)
CONDA_ROOT="${CONDA_ROOT:-/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-scdevelopment_r41}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

# Activate conda
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

echo "[INFO] CONDA_PREFIX=${CONDA_PREFIX}"
Rscript --version

# Shield user/system R libs: rely only on conda site library.
export R_LIBS_USER=""
export R_LIBS=""
export R_LIBS_SITE="${CONDA_PREFIX}/lib/R/library"
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

Rscript --vanilla development_script/2nd_fitdevelopmentalmodel/QC_check_gam_outputs_HCPD.R \
  --project_root="${project_root}" \
  --cvthr="${CVTHR}" \
  --force="${FORCE}"

