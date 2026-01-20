#!/usr/bin/env bash
#SBATCH -J plot_abcd_r2
#SBATCH -p q_fat_c
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --output=outputs/logs/combat_gam/plot_abcd_variance_decomposition_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment_r41

# Avoid user-level R libraries overriding conda packages (ABI mismatch).
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Avoid BLAS oversubscription when using mclapply.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript --vanilla combat_gam/scripts/plot_abcd_variance_decomposition.R

