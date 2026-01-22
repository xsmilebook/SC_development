#!/usr/bin/env bash
#SBATCH -J plot_hcpd_ch_r2
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --output=outputs/logs/combat_gam/plot_hcpd_chinese_variance_decomposition_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
# NOTE: `scdevelopment_r41` may fail on some compute nodes due to GLIBC mismatch
# (e.g., `GLIBC_2.32 not found` when loading `cli.so`). Default to the older
# `scdevelopment` env that previously worked on compute nodes; allow override.
conda activate "${CONDA_ENV:-scdevelopment}"

# Avoid user-level R libraries overriding conda packages (ABI mismatch).
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Avoid BLAS oversubscription when using mclapply.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript --vanilla combat_gam/scripts/plot_hcpd_chinese_variance_decomposition.R
