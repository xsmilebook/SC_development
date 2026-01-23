#!/usr/bin/env bash
#SBATCH -J combat_abcd_fc_base_gam
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --output=outputs/logs/combat_gam/combat_abcd_gam_fluidcomp_fc_baseline_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds}
OUTPUT_RDS=${2:-outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_fluidcomp_fc_baseline.rds}
TEST_N=${3:-0}

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment

# Avoid user-level R libraries overriding conda packages (ABI mismatch).
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
export R_LIBS=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Avoid BLAS oversubscription when using mclapply.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla combat_gam/scripts/run_abcd_nonlinear_combat_gam_fluidcomp_fc_baseline.R \
  "${INPUT_RDS}" \
  "${OUTPUT_RDS}" \
  "${TEST_N}"

