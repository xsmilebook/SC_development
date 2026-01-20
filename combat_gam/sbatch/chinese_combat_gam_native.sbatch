#!/usr/bin/env bash
#SBATCH -J combat_chinese_gam_native
#SBATCH -p q_fat_c
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --output=outputs/logs/combat_gam/combat_chinese_gam_native_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge_new.rds}
OUTPUT_RDS=${2:-outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_native.rds}
TEST_N=${3:-0}

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment

# Avoid BLAS oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python combat_gam/scripts/run_combat_gam_neuroharmonize_native.py \
  --input-rds "${INPUT_RDS}" \
  --output-rds "${OUTPUT_RDS}" \
  --batch-col study \
  --id-col subID \
  --age-col Age \
  --sex-col Sex \
  --mean-fd-col mean_fd \
  --test-n "${TEST_N}"

