#!/bin/bash
# ABCD | Nonlinear ComBat-GAM | SA7 | CV75
#
# Input (raw merge RDS; absolute path on GPFS):
# - /ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.merge.rds
#
# Output (written under SCDevelopment):
# - outputs/results/combat_gam/abcd/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.combatgam_*.rds
#
# Submit:
#   sbatch sbatch/run_abcd_combatgam_SA7_CV75_container.sbatch

#SBATCH -J abcd_combatgam_sa7_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH -o outputs/logs/combat_gam/abcd_combatgam_SA7_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

INPUT_RDS="${INPUT_RDS:-/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.merge.rds}"
OUTPUT_DIR="${OUTPUT_DIR:-outputs/results/combat_gam/abcd}"
TEST_N="${TEST_N:-0}"

cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  exit 1
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] OUTPUT_DIR=${OUTPUT_DIR}"
echo "[INFO] TEST_N=${TEST_N}"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
    "${sif}" "$@"
}

exec_sif Rscript --version
exec_sif Rscript --vanilla /work/combat_gam/scripts/run_abcd_nonlinear_combat_gam.R \
  "${INPUT_RDS}" \
  "/work/${OUTPUT_DIR}" \
  "${TEST_N}"

