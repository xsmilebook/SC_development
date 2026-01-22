#!/usr/bin/env bash
#SBATCH -J combat_chinese_gam
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --output=outputs/logs/combat_gam/combat_chinese_gam_%j.log

set -euo pipefail

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge_new.rds}
OUTPUT_RDS=${2:-outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}
TEST_N=${3:-0}

combat_gam/scripts/run_chinese_combat_gam.sh \
  "${OUTPUT_RDS}" \
  "${TEST_N}" \
  "${INPUT_RDS}"
