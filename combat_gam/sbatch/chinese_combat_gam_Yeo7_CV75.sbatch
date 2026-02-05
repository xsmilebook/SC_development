#!/usr/bin/env bash
#SBATCH -J combat_chinese_yeo7_gam
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --output=outputs/logs/combat_gam/combat_chinese_yeo7_gam_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"

INPUT_RDS=${1:-${project_root}/outputs/intermediate/1st_dataclean/chinese_cohort/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.merge.rds}
OUTPUT_RDS=${2:-${project_root}/outputs/results/combat_gam/chinese/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}
TEST_N=${3:-0}

mkdir -p "${project_root}/outputs/logs/combat_gam" "${project_root}/outputs/results/combat_gam/chinese"

combat_gam/scripts/run_chinese_combat_gam.sh \
  "${OUTPUT_RDS}" \
  "${TEST_N}" \
  "${INPUT_RDS}"
