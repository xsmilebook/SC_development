#!/usr/bin/env bash
#SBATCH -J combat_abcd_gam_nh_baseline
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --output=outputs/logs/combat_gam/abcd_combat_gam_neuroharmonize_baseline_%j.log

set -euo pipefail

combat_gam/scripts/run_abcd_combat_gam_neuroharmonize_baseline.sh \
  /ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds \
  outputs/results/combat_gam/abcd/baseline
