#!/usr/bin/env bash
#SBATCH -J combat_abcd_gam
#SBATCH -p q_fat_c
#SBATCH --cpus-per-task=72
#SBATCH --output=outputs/logs/combat_gam/combat_abcd_gam_%j.log

set -euo pipefail

combat_gam/scripts/run_abcd_nonlinear_combat_gam.sh \
  /ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds \
  outputs/results/combat_gam/abcd
