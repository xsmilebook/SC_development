#!/usr/bin/env bash
#SBATCH -J combat_hcpd_gam
#SBATCH -p q_fat_c
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=outputs/logs/combat_gam/combat_hcpd_gam_%j.log

set -euo pipefail

combat_gam/scripts/run_hcpd_combat_gam.sh \
  outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds
