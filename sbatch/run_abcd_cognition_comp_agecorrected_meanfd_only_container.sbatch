#!/bin/bash
#SBATCH -J abcd_cog_compac_mfd
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=180G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/abcd_cognition_comp_agecorrected_meanfd_only_%j.log

set -euo pipefail

# This is a thin wrapper of `sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch`
# to run the `meanfd_only` association variant (no `s(age)`, no `sex`, only `mean_fd`).

export COG_ASSOC_MODE="${COG_ASSOC_MODE:-meanfd_only}"
export COG_ASSOC_TAG="${COG_ASSOC_TAG:-meanfd_only}"

exec bash sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch

