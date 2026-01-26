#!/bin/bash
#SBATCH -J abcd_cog_compac_sexmfd
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=180G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/abcd_cognition_comp_agecorrected_sex_meanfd_%j.log

set -euo pipefail

# This is a thin wrapper of `sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch`
# to run the `sex_meanfd` association variant (no `s(age)`, includes `sex` + `mean_fd`).

export COG_ASSOC_MODE="${COG_ASSOC_MODE:-sex_meanfd}"
export COG_ASSOC_TAG="${COG_ASSOC_TAG:-sex_meanfd}"

exec bash sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch
