#!/bin/bash
#SBATCH -J cbcl_assoc_full
#SBATCH -p q_fat_c
#SBATCH -c 50
#SBATCH --mem=120G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/cbcl_assoc_full_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment

Rscript development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD.R
