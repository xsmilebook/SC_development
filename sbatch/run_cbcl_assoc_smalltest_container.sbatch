#!/bin/bash
#SBATCH -J cbcl_assoc_test
#SBATCH -p q_fat_c
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o outputs/logs/cbcl_assoc_test_%j.log

set -euo pipefail

module load singularity/3.7.0

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
cd "${project_root}"

singularity exec --bind /ibmgpfs:/ibmgpfs containers/scdevelopment_cbcl.sif \
  Rscript development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD_smalltest.R
