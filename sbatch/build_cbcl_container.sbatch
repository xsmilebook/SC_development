#!/bin/bash
#SBATCH -J build_cbcl_container
#SBATCH -p q_fat_c
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o outputs/logs/build_cbcl_container_%j.log

set -euo pipefail

module load singularity/3.7.0

singularity build --remote containers/scdevelopment_cbcl.sif containers/scdevelopment_cbcl.def
