#!/bin/bash
#SBATCH -J abcd_age_lmm_pf
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 40
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH -o outputs/logs/abcd_age_lmm_random_slope_pfactor_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

N_EDGES="${N_EDGES:-78}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_abcd_age_lmm_random_slope_pfactor_groups_container.sbatch"
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-40}"
export SINGULARITYENV_N_EDGES="${N_EDGES}"
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

singularity exec --cleanenv \
  --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
  --pwd /work \
  "${sif}" \
  bash -lc 'Rscript --vanilla /work/development_script/6th_pfactor/run_abcd_age_lmm_random_slope_pfactor_groups.R'
