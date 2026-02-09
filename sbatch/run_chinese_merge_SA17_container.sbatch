#!/bin/bash
# Chinese Cohort | Merge SC + phenotype | SA17 | container
#
# Submit:
#   sbatch sbatch/run_chinese_merge_SA17_container.sbatch

#SBATCH -J chinese_merge_sa17
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=80G
#SBATCH -t 24:00:00
#SBATCH -o outputs/logs/chinese_merge_SA17_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

OUT_DIR="${OUT_DIR:-${project_root}/outputs/intermediate/1st_dataclean/chinese_cohort}"
FIG_DIR="${FIG_DIR:-${project_root}/outputs/figures/1st_dataclean/chinese_cohort}"
MAKE_PLOTS="${MAKE_PLOTS:-0}"

mkdir -p "${project_root}/outputs/logs"
mkdir -p "${OUT_DIR}"
mkdir -p "${FIG_DIR}"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  exit 1
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-72}"
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
    --pwd /work \
    "${sif}" "$@"
}

if [[ "${OUT_DIR}" != "${project_root}/"* ]]; then
  echo "[ERROR] OUT_DIR must be under project_root for container binding: OUT_DIR=${OUT_DIR}"
  exit 1
fi
if [[ "${FIG_DIR}" != "${project_root}/"* ]]; then
  echo "[ERROR] FIG_DIR must be under project_root for container binding: FIG_DIR=${FIG_DIR}"
  exit 1
fi

exec_sif Rscript --vanilla /work/development_script/1st_dataclean/merge_Chinese_cohort/S2nd_mergedata_SA_ds_sumSC_ChineseCohort.R \
  --project_root=/work \
  --ds_res=17 \
  --output_dir="/work/${OUT_DIR#${project_root}/}" \
  --figure_dir="/work/${FIG_DIR#${project_root}/}" \
  --make_plots="${MAKE_PLOTS}"
