#!/bin/bash
# Chinese Cohort | Merge SC + phenotype | SA7 + SA17 + Yeo7 | container
#
# This job generates the pre-ComBat merged RDS files (CV75 and CV25; invnode and raw sumSC)
# under the current SCDevelopment project.
#
# Outputs (default):
# - outputs/intermediate/1st_dataclean/chinese_cohort/SCdata_SA7_*.merge.rds
# - outputs/intermediate/1st_dataclean/chinese_cohort/SCdata_SA17_*.merge.rds
# - outputs/intermediate/1st_dataclean/chinese_cohort/SCdata_Yeo7_*.merge.rds
#
# Submit:
#   sbatch sbatch/run_chinese_merge_SA7_SA17_Yeo7_container.sbatch

#SBATCH -J chinese_merge_sa7_sa17_yeo7
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH -t 24:00:00
#SBATCH -o outputs/logs/chinese_merge_SA7_SA17_Yeo7_%j.log

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

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-20}"
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

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] OUT_DIR=${OUT_DIR}"
echo "[INFO] FIG_DIR=${FIG_DIR}"
echo "[INFO] MAKE_PLOTS=${MAKE_PLOTS}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

exec_sif Rscript --version

exec_sif Rscript --vanilla /work/development_script/1st_dataclean/merge_Chinese_cohort/S2nd_mergedata_SA_ds_sumSC_ChineseCohort.R \
  --project_root=/work \
  --ds_res=7 \
  --output_dir="/work/${OUT_DIR#${project_root}/}" \
  --figure_dir="/work/${FIG_DIR#${project_root}/}" \
  --make_plots="${MAKE_PLOTS}"

exec_sif Rscript --vanilla /work/development_script/1st_dataclean/merge_Chinese_cohort/S2nd_mergedata_SA_ds_sumSC_ChineseCohort.R \
  --project_root=/work \
  --ds_res=17 \
  --output_dir="/work/${OUT_DIR#${project_root}/}" \
  --figure_dir="/work/${FIG_DIR#${project_root}/}" \
  --make_plots="${MAKE_PLOTS}"

exec_sif Rscript --vanilla /work/development_script/1st_dataclean/merge_Chinese_cohort/V2nd_mergedata_Yeo_sumSC_ChineseCohort.R \
  --project_root=/work \
  --yeo_res=7 \
  --output_dir="/work/${OUT_DIR#${project_root}/}" \
  --figure_dir="/work/${FIG_DIR#${project_root}/}" \
  --make_plots="${MAKE_PLOTS}"
