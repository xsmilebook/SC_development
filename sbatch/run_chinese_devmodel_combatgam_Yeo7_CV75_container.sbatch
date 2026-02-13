#!/bin/bash
# Chinese Cohort | Developmental model (S1-S4) | Yeo7 | CV75 | container
#
# Submit:
#   sbatch sbatch/run_chinese_devmodel_combatgam_Yeo7_CV75_container.sbatch

#SBATCH -J chinese_devmodel_yeo7_cg_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH -o outputs/logs/chinese_devmodel_combatgam_Yeo7_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
# Default to the "stable" SIF path; allow overriding for testing a newly built image.
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
YEO="${YEO:-7}"

INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/chinese/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"

N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"
MAKE_MATRIX_GRAPHS_YEO="${MAKE_MATRIX_GRAPHS_YEO:-1}"

mkdir -p "${project_root}/outputs/logs"
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

binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"

if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_CONTAINER="/work/${INPUT_RDS#${project_root}/}"
else
  INPUT_RDS_CONTAINER="${INPUT_RDS}"
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] YEO=${YEO}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS_YEO=${MAKE_MATRIX_GRAPHS_YEO}"
echo "[INFO] NOTE: S3/S4 always run with --force=1 to re-draw figures."

exec_sif() {
  singularity exec --cleanenv \
    --bind "${binds}" \
    --pwd /work \
    "${sif}" "$@"
}

exec_sif Rscript --version

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V1st_fitgammodels_Yeo_sumSCinvnode_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  --input_rds="${INPUT_RDS_CONTAINER}" \
  ${N_EDGES:+--n_edges=${N_EDGES}}

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V2nd_calculatederivative_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --skip_posterior="${SKIP_POSTERIOR}"

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V3rd_visualizationfitSCcurves_Yeo_sumSCinvnode_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  --force=1

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V4th_correlationTo_SArank_Yeo_sumSCinvnode_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  --make_matrix_graphs="${MAKE_MATRIX_GRAPHS_YEO}" \
  --force=1
