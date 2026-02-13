#!/bin/bash
# ABCD | Developmental model (S1-S4) | SA7 | CV75 | container
#
# Pipeline:
# - S1: fit GAMM models (raw + scaled)
# - S2: derivatives + posterior derivatives
# - S3: visualization + decile trajectories
# - S4: S-A rank correlation + scatter plots (optional distance control)
#
# Submit:
#   sbatch sbatch/run_abcd_devmodel_combatgam_SA7_CV75_container.sbatch

#SBATCH -J abcd_devmodel_sa7_cg_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/abcd_devmodel_combatgam_SA7_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-7}"
OUT_TAG="${OUT_TAG:-SA7}"

INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"
EUCLID_CSV="${EUCLID_CSV:-}"

N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"
MAKE_MATRIX_GRAPHS="${MAKE_MATRIX_GRAPHS:-0}"

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

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-50}"
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"

if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_CONTAINER="/work/${INPUT_RDS#${project_root}/}"
else
  INPUT_RDS_CONTAINER="${INPUT_RDS}"
fi

EUCLID_CONTAINER=""
if [ -n "${EUCLID_CSV:-}" ]; then
  if [[ "${EUCLID_CSV}" == "${project_root}/"* ]]; then
    EUCLID_CONTAINER="/work/${EUCLID_CSV#${project_root}/}"
  else
    EUCLID_CONTAINER="${EUCLID_CSV}"
  fi
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] OUT_TAG=${OUT_TAG}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] EUCLID_CSV=${EUCLID_CSV:-<unset>}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"
echo "[INFO] NOTE: S3/S4 always run with --force=1 to re-draw figures."

exec_sif() {
  singularity exec --cleanenv \
    --bind "${binds}" \
    --pwd /work \
    "${sif}" "$@"
}

exec_sif Rscript --version

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --input_rds="${INPUT_RDS_CONTAINER}" \
  ${N_EDGES:+--n_edges=${N_EDGES}}

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --skip_posterior="${SKIP_POSTERIOR}"

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --out_tag="${OUT_TAG}" \
  --force=1

cmd=(
  Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_ABCD_combatgam.R
  --project_root=/work
  --cvthr="${CVTHR}"
  --ds_res="${DS_RES}"
  --out_tag="${OUT_TAG}"
  --input_rds="${INPUT_RDS_CONTAINER}"
  --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}"
  --force=1
)
if [ -n "${EUCLID_CONTAINER:-}" ]; then
  cmd+=(--euclid_csv="${EUCLID_CONTAINER}")
fi
exec_sif "${cmd[@]}"

