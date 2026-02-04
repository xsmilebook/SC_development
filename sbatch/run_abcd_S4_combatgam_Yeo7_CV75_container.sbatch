#!/bin/bash
# ABCD | Developmental model | S4 (SCrank correlation + plots) | Yeo7 | CV75
#
# Note:
# - DS_RES is the number of nodes used to build the lower-triangle edge list.
#   For Yeo7 variants, confirm SC edge count matches DS_RES*(DS_RES+1)/2.
#
# Prerequisite:
# - S1 outputs under: outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/
#
# Submit:
#   sbatch sbatch/run_abcd_S4_combatgam_Yeo7_CV75_container.sbatch

#SBATCH -J abcd_S4_yeo7_cg_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 20
#SBATCH --mem=60G
#SBATCH -t 12:00:00
#SBATCH -o outputs/logs/abcd_S4_combatgam_Yeo7_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-7}"
OUT_TAG="${OUT_TAG:-Yeo7}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"
MAKE_MATRIX_GRAPHS="${MAKE_MATRIX_GRAPHS:-0}"
FORCE="${FORCE:-0}"

# Optional (default: unset -> skip control-distance residualization)
EUCLID_CSV="${EUCLID_CSV:-}"

cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  exit 1
fi

if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_ARG="/work/${INPUT_RDS#${project_root}/}"
else
  INPUT_RDS_ARG="${INPUT_RDS}"
fi

EUCLID_ARG=""
if [ -n "${EUCLID_CSV:-}" ]; then
  if [[ "${EUCLID_CSV}" == "${project_root}/"* ]]; then
    EUCLID_ARG="/work/${EUCLID_CSV#${project_root}/}"
  else
    EUCLID_ARG="${EUCLID_CSV}"
  fi
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] OUT_TAG=${OUT_TAG}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] EUCLID_CSV=${EUCLID_CSV:-<unset>}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"
echo "[INFO] FORCE=${FORCE}"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
    "${sif}" "$@"
}

cmd=(
  Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_ABCD_combatgam.R
  --project_root=/work
  --cvthr="${CVTHR}"
  --ds_res="${DS_RES}"
  --out_tag="${OUT_TAG}"
  --input_rds="${INPUT_RDS_ARG}"
  --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}"
  --force="${FORCE}"
)
if [ -n "${EUCLID_ARG}" ]; then
  cmd+=(--euclid_csv="${EUCLID_ARG}")
fi

exec_sif "${cmd[@]}"

