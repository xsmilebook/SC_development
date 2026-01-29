#!/bin/bash
#SBATCH -J hcpd_devmodel_sa7_cg_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/hcpd_devmodel_combatgam_SA7_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-7}"
OUT_TAG="${OUT_TAG:-SA7}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"
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

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] OUT_TAG=${OUT_TAG}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"
echo "[INFO] NOTE: S3/S4 always run with --force=1 to re-draw figures."

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
    "${sif}" "$@"
}

exec_sif Rscript --version

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --input_rds="/work/${INPUT_RDS#${project_root}/}" \
  ${N_EDGES:+--n_edges=${N_EDGES}}

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --skip_posterior="${SKIP_POSTERIOR}"

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --force=1

# No Euclidean-distance control for SA7 by default (pass EUCLID_CSV to enable).
# NOTE: avoid expanding an empty bash array under `set -u` (some environments treat it as unbound).
EUCLID_CSV="${EUCLID_CSV:-}"
if [ -n "${EUCLID_CSV}" ]; then
  if [[ "${EUCLID_CSV}" == "${project_root}/"* ]]; then
    EUCLID_CSV_ARG="/work/${EUCLID_CSV#${project_root}/}"
  else
    EUCLID_CSV_ARG="${EUCLID_CSV}"
  fi
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    --out_tag="${OUT_TAG}" \
    --sa_axis_mode=addsquare \
    --input_rds="/work/${INPUT_RDS#${project_root}/}" \
    --euclid_csv="${EUCLID_CSV_ARG}" \
    --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}" \
    --force=1
else
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    --out_tag="${OUT_TAG}" \
    --sa_axis_mode=addsquare \
    --input_rds="/work/${INPUT_RDS#${project_root}/}" \
    --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}" \
    --force=1
fi
