#!/bin/bash
#SBATCH -J hcpd_devmodel_yeo_trk_cg_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH -o outputs/logs/hcpd_devmodel_yeo_tractseg_combatgam_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
# Default to the "stable" SIF path; allow overriding for testing a newly built image.
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"

YEO7_INPUT_RDS="${YEO7_INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"
YEO17_INPUT_RDS="${YEO17_INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"
TRACTSEG_INPUT_RDS="${TRACTSEG_INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.TractSeg.combatgam.rds}"

N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"

# Matrix graphs are part of the original Yeo/TractSeg scripts; enable by default.
MAKE_MATRIX_GRAPHS_YEO="${MAKE_MATRIX_GRAPHS_YEO:-1}"
MAKE_MATRIX_GRAPHS_TRACTSEG="${MAKE_MATRIX_GRAPHS_TRACTSEG:-1}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_hcpd_devmodel_combatgam_CV75_yeo_tractseg_container.sbatch"
  exit 1
fi

# Prevent thread oversubscription / pthread_create failures (OpenBLAS/MKL/OMP).
# NOTE: exec_sif uses `singularity exec --cleanenv`, so pass env vars via
# SINGULARITYENV_* for them to be visible inside the container.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export SINGULARITYENV_OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}"
export SINGULARITYENV_OMP_NUM_THREADS="${OMP_NUM_THREADS}"
export SINGULARITYENV_MKL_NUM_THREADS="${MKL_NUM_THREADS}"
export SINGULARITYENV_NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS}"
export SINGULARITYENV_VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS}"
export SINGULARITYENV_BLIS_NUM_THREADS="${BLIS_NUM_THREADS}"

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] YEO7_INPUT_RDS=${YEO7_INPUT_RDS}"
echo "[INFO] YEO17_INPUT_RDS=${YEO17_INPUT_RDS}"
echo "[INFO] TRACTSEG_INPUT_RDS=${TRACTSEG_INPUT_RDS}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS_YEO=${MAKE_MATRIX_GRAPHS_YEO}"
echo "[INFO] MAKE_MATRIX_GRAPHS_TRACTSEG=${MAKE_MATRIX_GRAPHS_TRACTSEG}"
echo "[INFO] OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} OMP_NUM_THREADS=${OMP_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS}"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work" \
    "${sif}" "$@"
}

exec_sif Rscript --version

run_yeo() {
  local yeo="$1"
  local input_rds="$2"

  local force_s3="${3:-1}"
  local force_s4="${4:-1}"

  echo "[INFO] === HCPD Yeo${yeo} (ComBat-GAM) ==="
  echo "[INFO] input_rds=${input_rds}"
  echo "[INFO] FORCE_S3=${force_s3} FORCE_S4=${force_s4}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V1st_fitgammodels_Yeo_sumSCinvnode_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    --input_rds="/work/${input_rds#${project_root}/}" \
    ${N_EDGES:+--n_edges=${N_EDGES}}

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V2nd_calculatederivative_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior="${SKIP_POSTERIOR}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V3rd_visualizationfitSCcurves_Yeo_sumSCinvnode_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    --force="${force_s3}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V4th_correlationTo_SArank_Yeo_sumSCinvnode_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    --make_matrix_graphs="${MAKE_MATRIX_GRAPHS_YEO}" \
    --force="${force_s4}"
}

run_tractseg() {
  local input_rds="$1"
  local force_s3="${2:-1}"
  local force_s4="${3:-1}"

  echo "[INFO] === HCPD TractSeg (major-bundle, ComBat-GAM) ==="
  echo "[INFO] input_rds=${input_rds}"
  echo "[INFO] FORCE_S3=${force_s3} FORCE_S4=${force_s4}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_TractSeg_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --input_rds="/work/${input_rds#${project_root}/}" \
    ${N_EDGES:+--n_edges=${N_EDGES}}

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S2nd_calculatederivative_HCPD_TractSeg_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior="${SKIP_POSTERIOR}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_HCPD_TractSeg_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --force="${force_s3}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S4th_correlationTo_SArank_SA12sumSCinvnode_HCPD_TractSeg_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --make_matrix_graphs="${MAKE_MATRIX_GRAPHS_TRACTSEG}" \
    --force="${force_s4}"
}

echo "[INFO] NOTE: S3/S4 always run with --force=1 to re-draw figures."

run_yeo 7 "${YEO7_INPUT_RDS}" 1 1
run_yeo 17 "${YEO17_INPUT_RDS}" 1 1
run_tractseg "${TRACTSEG_INPUT_RDS}" 1 1
