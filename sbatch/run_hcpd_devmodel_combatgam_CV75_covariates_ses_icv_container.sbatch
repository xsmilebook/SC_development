#!/bin/bash
#SBATCH -J hcpd_devmodel_cov_ses_icv_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/hcpd_devmodel_covariates_ses_icv_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"
DEMO_CSV="${DEMO_CSV:-${project_root}/demopath/HCPD_demo_behav.csv}"

N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"

# Matrix graphs in S4 (tiff only); enable by default.
MAKE_MATRIX_GRAPHS="${MAKE_MATRIX_GRAPHS:-1}"

# Optional recompute for heavy steps:
FORCE_S1="${FORCE_S1:-0}"
FORCE_S2="${FORCE_S2:-0}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_hcpd_devmodel_combatgam_CV75_covariates_ses_icv_container.sbatch"
  exit 1
fi

# Prevent thread oversubscription / pthread_create failures (OpenBLAS/MKL/OMP).
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
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] DEMO_CSV=${DEMO_CSV}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"
echo "[INFO] FORCE_S1=${FORCE_S1} FORCE_S2=${FORCE_S2}"
echo "[INFO] NOTE: S3/S4 always run with --force=1 to re-draw figures."

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work" \
    "${sif}" "$@"
}

exec_sif Rscript --version

run_cov() {
  local cov_tag="$1"
  local add_covariate="$2"

  echo "[INFO] === HCPD SA12 devmodel (ComBat-GAM) | cov_tag=${cov_tag} add_covariate=${add_covariate} ==="

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_covariates.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --input_rds="/work/${INPUT_RDS#${project_root}/}" \
    --demo_csv="/work/${DEMO_CSV#${project_root}/}" \
    --add_covariate="${add_covariate}" \
    --cov_tag="${cov_tag}" \
    --force="${FORCE_S1}" \
    ${N_EDGES:+--n_edges=${N_EDGES}}

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V2nd_calculatederivative_HCPD_covariates.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --cov_tag="${cov_tag}" \
    --force="${FORCE_S2}" \
    --skip_posterior="${SKIP_POSTERIOR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}}

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V3rd_visualizationfitSCcurves_totalstrength_ICV_sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --cov_tag="${cov_tag}" \
    --force=1

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V4th_correlationTo_SArank_control_covariates_sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --cov_tag="${cov_tag}" \
    --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}" \
    --force=1
}

# Two variants:
# - SES: add_covariate=income.adj
# - ICV: add_covariate=ICV
run_cov "SES" "income.adj"
run_cov "ICV" "ICV"

