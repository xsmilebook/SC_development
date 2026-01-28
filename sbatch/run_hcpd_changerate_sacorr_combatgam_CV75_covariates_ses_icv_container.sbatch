#!/bin/bash
#SBATCH -J hcpd_sacorr_cov_ses_icv_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_hcpd_covariates_SES_ICV_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1: force recompute of 4th outputs (and optionally upstream derivatives if FORCE_UPSTREAM=1)
FORCE="${FORCE:-0}"
FORCE_UPSTREAM="${FORCE_UPSTREAM:-0}"

# Heavy step control
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"

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

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] N_DRAWS=${N_DRAWS}"
echo "[INFO] FORCE=${FORCE} FORCE_UPSTREAM=${FORCE_UPSTREAM}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

exec_sif() {
  binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"
  singularity exec --cleanenv \
    --bind "${binds}" \
    --pwd /work \
    "${sif}" "$@"
}

exec_sif Rscript --version

ensure_cov_derivatives() {
  local cov_tag="$1"
  local add_covariate="$2"

  local deriv_df="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/covariates/${cov_tag}/combat_gam/CV${CVTHR}/derivative.df78_CV${CVTHR}.rds"
  local deriv_post="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/covariates/${cov_tag}/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}.rds"

  if [ "${FORCE_UPSTREAM}" = "1" ] || [ ! -f "${deriv_df}" ] || [ ! -f "${deriv_post}" ]; then
    echo "[INFO] Generating covariate derivatives for cov_tag=${cov_tag} (add_covariate=${add_covariate})"

    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_covariates.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      --cov_tag="${cov_tag}" \
      --add_covariate="${add_covariate}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --force="${FORCE_UPSTREAM}"

    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Covariates/V2nd_calculatederivative_HCPD_covariates.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      --cov_tag="${cov_tag}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --skip_posterior="${SKIP_POSTERIOR}" \
      --force="${FORCE_UPSTREAM}"
  else
    echo "[INFO] Found derivative outputs for cov_tag=${cov_tag}; skipping upstream generation."
  fi
}

run_cov_alignment() {
  local cov_tag="$1"

  local posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/covariates/${cov_tag}/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}.rds"
  local derivative_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/covariates/${cov_tag}/combat_gam/CV${CVTHR}/derivative.df78_CV${CVTHR}.rds"

  echo "[INFO] === 4th_changerate_SAcorr | HCP-D cov_tag=${cov_tag} ==="
  exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_Covariates/V1st_SAcorr_alongAge_HCPD_add_covariates.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --cov_tag="${cov_tag}" \
    --n_draws="${N_DRAWS}" \
    --derivative_posterior_rds="${posterior_rds}" \
    --derivative_rds="${derivative_rds}" \
    --force="${FORCE}"
}

# SES: income.adj
ensure_cov_derivatives "SES" "income.adj"
run_cov_alignment "SES"

# ICV: ICV
ensure_cov_derivatives "ICV" "ICV"
run_cov_alignment "ICV"

