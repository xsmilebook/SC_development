#!/bin/bash
#SBATCH -J hcpd_sa_corr_age_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_hcpd_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1
FORCE="${FORCE:-0}"

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

# Ensure SLURM core count is visible inside the container under --cleanenv.
export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-50}"
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] N_DRAWS=${N_DRAWS}"
echo "[INFO] FORCE=${FORCE}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

exec_sif() {
  binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"
  singularity exec --cleanenv \
    --bind "${binds}" \
    --pwd /work \
    "${sif}" "$@"
}

exec_sif Rscript --version

# Ensure derivatives exist (skip if already present unless FORCE=1)
posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}.rds"
if [ "${FORCE}" = "1" ] || [ ! -f "${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV${CVTHR}/derivative.df78_CV${CVTHR}.rds" ] || [ ! -f "${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}.rds" ]; then
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior=0 \
    --force="${FORCE}"
else
  echo "[INFO] Found derivative outputs; skipping upstream derivative generation."
fi

exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --n_draws="${N_DRAWS}" \
  --derivative_posterior_rds="${posterior_rds}" \
  --force="${FORCE}"
