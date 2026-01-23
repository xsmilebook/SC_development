#!/bin/bash
#SBATCH -J abcd_sa_corr_age_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_abcd_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1
FORCE="${FORCE:-0}"

# Default to the longitudinal ComBat-GAM dataset with age+sex+meanFD protected.
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"

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

# Ensure INPUT_RDS is visible inside the container. If it is under project_root,
# access it via /work/...; otherwise add an extra bind for its directory.
binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"
if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_CONTAINER="/work/${INPUT_RDS#${project_root}/}"
else
  extra_dir="$(dirname "${INPUT_RDS}")"
  binds="${binds},${extra_dir}:${extra_dir}"
  INPUT_RDS_CONTAINER="${INPUT_RDS}"
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] N_DRAWS=${N_DRAWS}"
echo "[INFO] FORCE=${FORCE}"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${binds}" \
    --pwd /work \
    "${sif}" "$@"
}

exec_sif Rscript --version

# Ensure derivatives exist (skip if already present unless FORCE=1)
deriv_dir="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV${CVTHR}"
if [ "${FORCE}" = "1" ] || [ ! -f "${deriv_dir}/derivative.df78_CV${CVTHR}.rds" ] || [ ! -f "${deriv_dir}/derivative.posterior.df.SA12_CV${CVTHR}.rds" ]; then
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --input_rds="${INPUT_RDS_CONTAINER}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_ABCD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior=0 \
    --force="${FORCE}"
else
  echo "[INFO] Found derivative outputs; skipping upstream derivative generation."
fi

exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_ABCD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --n_draws="${N_DRAWS}" \
  --force="${FORCE}"
