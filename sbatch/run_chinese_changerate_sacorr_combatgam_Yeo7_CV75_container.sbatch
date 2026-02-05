#!/bin/bash
# Chinese Cohort | changerate S-A corr along age | Yeo7 | CV75 | container
#
# Submit:
#   sbatch sbatch/run_chinese_changerate_sacorr_combatgam_Yeo7_CV75_container.sbatch

#SBATCH -J chinese_sa_corr_age_yeo7_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_chinese_Yeo7_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
YEO="${YEO:-7}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1
FORCE="${FORCE:-0}"

INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/chinese/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"

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

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] YEO=${YEO}"
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

posterior_rds_host="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/chinese/yeo/Yeo${YEO}/combat_gam/CV${CVTHR}/derivative.posterior.df.Yeo${YEO}_CV${CVTHR}.rds"
if [ "${FORCE}" = "1" ] || [ ! -f "${posterior_rds_host}" ]; then
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V1st_fitgammodels_Yeo_sumSCinvnode_ChineseCohort_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${YEO}" \
    --input_rds="${INPUT_RDS_CONTAINER}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V2nd_calculatederivative_ChineseCohort_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${YEO}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior=0 \
    --force="${FORCE}"
else
  echo "[INFO] Found Yeo${YEO} posterior derivatives; skipping upstream derivative generation."
fi

exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_Yeo_network/V1st_SAcorr_alongAge_Yeo_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  --n_draws="${N_DRAWS}" \
  --derivative_posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/chinese/yeo/Yeo${YEO}/combat_gam/CV${CVTHR}/derivative.posterior.df.Yeo${YEO}_CV${CVTHR}.rds" \
  --force="${FORCE}"

exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_Yeo_network/V2nd_fitgammodels_Yeo_ageseperate_ChineseCohort_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --yeo="${YEO}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --force="${FORCE}"
