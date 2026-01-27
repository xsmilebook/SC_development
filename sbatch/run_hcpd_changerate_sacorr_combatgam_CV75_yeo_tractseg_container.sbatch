#!/bin/bash
#SBATCH -J hcpd_sa_corr_age_yeo_ts_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH -o outputs/logs/4th_changerate_SAcorr_hcpd_yeo_tractseg_CV75_%j.log

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

# Avoid OpenBLAS/OMP over-threading inside the container when using many cores.
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1
export SINGULARITYENV_VECLIB_MAXIMUM_THREADS=1
export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-72}"

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

run_yeo() {
  yeo="$1"
  echo "[INFO] ===== Yeo${yeo} ====="

  posterior_rds_host="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/yeo/Yeo${yeo}/combat_gam/CV${CVTHR}/derivative.posterior.df.Yeo${yeo}_CV${CVTHR}.rds"
  if [ "${FORCE}" = "1" ] || [ ! -f "${posterior_rds_host}" ]; then
    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V1st_fitgammodels_Yeo_sumSCinvnode_HCPD_combatgam.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      --yeo="${yeo}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --force="${FORCE}"

    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_Yeo_network/V2nd_calculatederivative_HCPD_combatgam.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      --yeo="${yeo}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --skip_posterior=0 \
      --force="${FORCE}"
  else
    echo "[INFO] Found Yeo${yeo} posterior derivatives; skipping upstream derivative generation."
  fi

  exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_Yeo_network/V1st_SAcorr_alongAge_Yeo_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    --n_draws="${N_DRAWS}" \
    --derivative_posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/yeo/Yeo${yeo}/combat_gam/CV${CVTHR}/derivative.posterior.df.Yeo${yeo}_CV${CVTHR}.rds" \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_Yeo_network/V2nd_fitgammodels_Yeo_ageseperate_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --yeo="${yeo}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"
}

run_tractseg() {
  echo "[INFO] ===== TractSeg ====="

  posterior_rds_host="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/tractseg/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}_TractSeg.rds"
  if [ "${FORCE}" = "1" ] || [ ! -f "${posterior_rds_host}" ]; then
    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD_TractSeg_combatgam.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --force="${FORCE}"

    exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_TractSeg/S2nd_calculatederivative_HCPD_TractSeg_combatgam.R \
      --project_root=/work \
      --cvthr="${CVTHR}" \
      ${N_EDGES:+--n_edges=${N_EDGES}} \
      --skip_posterior=0 \
      --force="${FORCE}"
  else
    echo "[INFO] Found TractSeg posterior derivatives; skipping upstream derivative generation."
  fi

  exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_TractSeg/S1st_SAcorr_alongAge_HCPD_TractSeg_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --n_draws="${N_DRAWS}" \
    --derivative_posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/tractseg/combat_gam/CV${CVTHR}/derivative.posterior.df.SA12_CV${CVTHR}_TractSeg.rds" \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/V_TractSeg/S2nd_fitgammodels_TractSeg_ageseperate_HCPD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"
}

run_yeo 7
run_yeo 17
run_tractseg

