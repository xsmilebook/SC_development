#!/bin/bash
#SBATCH -J hcpd_sa_corr_age_sa7_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_hcpd_SA7_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-7}"
OUT_TAG="${OUT_TAG:-SA7}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1
FORCE="${FORCE:-0}"

INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds}"

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
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] OUT_TAG=${OUT_TAG}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
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
deriv_dir="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV${CVTHR}"
if [ "${FORCE}" = "1" ] || [ ! -f "${deriv_dir}/derivative.df28_CV${CVTHR}.rds" ] || [ ! -f "${deriv_dir}/derivative.posterior.df.SA${DS_RES}_CV${CVTHR}.rds" ]; then
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    --input_rds="/work/${INPUT_RDS#${project_root}/}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_HCPD.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior=0 \
    --force="${FORCE}"
else
  echo "[INFO] Found derivative outputs; skipping upstream derivative generation."
fi

posterior_rds="/work/outputs/results/2nd_fitdevelopmentalmodel/hcpd/combat_gam/CV${CVTHR}/derivative.posterior.df.SA${DS_RES}_CV${CVTHR}.rds"
exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --out_tag="${OUT_TAG}" \
  --sa_axis_mode=addsquare \
  --do_control_distance=0 \
  --save_svg=0 \
  --n_draws="${N_DRAWS}" \
  --derivative_posterior_rds="${posterior_rds}" \
  --force="${FORCE}"

flip_csv="/work/outputs/results/4th_changerate_SAcorr/hcpd/combat_gam/CV${CVTHR}_${OUT_TAG}/alignment_summary_flip_age.csv"
exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/S2nd_fitgammodels_SA12sumSCinvnode_ageseperate_HCPD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --out_tag="${OUT_TAG}" \
  --sa_axis_mode=addsquare \
  --save_svg=0 \
  --flip_age_csv="${flip_csv}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --force="${FORCE}"

