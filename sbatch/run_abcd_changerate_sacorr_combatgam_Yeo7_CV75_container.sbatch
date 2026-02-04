#!/bin/bash
# ABCD | changerate S-A corr along age | Yeo7 | CV75 | container
#
# Note:
# - Default DS_RES=6 (edges=21). If your Yeo7 has 28 edges, set DS_RES=7.
#
# Submit:
#   sbatch sbatch/run_abcd_changerate_sacorr_combatgam_Yeo7_CV75_container.sbatch

#SBATCH -J abcd_sa_corr_age_yeo7_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/4th_changerate_SAcorr_abcd_Yeo7_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-6}"
OUT_TAG="${OUT_TAG:-Yeo7}"
N_EDGES="${N_EDGES:-}"
N_DRAWS="${N_DRAWS:-1000}"

# 0/1
FORCE="${FORCE:-0}"

INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"

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
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] OUT_TAG=${OUT_TAG}"
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

deriv_dir="${project_root}/outputs/results/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV${CVTHR}"
elementnum="$((DS_RES*(DS_RES+1)/2))"
deriv_df_expected="${deriv_dir}/derivative.df${elementnum}_CV${CVTHR}.rds"
deriv_post_expected="${deriv_dir}/derivative.posterior.df.SA${DS_RES}_CV${CVTHR}.rds"
if [ "${FORCE}" = "1" ] || [ ! -f "${deriv_df_expected}" ] || [ ! -f "${deriv_post_expected}" ]; then
  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    --input_rds="${INPUT_RDS_CONTAINER}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --force="${FORCE}"

  exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_ABCD_combatgam.R \
    --project_root=/work \
    --cvthr="${CVTHR}" \
    --ds_res="${DS_RES}" \
    ${N_EDGES:+--n_edges=${N_EDGES}} \
    --skip_posterior=0 \
    --force="${FORCE}"
else
  echo "[INFO] Found derivative outputs; skipping upstream derivative generation."
fi

exec_sif Rscript --vanilla /work/development_script/4th_changerate_SAcorr/S1st_SAcorr_alongAge_ABCD.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --out_tag="${OUT_TAG}" \
  --n_draws="${N_DRAWS}" \
  --force="${FORCE}"

