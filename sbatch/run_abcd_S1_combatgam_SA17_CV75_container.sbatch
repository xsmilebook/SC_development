#!/bin/bash
# ABCD | Developmental model | S1 (fit GAMM models) | SA17 | CV75
#
# Prerequisite:
# - ComBat-GAM output exists:
#   outputs/results/combat_gam/abcd/SCdata_SA17_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds
#
# Submit:
#   sbatch sbatch/run_abcd_S1_combatgam_SA17_CV75_container.sbatch

#SBATCH -J abcd_S1_sa17_cg_cv75
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH --mem=160G
#SBATCH -t 36:00:00
#SBATCH -o outputs/logs/abcd_S1_combatgam_SA17_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
DS_RES="${DS_RES:-17}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_SA17_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"
N_EDGES="${N_EDGES:-}"
FORCE="${FORCE:-0}"

cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  exit 1
fi

if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_ARG="/work/${INPUT_RDS#${project_root}/}"
else
  INPUT_RDS_ARG="${INPUT_RDS}"
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] DS_RES=${DS_RES}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] FORCE=${FORCE}"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work,/ibmgpfs:/ibmgpfs" \
    "${sif}" "$@"
}

exec_sif Rscript --version
exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --ds_res="${DS_RES}" \
  --input_rds="${INPUT_RDS_ARG}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --force="${FORCE}"

