#!/bin/bash
#SBATCH -J abcd_cog_pfac_nobeh
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=180G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/abcd_cognition_pfactor_nobeh_useexisting_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

# Existing ComBat outputs (no behavior covariates).
COG_INPUT_RDS="${COG_INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/baseline/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_neuroharmonize_baseline_age_sex_meanfd.rds}"
PFACTOR_INPUT_RDS="${PFACTOR_INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"

ABCD_EUCLID_CSV_HOST="${ABCD_EUCLID_CSV:-${project_root}/wd/interdataFolder_ABCD/average_EuclideanDistance_12.csv}"
ABCD_SA12_CSV_HOST="${ABCD_SA12_CSV:-${project_root}/wd/interdataFolder_ABCD/SA12_10.csv}"
ABCD_PLOTDATASUM_RDS_HOST="${ABCD_PLOTDATASUM_RDS:-${project_root}/outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds}"

COG_ASSOC_MODE="${COG_ASSOC_MODE:-original}"
COG_ASSOC_TAG="${COG_ASSOC_TAG:-nobehavior_covariates_existing}"
PFACTOR_TAG="${PFACTOR_TAG:-nobehavior_covariates_existing}"
FORCE="${FORCE:-0}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

if [ ! -f "${COG_INPUT_RDS}" ]; then
  echo "[ERROR] Missing COG_INPUT_RDS: ${COG_INPUT_RDS}"
  exit 1
fi
if [ ! -f "${PFACTOR_INPUT_RDS}" ]; then
  echo "[ERROR] Missing PFACTOR_INPUT_RDS: ${PFACTOR_INPUT_RDS}"
  exit 1
fi
if [ ! -f "${ABCD_EUCLID_CSV_HOST}" ]; then
  echo "[ERROR] Missing ABCD_EUCLID_CSV: ${ABCD_EUCLID_CSV_HOST}"
  exit 1
fi
if [ ! -f "${ABCD_SA12_CSV_HOST}" ]; then
  echo "[ERROR] Missing ABCD_SA12_CSV: ${ABCD_SA12_CSV_HOST}"
  exit 1
fi
if [ ! -f "${ABCD_PLOTDATASUM_RDS_HOST}" ]; then
  echo "[ERROR] Missing ABCD_PLOTDATASUM_RDS: ${ABCD_PLOTDATASUM_RDS_HOST}"
  exit 1
fi

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif"
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

binds="${project_root}:/work,/ibmgpfs:/ibmgpfs,/GPFS:/GPFS"

ABCD_EUCLID_CSV_CONTAINER="/work/${ABCD_EUCLID_CSV_HOST#${project_root}/}"
ABCD_SA12_CSV_CONTAINER="/work/${ABCD_SA12_CSV_HOST#${project_root}/}"
ABCD_PLOTDATASUM_RDS_CONTAINER="/work/${ABCD_PLOTDATASUM_RDS_HOST#${project_root}/}"
COG_INPUT_RDS_CONTAINER="/work/${COG_INPUT_RDS#${project_root}/}"
PFACTOR_INPUT_RDS_CONTAINER="/work/${PFACTOR_INPUT_RDS#${project_root}/}"

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-72}"
export SINGULARITYENV_ABCD_EUCLID_CSV="${ABCD_EUCLID_CSV_CONTAINER}"
export SINGULARITYENV_ABCD_SA12_CSV="${ABCD_SA12_CSV_CONTAINER}"
export SINGULARITYENV_ABCD_PLOTDATASUM_RDS="${ABCD_PLOTDATASUM_RDS_CONTAINER}"
export SINGULARITYENV_COG_INPUT_RDS="${COG_INPUT_RDS_CONTAINER}"
export SINGULARITYENV_PFACTOR_INPUT_RDS="${PFACTOR_INPUT_RDS_CONTAINER}"
export SINGULARITYENV_COG_ASSOC_MODE="${COG_ASSOC_MODE}"
export SINGULARITYENV_COG_ASSOC_TAG="${COG_ASSOC_TAG}"
export SINGULARITYENV_PFACTOR_TAG="${PFACTOR_TAG}"
export SINGULARITYENV_FORCE="${FORCE}"
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

singularity exec --cleanenv \
  --bind "${binds}" \
  --pwd /work \
  "${sif}" \
  bash -lc 'Rscript --vanilla /work/development_script/5th_cognition/run_abcd_cognition_comp_agecorrected_S1.R && Rscript --vanilla /work/development_script/6th_pfactor/run_abcd_pfactor_effect_continuous_S1.R'
