#!/bin/bash
#SBATCH -J abcd_pfactor
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=180G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/abcd_pfactor_effect_continuous_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

ABCD_EUCLID_CSV_HOST="${ABCD_EUCLID_CSV:-${project_root}/wd/interdataFolder_ABCD/average_EuclideanDistance_12.csv}"
ABCD_PLOTDATASUM_RDS_HOST="${ABCD_PLOTDATASUM_RDS:-${project_root}/outputs/intermediate/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV75/plotdatasum.df_SA12_sumSCinvnode_siteall_CV75.rds}"
FORCE="${FORCE:-0}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_abcd_pfactor_effect_continuous_container.sbatch"
  exit 1
fi

if [ ! -f "${ABCD_EUCLID_CSV_HOST}" ]; then
  echo "[ERROR] Missing ABCD_EUCLID_CSV: ${ABCD_EUCLID_CSV_HOST}"
  exit 1
fi

if [ ! -f "${ABCD_PLOTDATASUM_RDS_HOST}" ]; then
  echo "[ERROR] Missing ABCD_PLOTDATASUM_RDS: ${ABCD_PLOTDATASUM_RDS_HOST}"
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

binds="${project_root}:/work"
if [[ "${ABCD_EUCLID_CSV_HOST}" == "${project_root}/"* ]]; then
  ABCD_EUCLID_CSV_CONTAINER="/work/${ABCD_EUCLID_CSV_HOST#${project_root}/}"
else
  extra_dir="$(dirname "${ABCD_EUCLID_CSV_HOST}")"
  binds="${binds},${extra_dir}:${extra_dir}"
  ABCD_EUCLID_CSV_CONTAINER="${ABCD_EUCLID_CSV_HOST}"
fi

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-72}"
export SINGULARITYENV_ABCD_EUCLID_CSV="${ABCD_EUCLID_CSV_CONTAINER}"
export SINGULARITYENV_ABCD_PLOTDATASUM_RDS="${ABCD_PLOTDATASUM_RDS_HOST}"
export SINGULARITYENV_FORCE="${FORCE}"
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

singularity exec --cleanenv \
  --bind "${binds}" \
  --bind /ibmgpfs:/ibmgpfs \
  --bind /GPFS:/GPFS \
  --pwd /work \
  "${sif}" \
  bash -lc 'Rscript --vanilla /work/development_script/6th_pfactor/run_abcd_pfactor_effect_continuous_S1.R'
