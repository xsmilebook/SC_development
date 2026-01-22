#!/bin/bash
#SBATCH -J combat_chinese_gam_native
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --output=outputs/logs/combat_gam/combat_chinese_gam_native_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/congjing/double_check_scdevelopment/NC/interdataFolder_ChineseCohort/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge_new.rds}
OUTPUT_RDS=${2:-${project_root}/outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_native.rds}
TEST_N=${3:-0}

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch combat_gam/sbatch/chinese_combat_gam_native.sbatch"
  exit 1
fi

export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

map_path() {
  local p="$1"
  if [[ "${p}" == "${project_root}"/* ]]; then
    printf "/work/%s" "${p#${project_root}/}"
  else
    printf "%s" "${p}"
  fi
}

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work" \
    --bind /ibmgpfs:/ibmgpfs \
    --bind /GPFS:/GPFS \
    "${sif}" "$@"
}

exec_sif python3 --version

exec_sif python3 /work/combat_gam/scripts/run_combat_gam_neuroharmonize_native.py \
  --input-rds "${INPUT_RDS}" \
  --output-rds "$(map_path "${OUTPUT_RDS}")" \
  --batch-col study \
  --id-col subID \
  --age-col Age \
  --sex-col Sex \
  --mean-fd-col mean_fd \
  --smooth-df 4 \
  --test-n "${TEST_N}"
