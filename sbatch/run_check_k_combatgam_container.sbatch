#!/bin/bash
# Check GAM k values using ComBat-GAM data (container)
#
# Submit:
#   sbatch sbatch/run_check_k_combatgam_container.sbatch

#SBATCH -J check_k_combatgam
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 32
#SBATCH -o outputs/logs/check_k_combatgam_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

DATASET="${DATASET:-chinese}"
K_VALUES="${K_VALUES:-3,4,5,6}"
N_EDGES="${N_EDGES:-}"
INPUT_RDS="${INPUT_RDS:-}"

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

export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-32}"
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

if [ -z "${INPUT_RDS}" ]; then
  case "${DATASET}" in
    hcpd)
      INPUT_RDS="${project_root}/outputs/results/combat_gam/hcpd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
      ;;
    abcd)
      INPUT_RDS="${project_root}/outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds"
      ;;
    chinese)
      INPUT_RDS="${project_root}/outputs/results/combat_gam/chinese/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
      ;;
    *)
      echo "[ERROR] Unsupported DATASET: ${DATASET} (expected hcpd|abcd|chinese)"
      exit 1
      ;;
  esac
fi

if [[ "${INPUT_RDS}" == "${project_root}/"* ]]; then
  INPUT_RDS_CONTAINER="/work/${INPUT_RDS#${project_root}/}"
else
  INPUT_RDS_CONTAINER="${INPUT_RDS}"
fi

binds="${project_root}:/work,/ibmgpfs:/ibmgpfs"

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] DATASET=${DATASET}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] K_VALUES=${K_VALUES}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"

action=(
  Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/V_check_k/V1st_check_k.R
  --project_root=/work
  --dataset="${DATASET}"
  --input_rds="${INPUT_RDS_CONTAINER}"
  --k_values="${K_VALUES}"
)
if [ -n "${N_EDGES}" ]; then
  action+=(--n_edges="${N_EDGES}")
fi

singularity exec --cleanenv \
  --bind "${binds}" \
  --pwd /work \
  "${sif}" "${action[@]}"
