#!/bin/bash
#SBATCH -J plot_abcd_r2
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=72
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --output=outputs/logs/combat_gam/plot_abcd_variance_decomposition_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

mkdir -p "${project_root}/outputs/logs/combat_gam"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch combat_gam/sbatch/plot_abcd_variance_decomposition.sbatch"
  exit 1
fi

export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work" \
    --bind /ibmgpfs:/ibmgpfs \
    --bind /GPFS:/GPFS \
    "${sif}" "$@"
}

exec_sif Rscript --version
exec_sif Rscript --vanilla /work/combat_gam/scripts/plot_abcd_variance_decomposition.R
