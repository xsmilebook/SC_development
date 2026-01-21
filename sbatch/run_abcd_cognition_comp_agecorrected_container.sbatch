#!/bin/bash
#SBATCH -J abcd_cog_compac
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 72
#SBATCH --mem=180G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/abcd_cognition_comp_agecorrected_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

# External ABCD reference file used for Euclidean-distance control (existing pipeline path).
ABCD_EUCLID_CSV="${ABCD_EUCLID_CSV:-/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/average_EuclideanDistance_12.csv}"

# Force recomputation (0=use cached results when available).
FORCE="${FORCE:-0}"

mkdir -p "${project_root}/outputs/logs"
mkdir -p "${project_root}/outputs/reports/5th_cognition/abcd/comp_agecorrected"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_abcd_cognition_comp_agecorrected_container.sbatch"
  exit 1
fi

if [ ! -f "${ABCD_EUCLID_CSV}" ]; then
  echo "[ERROR] Missing ABCD_EUCLID_CSV: ${ABCD_EUCLID_CSV}"
  exit 1
fi

# Avoid BLAS oversubscription with parallel clusters.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Pass required env vars into the container even under --cleanenv.
export SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-72}"
export SINGULARITYENV_ABCD_EUCLID_CSV="${ABCD_EUCLID_CSV}"
export SINGULARITYENV_FORCE="${FORCE}"
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

singularity exec --cleanenv \
  --bind "${project_root}:/work" \
  --pwd /work \
  "${sif}" \
  bash -lc 'Rscript --vanilla - << "RS"
out_dir <- "/work/outputs/reports/5th_cognition/abcd/comp_agecorrected"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  "/work/development_script/5th_cognition/S1st_testCog_correlation_wholesample_ABCD.Rmd",
  output_dir = out_dir,
  quiet = TRUE
)
rmarkdown::render(
  "/work/development_script/5th_cognition/S2nd_compositescorePlot_scatterplot_ABCD.Rmd",
  output_dir = out_dir,
  quiet = TRUE
)
RS'
