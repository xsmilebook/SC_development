#!/bin/bash
#SBATCH -J cbcl_assoc_full
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH --mem=120G
#SBATCH -t 72:00:00
#SBATCH -o outputs/logs/cbcl_assoc_full_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"
abcd_interdata="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_cbcl_assoc_full.sbatch"
  exit 1
fi

binds="${project_root}:${project_root}"
if [ -d "${abcd_interdata}" ]; then
  binds="${binds},${abcd_interdata}:${abcd_interdata}"
else
  echo "[WARN] Missing ABCD interdata folder (may be OK if not needed): ${abcd_interdata}"
fi

singularity exec --cleanenv \
  --bind "${binds}" \
  --pwd "${project_root}" \
  "${sif}" \
  Rscript --vanilla "${project_root}/development_script/6th_pfactor/S2nd_cbcl_totalraw_effect_continuous_ABCD.R"
