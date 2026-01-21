#!/bin/bash
#SBATCH -J build_scdev_r41_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 8
#SBATCH -o outputs/logs/containers/build_scdev_r41_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
def_file="${project_root}/containers/scdevelopment_r41.def"
# Do NOT overwrite the existing SIF used by running pipelines.
# Generate a new image name per build; you can override by setting SIF_TAG.
tag="${SIF_TAG:-${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}}"
sif_out="${project_root}/outputs/containers/scdevelopment_r41_${tag}.sif"

mkdir -p "${project_root}/outputs/logs/containers"
mkdir -p "${project_root}/outputs/containers"
cd "${project_root}"

# Network: compute nodes require proxy to access internet (docker hub / CRAN / apt).
export http_proxy=http://10.11.100.5:3128
export HTTP_PROXY=http://10.11.100.5:3128
export https_proxy=http://10.11.100.5:3128
export HTTPS_PROXY=http://10.11.100.5:3128
export ftp_proxy=http://10.11.100.5:3128
export FTP_PROXY=http://10.11.100.5:3128
export all_proxy=http://10.11.100.5:3128
export ALL_PROXY=http://10.11.100.5:3128

# Also pass through to singularity runtime/build stages when applicable.
export SINGULARITYENV_http_proxy="${http_proxy}"
export SINGULARITYENV_HTTP_PROXY="${HTTP_PROXY}"
export SINGULARITYENV_https_proxy="${https_proxy}"
export SINGULARITYENV_HTTPS_PROXY="${HTTPS_PROXY}"
export SINGULARITYENV_ftp_proxy="${ftp_proxy}"
export SINGULARITYENV_FTP_PROXY="${FTP_PROXY}"
export SINGULARITYENV_all_proxy="${all_proxy}"
export SINGULARITYENV_ALL_PROXY="${ALL_PROXY}"

module load singularity/3.7.0
singularity --version

echo "[INFO] def_file=${def_file}"
echo "[INFO] sif_out=${sif_out}"
echo "[INFO] http_proxy=${http_proxy}"
echo "[INFO] https_proxy=${https_proxy}"

# This cluster does not grant root build privileges and fakeroot may be
# unavailable (no /etc/subuid mapping). Use remote build.
if ! singularity remote status 2>/dev/null | grep -qi "logged in"; then
  if [ -n "${SINGULARITY_REMOTE_TOKEN:-}" ]; then
    echo "[INFO] Logging into singularity remote (token via SINGULARITY_REMOTE_TOKEN)"
    token_file="$(mktemp)"
    printf "%s" "${SINGULARITY_REMOTE_TOKEN}" > "${token_file}"
    singularity remote login --tokenfile "${token_file}"
    rm -f "${token_file}"
  fi
fi

echo "[INFO] Building with: singularity build --remote ${sif_out} ${def_file}"
singularity build --remote "${sif_out}" "${def_file}"

echo "[INFO] Build done: ${sif_out}"
