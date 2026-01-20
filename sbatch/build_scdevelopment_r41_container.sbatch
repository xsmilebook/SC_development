#!/bin/bash
#SBATCH -J build_scdev_r41_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o outputs/logs/containers/build_scdev_r41_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
def_file="${project_root}/containers/scdevelopment_r41.def"
sif_out="${project_root}/outputs/containers/scdevelopment_r41.sif"

mkdir -p "${project_root}/outputs/logs/containers"
mkdir -p "${project_root}/outputs/containers"
cd "${project_root}"

# Cluster usage: module load singularity/3.7.0
if [ -f /usr/share/Modules/init/bash ]; then
  # shellcheck disable=SC1091
  source /usr/share/Modules/init/bash
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi

module load singularity/3.7.0
singularity --version

echo "[INFO] def_file=${def_file}"
echo "[INFO] sif_out=${sif_out}"

singularity build --fakeroot "${sif_out}" "${def_file}"

