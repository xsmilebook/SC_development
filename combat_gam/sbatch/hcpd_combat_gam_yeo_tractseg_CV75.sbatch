#!/usr/bin/env bash
#SBATCH -J combat_hcpd_gam_yeo_trk
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH --cpus-per-task=8
#SBATCH --output=outputs/logs/combat_gam/combat_hcpd_gam_yeo_tractseg_CV75_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
conda_root="/GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3"

# Optional: small-sample smoke test (0 means full data)
TEST_N="${TEST_N:-0}"

# Column mapping (override if needed)
BATCH_COL="${BATCH_COL:-site}"
ID_COL="${ID_COL:-subID}"
AGE_COL="${AGE_COL:-age}"
SEX_COL="${SEX_COL:-sex}"
MEAN_FD_COL="${MEAN_FD_COL:-mean_fd}"

out_dir="${OUT_DIR:-${project_root}/outputs/results/combat_gam/hcpd}"
log_dir="${project_root}/outputs/logs/combat_gam"

# Some derived HCPD tables (e.g., Yeo/TractSeg) may miss the `site` column.
# If so, merge `site` from the canonical SA12 merge table by `subID`.
DEFAULT_BATCH_SOURCE_RDS="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
BATCH_SOURCE_RDS="${BATCH_SOURCE_RDS:-${DEFAULT_BATCH_SOURCE_RDS}}"

mkdir -p "${out_dir}" "${log_dir}"
cd "${project_root}"

source "${conda_root}/etc/profile.d/conda.sh"
conda activate scdevelopment

export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[INFO] project_root=${project_root}"
echo "[INFO] out_dir=${out_dir}"
echo "[INFO] TEST_N=${TEST_N}"
echo "[INFO] columns: batch=${BATCH_COL} id=${ID_COL} age=${AGE_COL} sex=${SEX_COL} mean_fd=${MEAN_FD_COL}"
python --version || true

run_one() {
  local input_rds="$1"
  local output_rds="$2"
  local tag="$3"

  if [ ! -f "${input_rds}" ]; then
    echo "[ERROR] Missing input: ${input_rds}"
    exit 1
  fi

  echo "[INFO] === ${tag} ==="
  echo "[INFO] input_rds=${input_rds}"
  echo "[INFO] output_rds=${output_rds}"

  python combat_gam/scripts/run_combat_gam_neuroharmonize.py \
    --input-rds "${input_rds}" \
    --output-rds "${output_rds}" \
    --batch-col "${BATCH_COL}" \
    --batch-source-rds "${BATCH_SOURCE_RDS}" \
    --batch-source-id-col "${ID_COL}" \
    --batch-source-batch-col "${BATCH_COL}" \
    --id-col "${ID_COL}" \
    --age-col "${AGE_COL}" \
    --sex-col "${SEX_COL}" \
    --mean-fd-col "${MEAN_FD_COL}" \
    --test-n "${TEST_N}"

  if [ ! -f "${output_rds}" ]; then
    echo "[ERROR] Missing output after run: ${output_rds}"
    exit 1
  fi
  echo "[INFO] done: ${tag}"
}

# Inputs (absolute paths as provided)
in_yeo7="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
in_yeo17="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"
in_tractseg="/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.TractSeg.merge.rds"

# Outputs (project-local)
out_yeo7="${out_dir}/SCdata_Yeo7_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
out_yeo17="${out_dir}/SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"
out_tractseg="${out_dir}/SCdata_SA12_CV75_sumSCinvnode.TractSeg.combatgam.rds"

run_one "${in_yeo7}" "${out_yeo7}" "HCPD_Yeo7"
run_one "${in_yeo17}" "${out_yeo17}" "HCPD_Yeo17"
run_one "${in_tractseg}" "${out_tractseg}" "HCPD_TractSeg_major_bundle"

echo "[INFO] All done."
