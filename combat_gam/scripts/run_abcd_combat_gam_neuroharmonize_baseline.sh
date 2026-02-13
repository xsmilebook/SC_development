#!/usr/bin/env bash
set -euo pipefail

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds}
OUTPUT_DIR=${2:-outputs/results/combat_gam/abcd/baseline}
TEST_N=${3:-0}

set --
source /GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/bin/activate
conda activate scdevelopment

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant cognition \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant pfactor \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant baseline_age_sex_meanfd \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant cbcl \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant comp_agecorrected_baseline \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"

python combat_gam/scripts/run_abcd_combat_gam_neuroharmonize.py \
  --input-rds "${INPUT_RDS}" \
  --output-dir "${OUTPUT_DIR}" \
  --variant fluidcomp_fc_baseline \
  --id-col subID \
  --scan-col scanID \
  --batch-col siteID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --event-col eventname \
  --test-n "${TEST_N}"
