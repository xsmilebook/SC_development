#!/usr/bin/env bash
set -euo pipefail

OUTPUT_RDS=${1:?"output rds required"}
TEST_N=${2:-0}

set --
source /GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/bin/activate
conda activate ML

python combat_gam/scripts/run_combat_gam_neuroharmonize.py \
  --input-rds /ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds \
  --output-rds "${OUTPUT_RDS}" \
  --batch-col site \
  --id-col subID \
  --age-col age \
  --sex-col sex \
  --mean-fd-col mean_fd \
  --test-n "${TEST_N}"
