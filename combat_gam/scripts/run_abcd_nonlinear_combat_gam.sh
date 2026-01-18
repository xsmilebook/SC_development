#!/usr/bin/env bash
set -euo pipefail

INPUT_RDS=${1:-/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.merge.rds}
OUTPUT_DIR=${2:?"output dir required"}
TEST_N=${3:-0}

set --
source /GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/bin/activate
conda activate scdevelopment

export R_LIBS_USER=
export R_LIBS=

Rscript combat_gam/scripts/run_abcd_nonlinear_combat_gam.R "${INPUT_RDS}" "${OUTPUT_DIR}" "${TEST_N}"
