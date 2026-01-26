#!/bin/bash
#SBATCH -J abcd_devmodel_cg_cv75_sif
#SBATCH -p q_fat_c,q_fat,q_fat_l
#SBATCH -c 50
#SBATCH -o outputs/logs/abcd_devmodel_combatgam_CV75_container_%j.log

set -euo pipefail

project_root="/ibmgpfs/cuizaixu_lab/xuhaoshu/projects/SCDevelopment"
# Default to the "stable" SIF path; allow overriding for testing a newly built image.
sif="${SIF_PATH:-${project_root}/outputs/containers/scdevelopment_r41.sif}"

CVTHR="${CVTHR:-75}"
INPUT_RDS="${INPUT_RDS:-${project_root}/outputs/results/combat_gam/abcd/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds}"
EUCLID_CSV="${EUCLID_CSV:-${project_root}/wd/interdataFolder_ABCD/average_EuclideanDistance_12.csv}"
N_EDGES="${N_EDGES:-}"
SKIP_POSTERIOR="${SKIP_POSTERIOR:-0}"
MAKE_MATRIX_GRAPHS="${MAKE_MATRIX_GRAPHS:-0}"

mkdir -p "${project_root}/outputs/logs"
cd "${project_root}"

module load singularity/3.7.0
singularity --version

if [ ! -f "${sif}" ]; then
  echo "[ERROR] Missing SIF: ${sif}"
  echo "[HINT] Build it first: sbatch sbatch/build_scdevelopment_r41_container.sbatch"
  echo "[HINT] Or pass a versioned image via: SIF_PATH=/path/to/scdevelopment_r41_<tag>.sif sbatch sbatch/run_abcd_devmodel_combatgam_CV75_container.sbatch"
  exit 1
fi

echo "[INFO] project_root=${project_root}"
echo "[INFO] sif=${sif}"
echo "[INFO] CVTHR=${CVTHR}"
echo "[INFO] INPUT_RDS=${INPUT_RDS}"
echo "[INFO] EUCLID_CSV=${EUCLID_CSV}"
echo "[INFO] N_EDGES=${N_EDGES:-<unset>}"
echo "[INFO] SKIP_POSTERIOR=${SKIP_POSTERIOR}"
echo "[INFO] MAKE_MATRIX_GRAPHS=${MAKE_MATRIX_GRAPHS}"

FORCE_S3="${FORCE_S3:-0}"
FORCE_S4="${FORCE_S4:-0}"

decile_tiff="${project_root}/outputs/figures/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV${CVTHR}/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype10.tiff"
if [ ! -f "${decile_tiff}" ]; then
  FORCE_S3=1
fi

corr_tiff="${project_root}/outputs/figures/2nd_fitdevelopmentalmodel/abcd/combat_gam/CV${CVTHR}/correlation_sumSCinvnode_SCrank/meanpartialRsq2_SCrankcorr_n12.tiff"
if [ ! -f "${corr_tiff}" ]; then
  FORCE_S4=1
fi

echo "[INFO] FORCE_S3=${FORCE_S3} (set FORCE_S3=0/1 to override)"
echo "[INFO] FORCE_S4=${FORCE_S4} (set FORCE_S4=0/1 to override)"

exec_sif() {
  singularity exec --cleanenv \
    --bind "${project_root}:/work" \
    "${sif}" "$@"
}

exec_sif Rscript --version

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S1st_fitgammodels_SA_ds_sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --input_rds="/work/${INPUT_RDS#${project_root}/}" \
  ${N_EDGES:+--n_edges=${N_EDGES}}

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S2nd_calculatederivative_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  ${N_EDGES:+--n_edges=${N_EDGES}} \
  --skip_posterior="${SKIP_POSTERIOR}"

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S3rd_visualizationfitSCcurves_SA12sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --force="${FORCE_S3}"

exec_sif Rscript --vanilla /work/development_script/2nd_fitdevelopmentalmodel/S4th_correlationTo_SArank_SA12sumSCinvnode_ABCD_combatgam.R \
  --project_root=/work \
  --cvthr="${CVTHR}" \
  --input_rds="/work/${INPUT_RDS#${project_root}/}" \
  --euclid_csv="/work/${EUCLID_CSV#${project_root}/}" \
  --make_matrix_graphs="${MAKE_MATRIX_GRAPHS}" \
  --force="${FORCE_S4}"
