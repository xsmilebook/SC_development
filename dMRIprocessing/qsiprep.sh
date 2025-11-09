#!/bin/bash
#SBATCH -J HCPD-QSIPREP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu 20000
#SBATCH -p q_cn

#this script is only for preprocess

module load singularity/3.7.0
module load freesurfer

#User inputs:
bids_root_dir=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/BIDS
bids_root_dir_output=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/qsiprep
bids_root_dir_output_wd4singularity=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/wd
freesurfer_fatherfold=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/fmriresults01
templateflow=/ibmgpfs/cuizaixu_lab/xuxiaoyu/softwarepackages/templateflow
subj=$1
nthreads=8

freesurfer_dir=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/fmriresults01/${subj:4:14}_V1_MR/T1w
#Run qsiprep
echo ""
echo "Running qsiprep on participant: $subj"
echo ""

#Make qsiprep directory and participant directory in derivatives folder

if [ ! -d $bids_root_dir_output/${subj} ]; then
    mkdir $bids_root_dir_output/${subj}
fi

if [ ! -d $bids_root_dir_output_wd4singularity/qsiprep/${subj} ]; then
    mkdir $bids_root_dir_output_wd4singularity/qsiprep/${subj}
fi

#Run qsiprep_prep
export SINGULARITYENV_TEMPLATEFLOW_HOME='/ibmgpfs/cuizaixu_lab/xuxiaoyu/softwarepackages/templateflow'
export TMPDIR='/ibmgpfs/cuizaixu_lab/xuxiaoyu/tmp'
unset PYTHONPATH; singularity run --cleanenv --bind $bids_root_dir \
    -B $bids_root_dir_output_wd4singularity/qsiprep/${subj}:/wd \
    -B $bids_root_dir:/inputbids \
    -B $bids_root_dir_output/${subj}:/output \
    -B $bids_root_dir_output:/recon_input \
    -B $freesurfer_dir:/freesurfer \
    -B $templateflow:/ibmgpfs/cuizaixu_lab/xuxiaoyu/softwarepackages/templateflow \
    -B /ibmgpfs/cuizaixu_lab/xuxiaoyu/tmp:/tmp \
    -B /ibmgpfs/cuizaixu_lab/xuxiaoyu/freesurfer_license:/freesurfer_license \
    /ibmgpfs/cuizaixu_lab/xuxiaoyu/softwarepackages/qsiprep.sif \
    /inputbids /output \
    participant \
    --participant_label ${subj} \
    --unringing-method mrdegibbs \
    --output-resolution 1.5 \
    --distortion-group-merge average \
    --recon_input /recon_input \
    --recon_spec mrtrix_multishell_msmt_ACT-hsvs \
    --freesurfer-input /freesurfer \
    --skip-bids-validation \
    -w /wd \
    --verbose \
    --notrack \
    --nthreads $nthreads \
    --mem-mb 32000 \
    --fs-license-file /freesurfer_license/license.txt
