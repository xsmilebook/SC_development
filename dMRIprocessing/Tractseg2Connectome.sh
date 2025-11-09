#!/bin/bash
#SBATCH --nodes=1 # OpenMP requires a single node
#SBATCH -p q_cn
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=4
#SBATCH --mail-user xuxiaoyu@cibr.ac.nc
##### END OF JOB DEFINITION  #####

module load singularity
module load mrtrix3
subj=$1

processedpath=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/qsiprep
freesurfermri=/ibmgpfs/cuizaixu_lab/xuxiaoyu/HCPD/processed/fmriresults01/${subj:4:14}_V1_MR/T1w/${subj}/mri
echo $subj

#export http_proxy=10.11.100.5:3128
#export HTTP_PROXY=10.11.100.5:3128
#export https_proxy=10.11.100.5:3128
#export HTTPS_PROXY=10.11.100.5:3128
#export ftp_proxy=10.11.100.5:3128
#export FTP_PROXY=10.11.100.5:3128
#export all_proxy=10.11.100.5:3128
#export ALL_PROXY=10.11.100.5:3128

# 1. Convert the orientation of FOD to the Tractseg required orientation
# LPS to LAS
mrconvert ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd.mif.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.mif.gz -stride -1,2,3,4 -force

# 2. Conert .mif to .nii
mrconvert ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.mif.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.nii.gz -force

# 3. Compute the peak image
sh2peaks ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.nii.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2_peak.nii.gz -force

# 4. Prepare T1
mrconvert ${processedpath}/${subj}/qsiprep/${subj}/anat/${subj}_desc-preproc_brain.nii.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/T1w_acpc_dc_restore_brain.nii.gz -stride -1,2,3 -force

# 5. Segment bundle start and end regions
## TractSeg tracts
TractSeg -i ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2_peak.nii.gz \
    -o ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output --output_type tract_segmentation

TractSeg -i ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2_peak.nii.gz \
    -o ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output --output_type endings_segmentation

# 6. Tracking the 72 tracts
Tracking -i ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.nii.gz \
    -o ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output --track_FODs iFOD2 --nr_fibers 10000

# 7. Merge all the tck to build a global tractography
tckedit ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/FOD_iFOD2_trackings/*.tck ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/All_tracks720k.tck

# 8. Generate sift2 weights
mrconvert ${processedpath}/${subj}/qsirecon/${subj}/anat/${subj}_desc-preproc_space-fsnative_desc-hsvs_5tt.mif ${processedpath}/${subj}/qsirecon/${subj}/anat/${subj}_desc-preproc_space-fsnative_desc-hsvs_5tt_LAS.mif -stride -1,2,3,4 -force

tcksift2 -act ${processedpath}/${subj}/qsirecon/${subj}/anat/${subj}_desc-preproc_space-fsnative_desc-hsvs_5tt_LAS.mif -out_mu ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/sift_mu720k.txt -out_coeffs ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/sift_coeffs720k.txt ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/All_tracks720k.tck ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-wmFODmtnormed_msmtcsd_LAS2.mif.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/sift_720k.txt

# 9. Reconstruct connectome
mrconvert ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas_LAS.nii.gz -stride -1,2,3,4 -force

# 1) inverse node volume
tck2connectome -force -symmetric -nthreads 72 -assignment_radial_search 2 -scale_invnodevol \
    -tck_weights_in ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/sift_720k.txt \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/All_tracks720k.tck \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas_LAS.nii.gz \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/Tractseg_schaefer400_SC_invnode.csv \
    -out_assignment ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/assignments_720k_nodeS400.csv
# 2) without inverse node volume
tck2connectome -force -symmetric -nthreads 72 -assignment_radial_search 2 \
    -tck_weights_in ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/sift_720k.txt \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/All_tracks720k.tck \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas_LAS.nii.gz \
    ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/Tractseg_schaefer400_SC.csv \
    -out_assignment ${processedpath}/${subj}/qsirecon/${subj}/dwi/tractseg_output/assignments_720k_nodeS400.csv


