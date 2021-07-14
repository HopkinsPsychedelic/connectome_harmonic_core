#!/bin/bash

#bash transform_qsi_tck.sh /data/HCP_Raw/derivatives/tmp_105923/qsiprep_wf/single_subject_105923_wf/anat_preproc_wf/skullstrip_wf/rigid_acpc_align/transform0GenericAffine.mat /data/HCP_Raw/derivatives/tmp_105923/qsiprep_wf/single_subject_105923_wf/anat_preproc_wf/skullstrip_wf/t1_skull_strip/highres001_N4Corrected0.nii.gz /usr/local/connectome_harmonic_core/connectome_harmonic_core/mni_1mm_t1w_lps_brain.nii.gz /data/HCP_Raw/derivatives/qsirecon/sub-105923/dwi/sub-105923_space-T1w_desc-preproc_space-T1w_desc-tracks_ifod2.tck

#qsiprep generates native t1 to acpc transform matrix. 
#we want to get streamlines from acpc back to native

transform=$1
native_t1=$2 
acpc_template=$3
tracks=$4

#convert ANTs affine matrix (native T1 to acpc) to ITK txt format
ConvertTransformFile 3 ${transform} transform.txt

#get transform in mrtrix format
transformconvert transform.txt itk_import transform_mrtrix.txt

#get inverse of above (acpc to T1)
transformcalc transform_mrtrix.txt invert transform_mrtrix_inv.txt

#convert acpc (moving) and native_t1 (reference) to mrtrix format
mrconvert ${native_t1} native_t1.mif
mrconvert ${acpc_template} acpc_template.mif

#initialize identity warp for t1 
warpinit acpc_template.mif id_warp.mif

#transform identity warp to image warp with transform matrix
transformcompose id_warp.mif transform_mrtrix_inv.txt warp_streamlines_to_native.mif -template native_t1.mif

#check that the following two are equivalent (they are. both went from acpc to native space)
#mrtransform acpc_template.mif -warp acpc_to_native_warp.mif warped_acpc_warp.mif
#mrtransform acpc_template.mif -linear transform_mrtrix_inv.txt warped_acpc_linear.mif -template native_t1.mif



#transform tracks from acpc space to native space (INVERSE OF ABOVE TRANSFORM)
tcktransform ${tracks} warp_streamlines_to_native.mif transformed_tracks.tck

tckresample -downsample 15 transformed_tracks.tck ds_tf_tracks.tck

tckresample -endpoints ds_tf_tracks.tck ds_tf_endp.tck

tckconvert ds_tf_endp.tck ds_tf_endp.vtk

    