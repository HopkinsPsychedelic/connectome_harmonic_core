#!/bin/bash

#bash transform_qsi_tck.sh ~/native_t1_brain.nii.gz ~/sub-105923_desc-preproc_T1w_ras.nii.gz /data/HCP_Raw/derivatives/qsirecon/sub-105923/dwi/sub-105923_space-T1w_desc-preproc_space-T1w_desc-tracks_ifod2.tck

native_t1=$1
acpc_template=$2
tracks=$3

#register native t1 to acpc with antsRegistration
antsRegistration --collapse-output-transforms 1 --dimensionality 3 --initial-moving-transform [ ${acpc_template}, ${native_t1}, 1 ] --initialize-transforms-per-stage 0 --interpolation LanczosWindowedSinc --output [ transform, transform_Warped.nii.gz ] --transform Rigid[ 0.2 ] --metric Mattes[ ${acpc_template}, ${native_t1}, 1, 32, Random, 0.25 ] --convergence [ 10000x1000x10000x10000, 1e-06, 10 ] --smoothing-sigmas 7.0x3.0x1.0x0.0vox --shrink-factors 8x4x2x1 --use-histogram-matching 1 --winsorize-image-intensities [ 0.025, 0.975 ]  --write-composite-transform 0

#convert ANTs affine matrix (native t1 to acpc) to ITK txt format, then to mrtrix format
ConvertTransformFile 3 transform0GenericAffine.mat transform.txt
transformconvert transform.txt itk_import transform_mrtrix.txt

#convert native t1 (moving) and acpc template (reference) images to mrtrix format
mrconvert ${native_t1} native_t1.mif
mrconvert ${acpc_template} acpc_template.mif

#initialize identity warp for native t1 (moving)
warpinit native_t1.mif id_warp.mif

#transform identity warp to image warp with t1->acpc matrix
mrtransform id_warp.mif -linear transform_mrtrix.txt native_to_acpc_warp.mif -template acpc_template.mif -interp linear -nan

#check that the following two are equivalent (they are. both went from native to acpc space)
#mrtransform native_t1.mif -warp native_to_acpc_warp.mif warped_t1_warp.mif
#mrtransform native_t1.mif -linear transform_mrtrix.txt warped_t1_linear.mif -template acpc_template.mif
#mrconvert warped_t1_warp.mif warped_t1_warp.nii
#mrconvert warped_t1_linear.mif warped_t1_linear.nii

#transform tracks from acpc space to native space (INVERSE OF ABOVE TRANSFORM)
tcktransform ${tracks} native_to_acpc_warp.mif transformed_tracks.tck

tckresample -downsample 10 transformed_tracks.tck ds_tf_tracks_chap-warp.tck

#tckresample -endpoints ds_tf_tracks.tck ds_tf.tck

#tckconvert ds_tf_endp.tck ds_tf_endp.vtk

