#!/bin/bash

dMRIdatapath=$1
bvals=$2
bvecs=$3
T1path=$4
brainmask=$5
intermediary_output_path=$6
num_streamlines=$7

mkdir ${intermediary_output_path}
cd ${intermediary_output_path}

#3dresample -rmode NN -prefix resampled_mask.nii -input } 
#AFNI resample anatomical

5ttgen fsl ${T1path} 5TT.mif -premasked

mrconvert ${dMRIdatapath} DWI.mif -fslgrad ${bvecs} ${bvals} -datatype float32 -stride 0,0,0,1

dwi2response dhollander DWI.mif wm_txt.txt gm_txt.txt csf_txt.txt -force
#txt files are outputs

dwi2fod msmt_csd DWI.mif wm_txt.txt wm_fod.mif gm_txt.txt gm_fod.mif csf_txt.txt csf_fod.mif -lmax 4,8,8 -mask ${brainmask}
#takes outputs of 

mtnormalise wm_fod.mif wm_fod_norm.mif gm_fod.mif gm_fod_norm.mif csf_fod.mif csf_fod_norm.mif -mask ${brainmask}
#qsiprep step

tckgen wm_fod.mif ${num_streamlines}.tck -act 5TT.mif -backtrack -crop_at_gmwmi -seed_dynamic wm_fod.mif -maxlength 250 -minlength 30 -select ${num_streamlines} -power 0.33
#input and output are first two arguments

#tcksift2 ${num_streamlines}.tck wm_fod.mif -act 5TT.mif -out_mu out_mu.txt -out_weights sift_weights.txt
#are all the other options part of default qsi?

tckresample -endpoints ${num_streamlines}.tck ${num_streamlines}_endpoints.tck

tckconvert ${num_streamlines}_endpoints.tck ${num_streamlines}_endpoints.vtk