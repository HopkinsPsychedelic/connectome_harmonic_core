#!/bin/bash

dMRIdatapath=$1
bvals=$2
bvecs=$3
freesurfer_out=$4
brainmask=$5
intermediary_output_path=$6
num_streamlines=$7

mkdir ${intermediary_output_path}
cd ${intermediary_output_path}

5ttgen hsvs ${freesurfer_out} 5TT.mif 

mrconvert ${dMRIdatapath} DWI.mif -fslgrad ${bvecs} ${bvals} -datatype float32 -stride 0,0,0,1

mrconvert ${brainmask} brainmask.mif -datatype bit

dwi2response dhollander DWI.mif RF_WM.txt RF_GM.txt RF_CSF.txt -mask brainmask.mif

dwi2fod msmt_csd DWI.mif RF_WM.txt WM_FODs.mif RF_GM.txt GM.mif RF_CSF.txt CSF.mif -lmax 10,0,0  

mrconvert WM_FODs.mif - -coord 3 0 | mrcat FOD_CSF.mif FOD_GM.mif - tissues.mif -axis 3

tckgen WM_FODs.mif ${num_streamlines}.tck -act 5TT.mif -backtrack -crop_at_gmwmi -seed_dynamic WM_FODs.mif -maxlength 250 -select ${num_streamlines} -power 0.33

tckresample -endpoints ${num_streamlines}.tck ${num_streamlines}_endpoints.tck

tckconvert ${num_streamlines}_endpoints.tck ${num_streamlines}_endpoints.vtk
