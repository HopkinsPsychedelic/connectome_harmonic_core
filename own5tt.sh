#!/bin/bash

dMRIdatapath=$1
bvals=$2
bvecs=$3
T1path=$4
fivett=$5
intermediary_output_path=$6
num_streamlines=$7

mkdir ${intermediary_output_path}
cd ${intermediary_output_path}

mrconvert ${dMRIdatapath} DWI.mif -fslgrad ${bvecs} ${bvals} -datatype float32 -stride 0,0,0,1

dwi2response msmt_5tt DWI.mif ${fivett} RF_WM.txt RF_GM.txt RF_CSF.txt -voxels RF_voxels.mif

dwi2fod msmt_csd DWI.mif RF_WM.txt WM_FODs.mif RF_GM.txt GM.mif RF_CSF.txt CSF.mif -mask ${brainmask}

tckgen WM_FODs.mif ${num_streamlines}.tck -act ${fivett} -backtrack -crop_at_gmwmi -seed_dynamic WM_FODs.mif -maxlength 250 -select ${num_streamlines} -power 0.33

tckresample -endpoints ${num_streamlines}.tck ${num_streamlines}_endpoints.tck

tckconvert ${num_streamlines}_endpoints.tck ${num_streamlines}_endpoints.vtk