#!/bin/bash
tck_path=$1
tck_name=$2
output_dir=$3
tckresample -endpoints ${tck_path}/${tck_name}.tck ${output_dir}/${tck_name}_endpoints.tck
tckconvert ${output_dir}/${tck_name}_endpoints.tck ${output_dir}/${tck_name}_endpoints.vtk
