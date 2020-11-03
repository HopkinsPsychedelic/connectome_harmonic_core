#!/bin/bash
tck_path=$1
tck_name=$2
output_dir=$3
tckresample -endpoints ${tck_path}/${tck_name}.tck ${tck_path}/${tck_name}_endpoints.tck
tckconvert ${tck_path}/${tck_name}_endpoints.tck ${output_dir}/${tck_name}.vtk
