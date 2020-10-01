#!/bin/bash
tck_path=$1
tck_name=$2
tckconvert ${tck_path}/${tck_name}.tck /home/neuro/output/endpoints/${tck_name}_endpoints.vtk
