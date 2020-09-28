#!/bin/bash
num_streamlines = $1
mkdir /home/neuro/output/mrtrixqsi
cd /home/neuro/output/mrtrixqsi
tckconvert ${num_streamlines}_endpoints.tck ${num_streamlines}_endpoints.vtk
