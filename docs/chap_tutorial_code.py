#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 12:10:03 2021

@author: bwinston
"""

singularity run --cleanenv -B /data/hcp900:/data/hcp900 -B /bwinston:/bwinston fmriprep.sif \
    /data/hcp900/raw_data /bwinston/derivatives participant \
    --work /bwinston/derivatives/fmriprep_103818_work \
    --fs-license-file /bwinston/license.txt \
    --output-spaces anat \
    --participant-label 103818 
    
singularity run --cleanenv -B /data/hcp900:/data/hcp900 -B /bwinston:/bwinston mrtrix_connectome.simg \ 
    /data/hcp900/raw_data /bwinston/derivatives preproc \
    --scratch /bwinston/derivatives/mrtrix_103818_work \
    --t1w_preproc /bwinston/derivatives/fmriprep/sub-103818/anat/sub-103818_desc-preproc_T1w.nii.gz \ 
    --participant_label 103818
    
sudo docker run -it --rm \
    -v /data/hcp_test_retest_pp/:/data/hcp_test_retest_pp/  \
    -p 8888:8888 \
    winstonian3/connectome_harmonic:latest \
    /data/hcp_test_retest_pp/derivatives participant \
    --hcp_dir /data/hcp_test_retest_pp/source_data \
    --participant_label 105923

singularity run --cleanenv -B /data/hcp_test_retest_pp/:/data/hcp_test_retest_pp/ chap.sif \ 
    /data/hcp_test_retest_pp/derivatives participant \
    --hcp_dir /data/hcp_test_retest_pp/source_data \
    --participant_label 105923
    
sudo docker run -it --rm \
    -v /data/HCP_Raw:/data/HCP_Raw/ \
    -p 8888:8888 \
    winstonian3/connectome_harmonic:latest \
    /data/HCP_Raw/derivatives participant \
    --mrtrix_dir /data/HCP_Raw/derivatives/MRtrix3_connectome-preproc \
    --ciftify_dir /data/HCP_Raw/derivatives/ciftify \
    --freesurfer_dir /data/HCP_Raw/derivatives/freesurfer \
    --participant_label 105923

singularity run --cleanenv -B /data/HCP_Raw:/data/HCP_Raw/ chap.sif \
    /data/HCP_Raw/derivatives participant \
    --mrtrix_dir /data/HCP_Raw/derivatives/MRtrix3_connectome-preproc \
    --ciftify_dir /data/HCP_Raw/derivatives/ciftify \
    --freesurfer_dir /data/HCP_Raw/derivatives/freesurfer \
    --participant_label 105923
