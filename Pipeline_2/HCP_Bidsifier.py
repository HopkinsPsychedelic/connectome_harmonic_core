#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 17:19:49 2020

@author: bwinsto2

very lazy code to bidsify HCP unprocessed data - pls don't judge
acq messed up and BOLD -- don't use without fixing those things

make it create BIDS folders 
add sessions for test and retest

"""

#unzip and put into sub/func
from zipfile import ZipFile
import os, shutil
source_root = '/data2/HCP_Raw/source_data'
for sub in os.listdir(source_root):
    for ses in os.listdir(f'{source_root/{sub}')
        func = f'/data2/HCP_Raw/raw_data/sub-{sub}/func'
        for zipdir in os.listdir(f'{source_root}/{str(sub)}'):
            if 'md5' not in zipdir:
                for dtype in ['tfMRI_WM','rfMRI_REST2','rfMRI_REST1']:
                    if dtype in zipdir:              
                        with ZipFile(f'{source_root}/{sub}/{sub}_3T_{dtype}_unproc.zip', 'r') as zipObj:
                            zipObj.extractall(func)
                            for direction in ['LR', 'RL']:
                                shutil.copy(f'{func}/{sub}/unprocessed/3T/{dtype}_{direction}/{sub}_3T_{dtype}_{direction}.nii.gz',f'/data2/HCP_Raw/raw_data/sub-{sub}/func')
                        shutil.rmtree(f'/data2/HCP_Raw/raw_data/sub-{sub}/func/{sub}')


for sub in os.listdir(source_root):
    anat = f'/data2/HCP_Raw/raw_data/sub-{sub}/anat'
    for zipdir in os.listdir(f'{source_root}/{str(sub)}'):
        if 'md5' not in zipdir:
            if 'Struc' in zipdir:              
                with ZipFile(f'{source_root}/{sub}/{sub}_3T_Structural_unproc.zip', 'r') as zipObj:
                    zipObj.extractall(anat)
                    for direction in ['T1w_MPR1', 'T2w_SPC1']:
                        shutil.copy(f'{anat}/{sub}/unprocessed/3T/{direction}/{sub}_3T_{direction}.nii.gz',f'/data2/HCP_Raw/raw_data/sub-{sub}/anat')
                shutil.rmtree(f'/data2/HCP_Raw/raw_data/sub-{sub}/anat/{sub}')


for sub in os.listdir(source_root):
    diff = f'/data2/HCP_Raw/raw_data/sub-{sub}/dwi'
    for zipdir in os.listdir(f'{source_root}/{str(sub)}'):
        if 'md5' not in zipdir:
            if 'Diff' in zipdir:              
                with ZipFile(f'{source_root}/{sub}/{sub}_3T_Diffusion_unproc.zip', 'r') as zipObj:
                    zipObj.extractall(diff)
                    for direction in ['DWI_dir95_LR', 'DWI_dir95_RL', 'DWI_dir96_LR', 'DWI_dir96_RL', 'DWI_dir97_LR', 'DWI_dir97_RL']:
                        for dtype in ['nii.gz', 'bvec', 'bval']:
                            shutil.copy(f'{diff}/{sub}/unprocessed/3T/Diffusion/{sub}_3T_{direction}.{dtype}',f'/data2/HCP_Raw/raw_data/sub-{sub}/dwi')
                    shutil.rmtree(f'/data2/HCP_Raw/raw_data/sub-{sub}/dwi/{sub}')
        
bids_root = '/data2/HCP_Raw/raw_data'
for sub in os.listdir(bids_root):
    func = f'{bids_root}/{sub}/func'
    if 'sub' in sub:
        for img in os.listdir(func):
            if 'REST1_LR' in img:
               os.rename(f'{func}/{img}', f'{func}/{sub}_task-rest_acq-LR_run-1_bold.nii.gz')
            if 'REST1_RL' in img:
                os.rename(f'{func}/{img}', f'{func}/{sub}_task-rest_acq-RL_run-1_bold.nii.gz')
            if 'REST2_LR' in img:
                 os.rename(f'{func}/{img}', f'{func}/{sub}_task-rest_acq-LR_run-2)_bold.nii.gz')
            if 'REST2_RL' in img:
                 os.rename(f'{func}/{img}', f'{func}/{sub}_task-rest_acq-RL_run-2_bold.nii.gz')
            if 'WM_LR' in img:
                 os.rename(f'{func}/{img}', f'{func}/{sub}_task-WM_acq-LR_bold.nii.gz')
            if 'WM_RL' in img:
                 os.rename(f'{func}/{img}', f'{func}/{sub}_task-WM_acq-RL_bold.nii.gz')
                    
for sub in os.listdir(bids_root):
    anat = f'{bids_root}/{sub}/anat'
    if 'sub' in sub:
        for img in os.listdir(anat):
            for i in ['T1w','T2w']:
                if i in img:
                    os.rename(f'{anat}/{img}',f'{anat}/{sub}_{i}.nii.gz')
 
bids_root = '/data2/HCP_Raw/raw_data'                   
for sub in os.listdir(bids_root):
    dwi = f'{bids_root}/{sub}/dwi'
    if 'sub' in sub:
        for img in os.listdir(dwi):
            for ftype in ['nii.gz','bval','bvec']:
                if ftype in img:
                    for dire in ['95','96','97']:
                        if dire in img:
                            for acq in ['RL','LR']:
                                if acq in img:
                                    os.rename(f'{dwi}/{img}', f'{dwi}/{sub}_acq-{acq}_dir-{dire}_dwi.{ftype}')
                
                        
                
                



with ZipFile('/data2/HCP_Raw/source_data/105923/105923_3T_rfMRI_REST1_unproc.zip', 'r') as zipObj:
    zipObj.extractall('/data2/HCP_Raw/raw_data/sub-105923/func')
    
shutil.move('/data2/HCP_Raw/raw_data/sub-105923/func/105923/unprocessed/3T/tfMRI_WM_LR/105923_3T_tfMRI_WM_LR.nii.gz','/data2/HCP_Raw/raw_data/sub-105923/func')
shutil.rmtree('/data2/HCP_Raw/raw_data/sub-105923/func/105923')



