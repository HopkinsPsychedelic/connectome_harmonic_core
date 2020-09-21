#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 17:19:49 2020

@author: bwinsto2

very lazy and inneficient code to bidsify HCP unprocessed data - pls don't judge

"""

from zipfile import ZipFile
import os, shutil
source_root = '/data2/HCP_Raw/source_data'
bids_root = '/data2/HCP_Raw/raw_data'

#create BIDS folders
for sub in os.listdir(source_root):
    os.mkdir(f'{bids_root}/sub-{sub}')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-test')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-retest')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-test/anat')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-test/func')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-test/dwi')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-retest/anat')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-retest/func')
    os.mkdir(f'{bids_root}/sub-{sub}/ses-retest/dwi')
    
#unzip and put into sub/ses/func
for sub in os.listdir(source_root):
    for subj_dir in os.listdir(f'{source_root}/{sub}'):
        for ses in ['Test','Retest']:
            func = f'{bids_root}/sub-{sub}/ses-{ses.lower()}/func'
            for zipdir in os.listdir(f'{source_root}/{str(sub)}/Session{ses}'):
                if 'md5' not in zipdir:
                    for dtype in ['tfMRI_WM','rfMRI_REST2','rfMRI_REST1']:
                        if dtype in zipdir:              
                            with ZipFile(f'{source_root}/{sub}/Session{ses}/{sub}_3T_{dtype}_unproc.zip', 'r') as zipObj:
                                zipObj.extractall(func)
                                for direction in ['LR', 'RL']:
                                    shutil.copy(f'{func}/{sub}/unprocessed/3T/{dtype}_{direction}/{sub}_3T_{dtype}_{direction}.nii.gz',f'{bids_root}/sub-{sub}/ses-{ses.lower()}/func')
                            shutil.rmtree(f'{bids_root}/sub-{sub}/ses-{ses.lower()}/func/{sub}')
    
#unzip and put into anat
for sub in os.listdir(source_root):
    for subj_dir in os.listdir(f'{source_root}/{sub}'):
        for ses in ['Test','Retest']: 
            anat = f'{bids_root}/sub-{sub}/ses-{ses.lower()}/anat'
            for zipdir in os.listdir(f'{source_root}/{str(sub)}/Session{ses}'):
                if 'md5' not in zipdir:
                    if 'Struc' in zipdir:              
                        with ZipFile(f'{source_root}/{sub}/Session{ses}/{sub}_3T_Structural_unproc.zip', 'r') as zipObj:
                            zipObj.extractall(anat)
                            for direction in ['T1w_MPR1', 'T2w_SPC1']:
                                shutil.copy(f'{anat}/{sub}/unprocessed/3T/{direction}/{sub}_3T_{direction}.nii.gz',f'{bids_root}/sub-{sub}/ses-{ses.lower()}/anat')
                        shutil.rmtree(f'{bids_root}/sub-{sub}/ses-{ses.lower()}/anat/{sub}')

#dwi
for sub in os.listdir(source_root):
    for subj_dir in os.listdir(f'{source_root}/{sub}'):
        for ses in ['Test','Retest']:
            diff = f'{bids_root}/sub-{sub}/ses-{ses.lower()}/dwi'
            for zipdir in os.listdir(f'{source_root}/{str(sub)}/Session{ses}'):
                if 'md5' not in zipdir:
                    if 'Diff' in zipdir:              
                        with ZipFile(f'{source_root}/{sub}/Session{ses}/{sub}_3T_Diffusion_unproc.zip', 'r') as zipObj:
                            zipObj.extractall(diff)
                            for direction in ['DWI_dir95_LR', 'DWI_dir95_RL', 'DWI_dir96_LR', 'DWI_dir96_RL', 'DWI_dir97_LR', 'DWI_dir97_RL']:
                                for dtype in ['nii.gz', 'bvec', 'bval']:
                                    shutil.copy(f'{diff}/{sub}/unprocessed/3T/Diffusion/{sub}_3T_{direction}.{dtype}',f'{bids_root}/sub-{sub}/ses-{ses.lower()}/dwi')
                            shutil.rmtree(f'{bids_root}/sub-{sub}/ses-{ses.lower()}/dwi/{sub}')

#rename func                
for sub in os.listdir(bids_root):
    if 'sub' in sub:
        for ses in os.listdir(f'{bids_root}/{sub}'):
            func = f'{bids_root}/{sub}/{ses}/func'
            for img in os.listdir(func):
                if 'REST1_LR' in img:
                   os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-rest_acq-LR_run-1_bold.nii.gz')
                if 'REST1_RL' in img:
                    os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-rest_acq-RL_run-1_bold.nii.gz')
                if 'REST2_LR' in img:
                     os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-rest_acq-LR_run-2_bold.nii.gz')
                if 'REST2_RL' in img:
                     os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-rest_acq-RL_run-2_bold.nii.gz')
                if 'WM_LR' in img:
                     os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-WM_acq-LR_bold.nii.gz')
                if 'WM_RL' in img:
                     os.rename(f'{func}/{img}', f'{func}/{sub}_{ses}_task-WM_acq-RL_bold.nii.gz')

#rename anat                    
for sub in os.listdir(bids_root):
    if 'sub' in sub:
        for ses in os.listdir(f'{bids_root}/{sub}'):
            anat = f'{bids_root}/{sub}/{ses}/anat'
            for img in os.listdir(anat):
                for i in ['T1w','T2w']:
                    if i in img:
                        os.rename(f'{anat}/{img}',f'{anat}/{sub}_{ses}_{i}.nii.gz')
     
#rename dwi              
for sub in os.listdir(bids_root):
    if 'sub' in sub:
        for ses in os.listdir(f'{bids_root}/{sub}'):
            dwi = f'{bids_root}/{sub}/{ses}/dwi'
            for img in os.listdir(dwi):
                for ftype in ['nii.gz','bval','bvec']:
                    if ftype in img:
                        for dire in ['95','96','97']:
                            if dire in img:
                                for acq in ['RL','LR']:
                                    if acq in img:
                                        os.rename(f'{dwi}/{img}', f'{dwi}/{sub}_{ses}_acq-{acq}_dir-{dire}_dwi.{ftype}')
                    
    
