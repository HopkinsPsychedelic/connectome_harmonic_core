#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 00:21:55 2020

@author: bwinston
"""

from zipfile import ZipFile
import os
import input_output as inout
import construct_harmonics as ch
import shutil
import numpy as np

def hcp_chapper(args, sub, u):
    print(f'[CHAP] Creating directories for HCP subject {sub}') 
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}') #intermediate sub folder
    if os.path.exists(f'{args.hcp_dir}/ses-test') == False: #if regular HCP (one session)
         ses = '' #pass ses as empty parameter to next fxn
         multises = False
         prep_for_cs(args, sub, u, multises, ses)
    else:
         for ses in ['ses-test', 'ses-retest']: #test-retest subjects
             multises = True
             prep_for_cs(args, sub, u, multises, ses)
            
def prep_for_cs(args, sub, u, multises, ses):   
    u[f'{sub}_info'][ses], u[f'{sub}_info'][ses]['surfs'] = {}, {} #where hcp surface files will go
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}') #intermediate ses folder
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #chap output ses folder
    if any('REST' in x for x in os.listdir(f'{args.hcp_dir}/{ses}')): #if rsfc data is there
        inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/func') #hcp func folder
    hcp_types = ['Structural', 'Diffusion', 'REST1', 'REST2', 'WM'] #etc.
    for hcp_type in hcp_types:
        if os.path.exists(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/{hcp_type}'): #data was unzipped before
            hcp_types.remove(hcp_type)       
    for zipdir in os.listdir(f'{args.hcp_dir}/{ses}'):
        if sub in zipdir and 'md5' not in zipdir:
            for bids_type in hcp_types:
                if bids_type in zipdir:
                    with ZipFile(f'{args.hcp_dir}/{ses}/{zipdir}', 'r') as zipObj:
                        print(f'[CHAP] Unzipping {sub} {ses} session {bids_type} directory')
                        zipObj.extractall(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/{bids_type}')      
    diffusion_dir = f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/Diffusion/{sub}/T1w/Diffusion' #diffusion path in intermediate dir
    struc_dir = f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/Structural/{sub}/T1w' #struc path in intermediate dir
    u[f'{sub}_info'][ses]['surfs']['lh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.L.white.32k_fs_LR.surf.gii' #hcp left hem
    u[f'{sub}_info'][ses]['surfs']['rh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.R.white.32k_fs_LR.surf.gii' #hcp right hem
    u[f'{sub}_info'][ses]['surfs']['lh_inf'] = f'{struc_dir}/fsaverage_LR32k/{sub}.L.very_inflated.32k_fs_LR.surf.gii' #hcp left hem inflated
    u[f'{sub}_info'][ses]['surfs']['rh_inf'] = f'{struc_dir}/fsaverage_LR32k/{sub}.R.very_inflated.32k_fs_LR.surf.gii' #hcp right hem inflated
    for n in ['1','2']:
        u[f'{sub}_info'][ses][f'rest{n}_lr'] = f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/REST{n}/{sub}/MNINonLinear/Results/rfMRI_REST{n}_LR/rfMRI_REST{n}_LR_Atlas_hp2000_clean.dtseries.nii'
        u[f'{sub}_info'][ses][f'rest{n}_rl'] = f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}/REST{n}/{sub}/MNINonLinear/Results/rfMRI_REST{n}_RL/rfMRI_REST{n}_RL_Atlas_hp2000_clean.dtseries.nii'
    if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/10000000_endpoints.vtk'): #endpoints have been generated previously, skip mrtrix pipeline
        print('[CHAP] Endpoints already detected')
    else: #streamlines haven't been generated before, so run mrtrix diffusion pipeline with 10 million streamlines
        print('[CHAP] Running MRtrix commands to generate streamline endpoints...')
        os.system(f'bash /home/neuro/repo/run_mrtrix_diffusion_pipeline.sh {diffusion_dir}/data.nii.gz {diffusion_dir}/bvals {diffusion_dir}/bvecs  {struc_dir}/T1w_acpc_dc_restore_brain.nii.gz {diffusion_dir}/nodif_brain_mask.nii.gz {args.output_dir}/chap/sub-{sub}/{ses}/mrtrix 10000000')
        print('[CHAP] Removing intermediate files...')
        for file in ['10000000.tck', 'DWI.mif', '5TT.mif', 'WM_FODs.mif', '10000000_endpoints.tck']: #remove large intermediate files from chap mrtrix dir. won't delete endpoints.vtk, which is needed for harmonics. 
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{file}')
    u[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/10000000_endpoints.vtk' #streamline endpoints
    ch.construct_harmonics_calculate_spectra(args, sub, ses, u, multises) #run chcs function
    shutil.rmtree(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}') #remove intermediate ses folder recursively
    
#/Users/bwinston/Documents/fMRI/BIDS/HCP_Preproc/









