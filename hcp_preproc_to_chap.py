#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 00:21:55 2020

@author: bwinston
"""
#need to fix qsi method in terms of user_info and session stuff


from zipfile import ZipFile
import os
import input_output as inout
import construct_harmonics as ch

def hcp_chapper(args, sub, user_info):
    print(f'[CHAP] Creating directories for HCP subject {sub}')
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}')
    for ses in ['test', 'retest']:
        user_info[f'{sub}_info'][ses], user_info[f'{sub}_info'][ses]['surfs'] = {}, {} 
        inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}')
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/ses-{ses}')
        if args.fprep_dir:
            inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/{sub}/ses-{ses}/func')         
        for zipdir in os.listdir(f'{args.hcp_dir}/ses-{ses}'):
            if sub in zipdir and 'md5' not in zipdir:
                for bids_type in ['Structural', 'Diffusion']:
                    if bids_type in zipdir:
                        with ZipFile(f'{args.hcp_dir}/ses-{ses}/{zipdir}', 'r') as zipObj:
                            print(f'[CHAP] Unzipping {sub} {ses} {bids_type} directory')
                            zipObj.extractall(f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}/{bids_type}')       
        diffusion_dir = f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}/Diffusion/{sub}/T1w/Diffusion' 
        struc_dir = f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}/Structural/{sub}/T1w'
        user_info[f'{sub}_info'][ses]['surfs']['lh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.lh.white.32k_fs_LR.surf.gii'
        user_info[f'{sub}_info'][ses]['surfs']['rh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.rh.white.32k_fs_LR.surf.gii'
        if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix/10000000_endpoints.vtk'):
            print(f'[CHAP] Endpoints already detected')
            user_info[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix/10000000_endpoints.vtk'
        else: 
            print(f'[CHAP] Running MRtrix commands')
            os.system(f'bash /home/neuro/repo/run_mrtrix_diffusion_pipeline.sh {diffusion_dir}/data.nii.gz {diffusion_dir}/bvals {diffusion_dir}/bvecs  {struc_dir}/T1w_acpc_dc_restore_brain.nii.gz {diffusion_dir}/nodif_brain_mask.nii.gz {args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix 10000000')
            os.remove(f'{args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix/10000000.tck')
            os.remove(f'{args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix/DWI.mif')
            user_info[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/ses-{ses}/mrtrix/10000000_endpoints.vtk'
        ch.construct_harmonics_calculate_spectra(args, sub, ses, user_info, multises = True)   



        
'''               
   getattr(args, f'hcp_{ses}_dir')     
#run mrtrix pipeline on user_info[f'{sub}_info'][ses][dwi]
#do surface stuff for surfs
        
        
        if any('ses' in x for x in os.listdir(f'{args.hcp_dir}/sub-{sub}')):
            multises = True
            print(f'[CHAP] Detected multiple sessions for {sub}')
            for ses in os.listdir(f'{args.hcp_dir}/sub-{sub}'): 
                print(f'[CHAP] Locating files for {ses}...')
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}')
                for file in os.listdir(f'{args.hcp_dir}/sub-{sub}/{ses}/{sub}'):
                    
                    
                    #user_info[f'{sub}_info'][ses]['dwi'] = f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}/Diffusion/{sub}/T1w/Diffusion/data.nii.gz'
        #print(f'[CHAP] Found dwi data')
        #user_info[f'{sub}_info'][ses]['surfs'] = [f'{args.output_dir}/hcp_preproc/sub-{sub}/ses-{ses}/Structural/{sub}/T1w/fsaverage_LR32k/{sub}.{hem}.white.32k_fs_LR.surf.gii' for hem in ['L', 'R']]
        #print(f'[CHAP] Found surface data')
                
'''