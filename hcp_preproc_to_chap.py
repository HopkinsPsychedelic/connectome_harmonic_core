#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 00:21:55 2020

@author: bwinston
"""
#only works for test-retest two sessions right now, only diffusion and structural (no func yet)
#TODO save out very inflated surfaces also

from zipfile import ZipFile
import os
import input_output as inout
import construct_harmonics as ch
import shutil

def hcp_chapper(args, sub, user_info):
    print(f'[CHAP] Creating directories for HCP subject {sub}')
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}') #intermediate sub folder
    for ses in ['ses-test', 'ses-retest']:
        user_info[f'{sub}_info'][ses], user_info[f'{sub}_info'][ses]['surfs'] = {}, {} #where hcp surface files will go
        inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}') #intermediate ses folder
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #chap output ses folder
        if args.fprep_dir: 
            inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc/{sub}/{ses}/func') #hcp func folder
        hcp_types = ['Structural', 'Diffusion', 'REST1']
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
        user_info[f'{sub}_info'][ses]['surfs']['lh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.L.white.32k_fs_LR.surf.gii' #hcp left hem
        user_info[f'{sub}_info'][ses]['surfs']['rh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.R.white.32k_fs_LR.surf.gii' #hcp right hem
        if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/10000000_endpoints.vtk'): #endpoints have been generated previously, skip mrtrix pipeline
            print('[CHAP] Endpoints already detected')
        else: #streamlines haven't been generated before, so run mrtrix diffusion pipeline with 10 million streamlines
            print('[CHAP] Running MRtrix commands to generate streamline endpoints...')
            os.system(f'bash /home/neuro/repo/run_mrtrix_diffusion_pipeline.sh {diffusion_dir}/data.nii.gz {diffusion_dir}/bvals {diffusion_dir}/bvecs  {struc_dir}/T1w_acpc_dc_restore_brain.nii.gz {diffusion_dir}/nodif_brain_mask.nii.gz {args.output_dir}/chap/sub-{sub}/{ses}/mrtrix 10000000')
            print('[CHAP] Removing intermediate files...')
            for file in ['10000000.tck', 'DWI.mif', '5TT.mif', 'WM_FODs.mif', '10000000_endpoints.tck']: #remove large intermediate files from chap mrtrix dir. won't delete endpoints.vtk, which is needed for harmonics. 
                os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{file}')
        user_info[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/10000000_endpoints.vtk' #streamline endpoints
        ch.construct_harmonics_calculate_spectra(args, sub, ses, user_info, multises = True) #run chcs function
        shutil.rmtree(f'{args.output_dir}/hcp_preproc/sub-{sub}/{ses}') #remove intermediate ses folder recursively



        
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