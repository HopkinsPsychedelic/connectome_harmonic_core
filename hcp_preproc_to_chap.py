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
import compute_spectra as cs
import shutil
from itertools import product
import utility_functions as uts
import numpy as np

def hcp_chapper(args, sub, u):
    inout.if_not_exist_make(f'{args.output_dir}/chap_work/sub-{sub}') #intermediate sub folder
    if os.path.exists(f'{args.hcp_dir}/ses-test') == False: #if regular HCP (one session)
         ses = '' #pass ses as empty parameter to next fxn
         multises = False
         hcp_prep_for_ch(args, sub, u, multises, ses)
    else:
         for ses in ['ses-test', 'ses-retest']: #test-retest subjects
             multises = True
             hcp_prep_for_ch(args, sub, u, multises, ses)
            
def hcp_prep_for_ch(args, sub, u, multises, ses):  
    #create directories/dicts
    u[f'{sub}_info'][ses], u[f'{sub}_info'][ses]['surfs'] = {}, {} #where hcp surface files will go
    inout.if_not_exist_make(f'{args.output_dir}/chap_work/sub-{sub}/{ses}') #intermediate ses folder
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #chap output ses folder
    #check if there are HCP functional data
    u[f'{sub}_info'][ses]['hcp_types'] = ['REST1', 'REST2', 'WM','MOTOR','LANGUAGE','EMOTION','GAMBLING','SOCIAL','RELATIONAL']   
    for hcp_type in u[f'{sub}_info'][ses]['hcp_types']:
        if any(hcp_type in x for x in os.listdir(f'{args.hcp_dir}/{ses}')): #if func data are downloaded
            inout.if_not_exist_make(f'{args.output_dir}/chap_work/sub-{sub}/{ses}/func') #hcp func folder
            u[f'{sub}_info'][ses]['is_func'] = 'hcp'
            break
    if 'is_func' not in u[f'{sub}_info'][ses]: #if no functional 
        u[f'{sub}_info'][ses]['hcp_types'].clear()
    #rename Structural_preproc_extended dirs Freesurfer
    for download_dir in os.listdir(f'{args.hcp_dir}/{ses}'):
        if sub in download_dir and 'Structural_preproc_extended' in download_dir and 'md5' not in download_dir:
            os.rename(f'{args.hcp_dir}/{ses}/{download_dir}',f'{args.hcp_dir}/{ses}/{sub}_3T_Freesurfer.zip')
    u[f'{sub}_info'][ses]['hcp_types'].extend(['Structural','Diffusion','Freesurfer']) #add structural and diffusion (required) to hcp_types
    add_back = [] 
    hcp_types = u[f'{sub}_info'][ses]['hcp_types'].copy()
    for hcp_type in u[f'{sub}_info'][ses]['hcp_types']: #check if there are prev. data computed
        if os.path.exists(f'{args.output_dir}/chap_work/sub-{sub}/{ses}/{hcp_type}'): #data were unzipped before
            hcp_types.remove(hcp_type) 
            add_back.append(hcp_type)
    u[f'{sub}_info'][ses]['hcp_types'] = hcp_types #already unzipped ones are removed from hcp_types
    #unzip remaining HCP data
    for zipdir in os.listdir(f'{args.hcp_dir}/{ses}'):
        if sub in zipdir and 'md5' not in zipdir:
            for hcp_type in u[f'{sub}_info'][ses]['hcp_types']:
                if hcp_type in zipdir:
                    with ZipFile(f'{args.hcp_dir}/{ses}/{zipdir}', 'r') as zipObj:
                        print(f'[CHAP] Unzipping {sub} {ses} session {hcp_type} directory')
                        zipObj.extractall(f'{args.output_dir}/chap_work/sub-{sub}/{ses}/{hcp_type}') #extract to intermediate
    hcp_types = u[f'{sub}_info'][ses]['hcp_types'].copy()
    for hcp_type in u[f'{sub}_info'][ses]['hcp_types']: 
        if not os.path.exists(f'{args.output_dir}/chap_work/sub-{sub}/{ses}/{hcp_type}'): #if they don't have a task downloaded
            hcp_types.remove(hcp_type)
    u[f'{sub}_info'][ses]['hcp_types'] = hcp_types
    u[f'{sub}_info'][ses]['hcp_types'].extend(add_back) #add stuff back to hcp_types that were previously unzipped
    #define paths
    diffusion_dir = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/Diffusion/{sub}/T1w/Diffusion' #diffusion path in intermediate dir
    struc_dir = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/Structural/{sub}/T1w' #struc path in intermediate dir
    freesurfer_dir = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/Freesurfer/{sub}/T1w/{sub}' #freesurfer path in intermediate dir
    #save hcp surfaces in dict
    u[f'{sub}_info'][ses]['surfs']['lh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.L.white.32k_fs_LR.surf.gii' #hcp left hem
    u[f'{sub}_info'][ses]['surfs']['rh'] = f'{struc_dir}/fsaverage_LR32k/{sub}.R.white.32k_fs_LR.surf.gii' #hcp right hem
    #check for streamlines; if not, run mrtrix reconstruction script
    mrtrix_recon(u,sub,ses,args,f'{diffusion_dir}/data.nii.gz',f'{diffusion_dir}/bvals',f'{diffusion_dir}/bvecs',freesurfer_dir,f'{diffusion_dir}/nodif_brain_mask.nii.gz')
    #send to chcs fxn
    if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs.npy'):
        print('[CHAP] Harmonics already detected. Checking for spectra...')
        ch.check_func(args,sub,ses,u,np.load(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs.npy'),np.load(f'{args.output_dir}/chap/sub-{sub}/{ses}/vals.npy'))
    else:
        ch.construct_harmonics(args, sub, ses, u, multises) #run ch function
    shutil.rmtree(f'{args.output_dir}/chap_work/sub-{sub}/{ses}') #remove intermediate ses folder recursively

def hcp_spectra_prep(args,sub,ses,u,vecs,vals):  
    #func prep stuff
    func_dir = f'{args.output_dir}/chap/sub-{sub}/{ses}/func'  
    u[f'{sub}_info'][ses]['hcp_types'] = [i for i in u[f'{sub}_info'][ses]['hcp_types'] if i not in ('Structural', 'Diffusion', 'Freesurfer')] #get hcp types that aren't structural or diffusion or freesrufr
    if 'REST1' in u[f'{sub}_info'][ses]['hcp_types']:
        u[f'{sub}_info'][ses]['hcp_types'].remove('REST2') #don't need to run below twice
    for hcp_type in u[f'{sub}_info'][ses]['hcp_types']: #for each functional scan
        #rest stuff
        if hcp_type == 'REST1': 
            for n in ['1','2']:  
                u[f'{sub}_info'][ses][f'rest{n}_lr'] = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/REST{n}/{sub}/MNINonLinear/Results/rfMRI_REST{n}_LR/rfMRI_REST{n}_LR_Atlas_hp2000_clean.dtseries.nii'
                u[f'{sub}_info'][ses][f'rest{n}_rl'] = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/REST{n}/{sub}/MNINonLinear/Results/rfMRI_REST{n}_RL/rfMRI_REST{n}_RL_Atlas_hp2000_clean.dtseries.nii'      
                for dire in ['lr', 'rl']:
                    scan = u[f'{sub}_info'][ses][f'rest{n}_{dire}']
                    bids_stuff = f'sub-{sub}_{ses}_task-rest{n}_acq-{dire}'
                    print(f'[CHAP] Extracting timeseries from REST{n} {dire} direction dtseries...')
                    os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {scan} COLUMN -metric CORTEX_LEFT {func_dir}/{bids_stuff}_hem-l.func.gii')
                    os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {scan} COLUMN -metric CORTEX_RIGHT {func_dir}/{bids_stuff}_hem-r.func.gii')
                    u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'] = cs.read_functional_timeseries(f'{func_dir}/{bids_stuff}_hem-l.func.gii', f'{func_dir}/{bids_stuff}_hem-r.func.gii')
                    u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'], u['mask'])
                print(f'[CHAP] Concatenating LR and RL PE direction scans for REST{n}...')
                u[f'{sub}_info'][ses][f'rest{n}_comb'] = inout.combine_pe(u[f'{sub}_info'][ses][f'timeseries_rest{n}_lr'], u[f'{sub}_info'][ses][f'timeseries_rest{n}_rl'])  
                ch.func_spectra(args, sub, ses, u[f'{sub}_info'][ses][f'rest{n}_comb'], f'REST{n}', bids_stuff, vecs, vals)
            for n, dire, hem in product(('1','2'), ('lr','rl'), ('l','r')): #remove giftis
                os.remove(f'{func_dir}/sub-{sub}_{ses}_task-rest{n}_acq-{dire}_hem-{hem}.func.gii')
        #tasks
        else: #e.g. MOTOR, LANGUAGE
            u[f'{sub}_info'][ses][hcp_type] = {}
            results_dir = f'{args.output_dir}/chap_work/sub-{sub}/{ses}/{hcp_type}/{sub}/MNINonLinear/Results'
            inout.if_not_exist_make(f'{func_dir}/{hcp_type}')
            inout.if_not_exist_make(f'{func_dir}/{hcp_type}/movement_regressors')
            for dire in ['LR','RL']:
                #save behavioral files
                if not os.path.exists(f'{func_dir}/{hcp_type}/{dire}_EVs'):
                    shutil.copytree(f'{results_dir}/tfMRI_{hcp_type}_{dire}/EVs',f'{func_dir}/{hcp_type}/{dire}_EVs')
                for reg_file in os.listdir(results_dir):
                    if 'Movement' in reg_file:
                        shutil.copyfile(f'{results_dir}/{reg_file}',f'{func_dir}/{hcp_type}/movement_regressors/{dire}_reg_file')
                #process timeseries
                u[f'{sub}_info'][ses][hcp_type][dire] = f'{results_dir}/tfMRI_{hcp_type}_{dire}/tfMRI_{hcp_type}_{dire}_Atlas_MSMAll.dtseries.nii'
                bids_stuff = f'sub-{sub}_{ses}_task-{hcp_type}_acq-{dire}'
                inout.dts_to_func_gii(u[f'{sub}_info'][ses][hcp_type][dire], f'{func_dir}/{bids_stuff}')
                u[f'{sub}_info'][ses][hcp_type][dire] = cs.read_functional_timeseries(f'{func_dir}/{bids_stuff}_hem-l.func.gii', f'{func_dir}/{bids_stuff}_hem-r.func.gii')
                u[f'{sub}_info'][ses][hcp_type][dire] = uts.mask_timeseries(u[f'{sub}_info'][ses][hcp_type][dire],u['mask'])
                #np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{hcp_type}_{dire}.npy',u[f'{sub}_info'][ses][hcp_type][dire])
                os.remove(f'{func_dir}/{bids_stuff}_hem-l.func.gii')
                os.remove(f'{func_dir}/{bids_stuff}_hem-r.func.gii')
            #concatenate PE directions
            print(f'[CHAP] Concatenating LR and RL PE direction scans for {sub} {ses} {hcp_type} scan...')
            u[f'{sub}_info'][ses][hcp_type]['ts'] = inout.combine_pe(u[f'{sub}_info'][ses][hcp_type]['LR'],u[f'{sub}_info'][ses][hcp_type]['RL'])  
            ch.func_spectra(args,sub,ses,u[f'{sub}_info'][ses][hcp_type]['ts'],hcp_type,bids_stuff,vecs,vals)
                 
def mrtrix_recon(u,sub,ses,args,diff_preproc,bvals,bvecs,freesurfer_dir,diff_mask):                
    #check if endpoints already computed, if not run diffusion pipeline
    if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{args.streamlines}_endpoints.vtk'): #endpoints have been generated previously, skip mrtrix pipeline
        print('[CHAP] Endpoints already detected')
    else: #streamlines haven't been generated before, so run mrtrix diffusion pipeline with 10 million streamlines
        print('[CHAP] Running mrtrix commands to generate streamline endpoints. Could take up to a couple hours...')
        if args.diff_pipeline == 'msmt':
            os.system(f'bash /home/neuro/repo/msmt_5tt_mrtrix_diffusion_pipeline.sh {diff_preproc} {bvals} {bvecs} {freesurfer_dir} {diff_mask} {args.output_dir}/chap/sub-{sub}/{ses}/mrtrix {args.streamlines}')
        else:
            os.system(f'bash /home/neuro/repo/dhollander_mrtrix_pipeline.sh {diff_preproc} {bvals} {bvecs} {freesurfer_dir} {diff_mask} {args.output_dir}/chap/sub-{sub}/{ses}/mrtrix {args.streamlines}')
        print('[CHAP] Removing intermediate files...')
        for file in ['DWI.mif', '5TT.mif', 'WM_FODs.mif', f'{args.streamlines}_endpoints.tck', f'{args.streamlines}.tck']: #remove large intermediate files from chap mrtrix dir. won't delete endpoints.vtk, which is needed for harmonics. 
        #for file in ['DWI.mif', 'WM_FODs.mif', f'{args.streamlines}_endpoints.tck']: #remove large intermediate files from chap mrtrix dir. won't delete endpoints.vtk, which is needed for harmonics. 
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{file}')
        for item in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/'):
            if 'dwi2response' in item:
                shutil.rmtree(f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{item}')               
    u[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/mrtrix/{args.streamlines}_endpoints.vtk' #save streamline endpoints in dict
        
                
                
                
                
                
                
                
                
                
                
                
                
                
                




