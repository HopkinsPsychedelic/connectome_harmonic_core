#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 16:01:27 2020
@author: bwinston
used in CHAP entrypoint_script

    sc - array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    si - array of surface indices of size (N_verts*2, 3) where each row defines indices of triange in surface mesh 
    ec - array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]
    tol(tolerance) - search radius of nearest neighbor search for matching endpoints to surface vertices
    NNnum - number of nearest neighboring surface vertices to assign to each endpoint
"""
import decomp as dcp
import input_output as inout
import utility_functions as util
import subprocess
import numpy as np
import os
import matrix_methods as mm
import compute_spectra as cs
from scipy import sparse

def qsi_chap(user_info, args, sub):
    user_info[f'{sub}_info']['streamlines'] = [] #where streamlines files will go
    print(f'[CHAP] Reconstructing surfaces for {sub}...')
    os.system(f'bash /home/neuro/repo/surface_resample.sh {args.surf_dir}/sub-{sub}/surf /home/neuro/repo') #run surface reconstruction script on freesurfer output
    if args.fprep_dir:
        user_info[f'{sub}_info']['func'] = [] #where functional files will go
    if any('ses' in x for x in os.listdir(f'{args.qsi_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.qsi_dir}/sub-{sub}'): 
            if 'ses' in ses:
                print(f'[CHAP] Locating files for {ses}...')
                user_info[f'{sub}_info'][ses] = {}
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'): #look in qsirecon output dir for tck
                    if 'tck' in file:
                        user_info[f'{sub}_info'][ses]['streamlines'] = file #streamlines list with each session's .tck
                        print('[CHAP] Located streamlines')
                if args.fprep_dir and os.path.exists(f'{args.fprep_dir}/sub-{sub}/{ses}/func'):
                    print(f'[CHAP] Detected functional images for {sub}')
                    for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/{ses}/func'): #look for preprocessed images in fmriprep dir 
                        if 'space-T1w_desc-preproc_bold.nii.gz' in file:
                            user_info[f'{sub}_info'][ses]['func'].append(file)  
            cs.prep_harmonics_bids(args, sub, user_info, multises, ses)                                  
    else: #if sub has just one session
        print('[CHAP] Detected only one session')
        multises = False
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                user_info[f'{sub}_info']['streamlines'] = file
                print('[CHAP] Located streamlines')
        if args.fprep_dir:
            for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/func'):
                if 'space-T1w_desc-preproc_bold.nii.gz' in file: 
                    user_info[f'{sub}_info']['func'].append(file)   
        cs.prep_harmonics_bids(args, sub, user_info, multises, ses = '')

def prep_harmonics_bids(args, sub, user_info, multises, ses):
    tck_name = user_info[f'{sub}_info'][ses]['streamlines'].split('/')[-1][:-4]
    print('[CHAP] Saving streamline endpoints and converting to vtk...')
    subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
    user_info[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk' #save endpoints to user_info
    for file in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}'):
        if '_endpoints.tck' in file:
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
    print('[CHAP] Finished MRtrix commands')
    #save output of surface reconstruction to user_info
    user_info[f'{sub}_info'][ses]['surfs']['lh'] = f'{args.surf_dir}/sub-{sub}/surf/lh.white.corresponded.vtk'
    user_info[f'{sub}_info'][ses]['surfs']['rh'] = f'{args.surf_dir}/sub-{sub}/surf/rh.white.corresponded.vtk'
    #run chcs function
    construct_harmonics_calculate_spectra(args, sub, ses, user_info, multises)

def construct_harmonics_calculate_spectra(args, sub, ses, user_info, multises): 
    if 'vtk' in user_info[f'{sub}_info'][ses]['surfs']['lh']: #probably used bids method, expects vtk surfaces
        sc,si=inout.read_vtk_surface_both_hem(user_info[f'{sub}_info'][ses]['surfs']['lh'], user_info[f'{sub}_info'][ses]['surfs']['rh']) #generate sc and si (see top of file) 
    else: #probably used hcp method, expects gii surfaces
        sc,si=inout.read_gifti_surface_both_hem(user_info[f'{sub}_info'][ses]['surfs']['lh'], user_info[f'{sub}_info'][ses]['surfs']['rh'], hcp = True)
        sc_inf,si_inf = inout.read_gifti_surface_both_hem(user_info[f'{sub}_info'][ses]['surfs']['lh_inf'], user_info[f'{sub}_info'][ses]['surfs']['rh_inf'], hcp = True)
    print('[CHAP] Saved surface coordinates and surface indices')
    ec=inout.read_streamline_endpoints(user_info[f'{sub}_info'][ses]['endpoints']) #read endpoint locations into numpy array (see top of file for definition of ec)
    print('[CHAP] Saved endpoint coordinates')
    print('[CHAP] Constructing surface matrix...')
    surf_mat=mm.construct_surface_matrix(sc,si) #construct surface matrix from sc and si    
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/surf_mat', surf_mat) #save out surface matrix
    print('[CHAP] Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc, ec, tol = args.tol, NNnum = args.nnum) #construct struc conn matrix from ec and sc (see matrix methods comments) 
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/struc_conn_mat', struc_conn_mat) #save out structural connectivity matrix
    print('[CHAP] Saved structural connectivity matrix')
    connectome = struc_conn_mat + surf_mat #sum connections and surface
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/connectome', connectome) #save out connectome 
    print('[CHAP] Saved connectome (surface + connections)')
    print('[CHAP] Computing harmonics...')
    vals,vecs=dcp.lapDecomp(connectome, args.evecs) #laplacian decomposition, returns eigenvals and eigenvectors (see decomp.py)
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis') #create visualization output directory
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vals',vals) #save np array eigenvals
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs',vecs) #save np array eigenvecs
    if args.mask_med_wall == True:
        masked_connectome = util.mask_connectivity_matrix(connectome, user_info[f'{sub}_info'][ses]['mask'])
        print('[CHAP] Masked out medial wall vertices; computing harmonics...')
        m_vals,m_vecs=dcp.lapDecomp(masked_connectome, args.evecs) #laplacian decomposition, returns eigenvals and eigenvectors (see decomp.py)
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/masked_vals',m_vals)
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/masked_vecs',m_vecs)   
    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,vecs) #harmonics.vtk
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_infl_harmonics.vtk',sc_inf,si_inf,vecs) #inflated harmonics.vtk
        print(f'[CHAP] Saved harmonics for {sub} {ses}')
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,vecs)
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_infl_harmonics.vtk',sc_inf,si_inf,vecs)
        print(f'[CHAP] Saved harmonics for {sub}')
    if args.fprep_dir: #if functional images are specified
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func') #func output folder
        for vol in user_info[f'{sub}_info']['func'][ses]: #each functional volume      
            full_path_lh = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/surfmapped_vol_lh.func.gii'
            full_path_rh = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/surfmapped_vol_rh.func.gii'
            task = inout.get_task(vol) #get taskname
            full_path_lh = full_path_lh[:-11] + f'task-{task}_' + full_path_lh[-11:]
            full_path_rh = full_path_rh[:-11] + f'task-{task}_' + full_path_rh[-11:]
            if 'acq' in vol:
                acq = inout.get_acq(vol)
                full_path_lh = full_path_lh[:-11] + f'acq-{acq}_' + full_path_lh[-11:]
                full_path_rh = full_path_rh[:-11] + f'acq-{acq}_' + full_path_rh[-11:]
            if 'run' in vol:
                run = inout.get_run(vol)
                full_path_lh = full_path_lh[:-11] + f'run-{run}_' + full_path_lh[-11:]
                full_path_rh = full_path_rh[:-11] + f'run-{run}_' + full_path_rh[-11:]
            print(f'[CHAP] Mapping {vol} to cortical surface') 
            os.system(f'bash /home/neuro/repo/volume_to_surface_map_fMRI.sh {args.surf_dir}/sub-{sub}/surf {args.fprep_dir}/sub-{sub}/{ses}/func/{vol} {full_path_lh} {full_path_rh}')
            bids_stuff = inout.get_bids_stuff(full_path_lh) #part of filename
            func_spectra(args, sub, ses, full_path_lh, full_path_rh, bids_stuff, vecs, vals)
    elif any('REST' in x for x in os.listdir(f'{args.hcp_dir}/{ses}')):
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func')
        func_dir = f'{args.output_dir}/chap/sub-{sub}/{ses}/func'
        if 'rest1' in user_info[f'{sub}_info'][ses]: #if rest1 data
            rest1_lr = user_info[f'{sub}_info'][ses]['rest1_lr']
            rest1_rl = user_info[f'{sub}_info'][ses]['rest1_rl']
            bids_stuff = f'sub-{sub}_{ses}_task-rest1'
            print('[CHAP] Extracting timecourse from HCP surface files...')
            os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {rest1_lr} COLUMN -metric CORTEX_LEFT {func_dir}/{bids_stuff}_acq-lr_hem-l.func.gii')
            os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {rest1_lr} COLUMN -metric CORTEX_RIGHT {func_dir}/{bids_stuff}_acq-lr_hem-r.func.gii')
            os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {rest1_rl} COLUMN -metric CORTEX_LEFT {func_dir}/{bids_stuff}_acq_-rl_hem-l.func.gii')
            os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {rest1_rl} COLUMN -metric CORTEX_RIGHT {func_dir}/{bids_stuff}_acq_-rl_hem-r.func.gii')
            timeseries_lr = cs.read_functional_timeseries(f'{func_dir}/{bids_stuff}_acq_-lr_hem-l.func.gii', f'{func_dir}/{bids_stuff}_acq-lr_hem-r.func.gii')
            timeseries_rl = cs.read_functional_timeseries(f'{func_dir}/{bids_stuff}_acq_-rl_hem-l.func.gii', f'{func_dir}/{bids_stuff}_acq-rl_hem-r.func.gii')
            timeseries_lr = (timeseries_lr - np.mean(timeseries_lr)) / np.std(timeseries_lr)
            timeseries_rl = (timeseries_rl - np.mean(timeseries_rl)) / np.std(timeseries_rl)
            print('[CHAP] Combining LR and RL PE direction scans...')
            timeseries = np.concatenate((timeseries_lr, timeseries_rl), axis=1)  
            func_spectra(args, sub, ses, timeseries, bids_stuff, vecs, vals)
    print(f'[CHAP] Finished session: {ses}')

def func_spectra(args, sub, ses, timeseries, bids_stuff, vecs, vals):
    #read functional timeseries 
    if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra'):
        print('[CHAP] Spectra already computed. If you want to run again, delete the old stuff, chap')
    else:
        #power spectra
        print(f'[CHAP] Computing mean, dynamic, and normalized power spectra for {bids_stuff} scan...')
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra')
        mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs) #average power over the whole scan (average of dynamic for each harmonic)
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_mean_power_spectrum', mean_power_spectrum)
        dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals) #power at each TR
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_dynamic_power_spectrum', dynamic_power_spectrum)
        normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_normalized_power_spectrum', normalized_power_spectrum)
        print(f'[CHAP] Saved power spectra for {bids_stuff} scan')
        #energy spectra
        print(f'[CHAP] Computing mean and dynamic energy spectra for {bids_stuff} scan...')
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra')
        mean_energy_spectrum = cs.mean_energy_spectrum(timeseries, vecs, vals) #average energy over the whole scan (average of dynamic for each harmonic)
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra/{bids_stuff}_mean_energy_spectrum', mean_energy_spectrum)
        dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals) #energy at each TR
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra/{bids_stuff}_dynamic_energy_spectrum', dynamic_energy_spectrum)
        print(f'[CHAP] Saved mean and dynamic energy spectra for {bids_stuff} scan')
        #reconstruction spectrum
        print(f'[CHAP] Computing dynamic reconstruction spectrum for {bids_stuff} scan...')
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/reconspectra')
        dynamic_reconstruction_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals) #takes on negative values
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/reconspectra/{bids_stuff}_dynamic_reconstruction_spectrum', dynamic_reconstruction_spectrum)
        print(f'[CHAP] Saved dynamic reconstruction spectrum for {bids_stuff} scan')

                
'''
wb_command -cifti-separate in.dtseries.nii COLUMN -metric CORTEX_LEFT out.func.gii 
''' 

  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
                
                
                