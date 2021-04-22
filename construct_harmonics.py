
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
import utility_functions as uts
import subprocess
import numpy as np
import os
import matrix_methods as mm
import compute_spectra as cs
from scipy import sparse
from itertools import product

def qsi_chap(u, args, sub):
    u[f'{sub}_info']['streamlines'] = [] #where streamlines files will go
    print(f'[CHAP] Reconstructing surfaces for {sub}...')
    #for local testing, TODO: comment out/remove for production
    #Tring to remove volume to surface reconstruction for fsLR space 
    #os.system(f'bash /home/quintin_frerichs/connectome_harmonic_core/surface_resample.sh {args.surf_dir}sub-{sub}/surf /home/quintin_frerichs/connectome_harmonic_core') #run surface reconstruction script on freesurfer output
    #os.system(f'bash /home/neuro/repo/surface_resample.sh {args.surf_dir}/sub-{sub}/surf /home/neuro/repo') #run surface reconstruction script on freesurfer output
    #print(f'[CHAP] Reconstructing surfaces for {sub} complete...')
    if args.fprep_dir:
        u[f'{sub}_info']['func'] = [] #where functional files will go
    if any('ses' in x for x in os.listdir(f'{args.qsi_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.qsi_dir}/sub-{sub}'): 
            if 'ses' in ses:
                print(f'[CHAP] Locating files for {ses}...')
                u[f'{sub}_info'][ses] = {}
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'): #look in qsirecon output dir for tck
                    if 'tck' in file:
                        u[f'{sub}_info'][ses]['streamlines'] = file #streamlines list with each session's .tck
                        print('[CHAP] Located streamlines')
                if args.fprep_dir and os.path.exists(f'{args.fprep_dir}/sub-{sub}/{ses}/func'):
                    print(f'[CHAP] Detected functional images for {sub}')
                    for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/{ses}/func'): #look for preprocessed images in fmriprep dir 
                        if 'space-T1w_desc-preproc_bold.nii.gz' in file:
                            u[f'{sub}_info'][ses]['func'].append(file)  
            prep_harmonics_bids(args, sub, u, multises, ses)                                  
    else: #if sub has just one session
        print('[CHAP] Detected only one session')
        multises = False
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                u[f'{sub}_info']['streamlines'] = file
                print('[CHAP] Located streamlines')
        if args.fprep_dir:
            for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/func'):
                if 'space-T1w_desc-preproc_bold.nii.gz' in file: 
                    u[f'{sub}_info']['func'].append(file)   
        prep_harmonics_bids(args, sub, u, multises, ses = '')

def prep_harmonics_bids(args, sub, u, multises, ses):
    if(ses == ''):
        tck_name = u[f'{sub}_info']['streamlines'].split('/')[-1][:-4]
        print('[CHAP] Saving streamline endpoints and converting to vtk...')
        #for local testing, TODO: comment out 
        subprocess.check_call("/home/quintin_frerichs/connectome_harmonic_core/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}sub-{sub}/dwi', tck_name, f'{args.output_dir}chap/sub-{sub}'), shell=True) #run mrtrix bash script
        #for dcoker version
        #subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/'), shell=True) #run mrtrix bash script
        u[f'{sub}_info']['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{tck_name}_endpoints.vtk' #save endpoints to u
        for file in os.listdir(f'{args.output_dir}chap/sub-{sub}/'):
            if '_endpoints.tck' in file:
                os.remove(f'{args.output_dir}/chap/sub-{sub}/{file}') #remove endpoints tck
        print('[CHAP] Finished MRtrix commands')
        #save output of surface reconstruction to u
        u[f'{sub}_info']['surfs'] = {}
        u[f'{sub}_info']['surfs']['lh'] = f'{args.surf_dir}sub-{sub}/surf/lh.white.corresponded.vtk'
        u[f'{sub}_info']['surfs']['rh'] = f'{args.surf_dir}sub-{sub}/surf/rh.white.corresponded.vtk'
        #run chcs function
        construct_harmonics_calculate_spectra(args, sub, ses, u, multises)    
    else:
        tck_name = u[f'{sub}_info']['streamlines'].split('/')[-1][:-4]
        print('[CHAP] Saving streamline endpoints and converting to vtk...')
        #for local testing, TODO: comment out 
        subprocess.check_call("/home/quintin_frerichs/connectome_harmonic_core/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
        #for dcoker version
        #subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/'), shell=True) #run mrtrix bash script
        u[f'{sub}_info']['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk' #save endpoints to u
        for file in os.listdir(f'{args.output_dir}chap/sub-{sub}/{ses}'):
            if '_endpoints.tck' in file:
                os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
        print('[CHAP] Finished MRtrix commands' + f'{args.surf_dir}')
        #save output of surface reconstruction to u
        u[f'{sub}_info']['surfs'] = {}
        u[f'{sub}_info']['surfs']['lh'] = f'{args.surf_dir}sub-{sub}/surf/lh.white.corresponded.vtk'
        u[f'{sub}_info']['surfs']['rh'] = f'{args.surf_dir}sub-{sub}/surf/rh.white.corresponded.vtk'
        #run chcs function
        construct_harmonics_calculate_spectra(args, sub, ses, u, multises)

def construct_harmonics_calculate_spectra(args, sub, ses, u, multises): 
    #TODO - implementation for multises, removed [ses]
    is_hcp = False
    if 'vtk' in u[f'{sub}_info']['surfs']['lh']: #probably used bids method, expects vtk surfaces
        print("Surfaces:" + u[f'{sub}_info']['surfs']['lh'])
        sc,si=inout.read_vtk_surface_both_hem(u[f'{sub}_info']['surfs']['lh'], u[f'{sub}_info']['surfs']['rh']) #generate sc and si (see top of file)
    else: #probably used hcp method, expects gii surfaces
        sc,si=inout.read_gifti_surface_both_hem(u[f'{sub}_info']['surfs']['lh'], u[f'{sub}_info']['surfs']['rh'], hcp = True)
        sc_inf,si_inf = inout.read_gifti_surface_both_hem(u[f'{sub}_info']['surfs']['lh_inf'], u[f'{sub}_info']['surfs']['rh_inf'], hcp = True)
        is_hcp  = True
    print('[CHAP] Saved surface coordinates and surface indices')
    ec=inout.read_streamline_endpoints(u[f'{sub}_info']['endpoints']) #read endpoint locations into numpy array (see top of file for definition of ec)
    print('[CHAP] Saved endpoint coordinates')
    surf_mat=mm.construct_surface_matrix(sc,si) #construct surface matrix from sc and si    
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/surf_mat', surf_mat) #save out surface matrix
    print('[CHAP] Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc, ec, tol = args.tol, NNnum = args.nnum) #construct struc conn matrix from ec and sc (see matrix methods comments) 
    connectome = struc_conn_mat + surf_mat #sum connections and surface
    if args.hcp_dir:
        connectome = uts.mask_connectivity_matrix(connectome, u['mask']) #mask medial wall
        print('[CHAP] Masked out medial wall vertices; computing harmonics...')
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/connectome', connectome) #save out connectome 
    print('[CHAP] Saved connectome (surface + connections)')
    print('[CHAP] Computing harmonics...')
    vals,vecs=dcp.lapDecomp(connectome, args.evecs) #laplacian decomposition, returns eigenvals and eigenvectors (see decomp.py)
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/vis') #create visualization output directory
    np.save(f'{args.output_dir}/chap/sub-{sub}/vals',vals) #save np array eigenvals
    np.save(f'{args.output_dir}/chap/sub-{sub}/vecs',vecs) #save np array eigenvecs
    np.save(f'{args.output_dir}/chap/sub-{sub}/surface_indices',si) #save np array eigenvecs
    np.save(f'{args.output_dir}/chap/sub-{sub}/surface_coordinates',sc) #save np array eigenvecs


    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,vecs) #harmonics.vtk
        inout.visualize_harmonics(sub, ses, si.astype(int), sc.astype(int), vecs, f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/')

        if is_hcp:
            inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_infl_harmonics.vtk',sc_inf,si_inf,vecs) #inflated harmonics.vtk
        print(f'[CHAP] Saved harmonics for {sub} {ses}')
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,vecs)
        inout.visualize_harmonics(sub, 'single', si.astype(int), sc.astype(int), vecs, f'{args.output_dir}/chap/sub-{sub}/vis/')
        if is_hcp:
            inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_infl_harmonics.vtk',sc_inf,si_inf,vecs)
        print(f'[CHAP] Saved harmonics for {sub}')
    if args.fprep_dir: #if functional images are specified (BIDS method)
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func') #func output folder
        for vol in u[f'{sub}_info']['func'][ses]: #each functional volume      
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
    if is_hcp:        
        if any('REST' in x for x in os.listdir(f'{args.hcp_dir}/{ses}')): #functional stuff, HCP method
            func_dir = f'{args.output_dir}/chap/sub-{sub}/{ses}/func'    
            inout.if_not_exist_make(func_dir)
            if 'rest1_lr' in u[f'{sub}_info'][ses]: #if rest1 data, assuming also rest2
                for n in ['1','2']:
                    for dire in ['lr', 'rl']:
                        scan = u[f'{sub}_info'][ses][f'rest{n}_{dire}']
                        bids_stuff = f'sub-{sub}_{ses}_task-rest{n}_acq-{dire}'
                        print('[CHAP] Extracting timecourse from HCP surface files...')
                        os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {scan} COLUMN -metric CORTEX_LEFT {func_dir}/{bids_stuff}_hem-l.func.gii')
                        os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {scan} COLUMN -metric CORTEX_RIGHT {func_dir}/{bids_stuff}_hem-r.func.gii')
                        u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'] = cs.read_functional_timeseries(f'{func_dir}/{bids_stuff}_hem-l.func.gii', f'{func_dir}/{bids_stuff}_hem-r.func.gii')
                        u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'timeseries_rest{n}_{dire}'], u['mask'])
                    print(f'[CHAP] Combining LR and RL PE direction scans for REST{n}...')
                    u[f'{sub}_info'][ses][f'rest{n}_comb'] = inout.combine_pe(u[f'{sub}_info'][ses][f'timeseries_rest{n}_lr'], u[f'{sub}_info'][ses][f'timeseries_rest{n}_rl'])  
                    func_spectra(args, sub, ses, u[f'{sub}_info'][ses][f'rest{n}_comb'], f'REST{n}', bids_stuff, vecs, vals)
                for n, dire, hem in product(('1','2'), ('lr','rl'), ('l','r')): #remove giftis
                    os.remove(f'{func_dir}/sub-{sub}_{ses}_task-rest{n}_acq-{dire}_hem-{hem}.func.gii')
    print(f'[CHAP] Finished session: {ses}')

def func_spectra(args, sub, ses, timeseries, task, bids_stuff, vecs, vals):
    #read functional timeseries
    task_dir = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{task}'
    if os.path.exists(f'{task_dir}'):
        print('[CHAP] Spectra computed previously. If you want to run again, delete the old stuff, chap')
    else:
        inout.if_not_exist_make(f'{task_dir}')
        bids_stuff = bids_stuff[:-7]
        for spec in ['powerspectra', 'energyspectra','reconspectra']:
            inout.if_not_exist_make(f'{task_dir}/{spec}') 
        #power spectra
        print(f'[CHAP] Computing mean, dynamic, and normalized power spectra for {bids_stuff} scan...')
        mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs) #average power over the whole scan (average of dynamic for each harmonic)
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_mean_power_spectrum', mean_power_spectrum)
        dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals) #power at each TR
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_dynamic_power_spectrum', dynamic_power_spectrum)
        normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_normalized_power_spectrum', normalized_power_spectrum)
        print(f'[CHAP] Saved power spectra for {bids_stuff} scan')
        #energy spectra
        mean_energy_spectrum = cs.mean_energy_spectrum(timeseries, vecs, vals) #average energy over the whole scan (average of dynamic for each harmonic)
        np.save(f'{task_dir}/energyspectra/{bids_stuff}_mean_energy_spectrum', mean_energy_spectrum)
        dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals) #energy at each TR
        np.save(f'{task_dir}/energyspectra/{bids_stuff}_dynamic_energy_spectrum', dynamic_energy_spectrum)
        print(f'[CHAP] Saved energy spectra for {bids_stuff} scan')
        #reconstruction spectrum
        dynamic_reconstruction_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals) #takes on negative values
        np.save(f'{task_dir}/reconspectra/{bids_stuff}_dynamic_reconstruction_spectrum', dynamic_reconstruction_spectrum)
        print(f'[CHAP] Saved dynamic reconstruction spectrum for {bids_stuff} scan')

                
'''
wb_command -cifti-separate in.dtseries.nii COLUMN -metric CORTEX_LEFT out.func.gii 
''' 

    
  
    
  
    
  
                
                
                
