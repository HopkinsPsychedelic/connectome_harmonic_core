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

def construct_harmonics_calculate_spectra(args, sub, ses, u, multises): 
    sc,si=inout.read_gifti_surface_both_hem(u[f'{sub}_info'][ses]['surfs']['lh'], u[f'{sub}_info'][ses]['surfs']['rh'], hcp = True)
    print('[CHAP] Saved surface coordinates and surface indices')
    ec=inout.read_streamline_endpoints(u[f'{sub}_info'][ses]['endpoints']) #read endpoint locations into numpy array (see top of file for definition of ec)
    print('[CHAP] Saved endpoint coordinates')
    surf_mat=mm.construct_surface_matrix(sc,si) #construct surface matrix from sc and si    
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/surf_mat', surf_mat) #save out surface matrix
    print('[CHAP] Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc, ec, tol = args.tol, NNnum = args.nnum) #construct struc conn matrix from ec and sc (see matrix methods comments) 
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/struc_conn_mat', struc_conn_mat)
    connectome = struc_conn_mat + surf_mat #sum connections and surface
    if args.mask_med_wall==True:
        connectome = uts.mask_connectivity_matrix(connectome, u['mask']) #mask medial wall
        print('[CHAP] Masked out medial wall vertices; computing harmonics...')
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/connectome', connectome) #save out connectome 
    print('[CHAP] Saved connectome (surface + connections)')
    print('[CHAP] Computing harmonics...')
    vals,vecs=dcp.lapDecomp(connectome, args.evecs) #laplacian decomposition, returns eigenvals and eigenvectors (see decomp.py)
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis') #create visualization output directory
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vals',vals) #save np array eigenvals
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs',vecs) #save np array eigenvecs
    if args.mask_med_wall==True:
        unmasked_vecs = np.empty([64984,args.evecs])
        for ev in range(args.evecs):
            unmasked_vecs[:,ev]=uts.unmask_medial_wall(vecs[:,ev],u['mask'])
    else:
        unmasked_vecs = vecs
    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,unmasked_vecs) #harmonics.vtk
        print(f'[CHAP] Saved harmonics for {sub} {ses}')
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,unmasked_vecs)
        print(f'[CHAP] Saved harmonics for {sub}')
    if u[f'{sub}_info'][ses]['is_func']: #if functional images are specified (and BIDS method)
        inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func') #func output folder
        for dts in u[f'{sub}_info'][ses]['func']: #each functional volume
            bids_stuff = f'sub-{sub}_{inout.get_bids_stuff(dts)}' #e.g. sub-{sub}_ses-{ses}_task-{task}
            inout.dts_to_func_gii(dts, f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}') #extract cortex timeseries with connectome workbench
            u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = cs.read_functional_timeseries(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-l.func.gii', f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-r.func.gii') #func.gii to timeseries
            u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], u['mask'])
            func_spectra(args, sub, ses, u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], inout.get_task(dts), bids_stuff, vecs, vals)
    elif any('REST' in x for x in os.listdir(f'{args.hcp_dir}/{ses}')): #functional stuff, HCP method
        if args.skip_func == False:
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
    if os.path.exists(f'{task_dir}/powerspectra/{bids_stuff}_mean_power_spectrum.npy'):
        print(f'[CHAP] Spectra computed previously for {bids_stuff}. If you want to run again, delete the old stuff, chap')
    else:
        inout.if_not_exist_make(f'{task_dir}')
        for spec in ['powerspectra', 'energyspectra','reconspectra']:
            inout.if_not_exist_make(f'{task_dir}/{spec}') 
        #power spectra
        print(f'[CHAP] Computing mean, dynamic, and normalized power spectra for {bids_stuff}...')
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

                
    
  
    
  
    
  
                
                
                