#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 16:01:27 2020
@author: bwinston

    sc - array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    si - array of surface indices of size (N_verts*2, 3) where each row defines indices of triange in surface mesh 
    ec - array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]
    tol(tolerance) - search radius of nearest neighbor search for matching endpoints to surface vertices
    NNnum - number of nearest neighboring surface vertices to assign to each endpoint
"""
import decomp as dcp
import input_output as inout
import utility_functions as uts
import numpy as np
import os
import matrix_methods as mm
import compute_spectra as cs
from scipy import sparse
import bids_to_ch as bids
import hcp_preproc_to_chap as hcp_prep

def construct_harmonics(args, sub, ses, u, multises): 
    sc,si = inout.read_gifti_surface_both_hem(u[f'{sub}_info'][ses]['surfs']['lh'], u[f'{sub}_info'][ses]['surfs']['rh'], hcp = True)
    ec = inout.read_streamline_endpoints(u[f'{sub}_info'][ses]['endpoints']) #read endpoint locations into numpy array (see top of file for definition of ec)
    surf_mat = uts.mask_connectivity_matrix(mm.construct_surface_matrix(sc,si),u['mask']) #construct surface matrix from sc and si and mask
    #construct struc conn matrix from ec and sc (see matrix methods comments)
    struc_conn_mat = mm.construct_smoothed_connectivity_matrix(sc,si,ec,u['mask'], args.tol, args.sigma, args.epsilon, binarize=args.binarize)   
    zeromask = uts.get_zero_mask_from_connectivity(struc_conn_mat)
    struc_conn_mat = uts.mask_connectivity_matrix(struc_conn_mat,zeromask)
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/struc_conn_mat', struc_conn_mat)
    surf_mat = uts.mask_connectivity_matrix(surf_mat,zeromask)
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/surf_mat', surf_mat) #save out surface matrix
    #sum connections and surface
    connectome = struc_conn_mat + surf_mat 
    #save out connectome 
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/connectome', connectome) 
    print('[CHAP] Saved connectome (surface + long-range connections)')
    #compute harmonics
    vals,vecs = dcp.lapDecomp(connectome, args.evecs) #laplacian decomposition, returns eigenvals and eigenvecs (see decomp.py)
    vecs = uts.unmask_medial_wall_vecs(vecs,zeromask)
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vals',vals) #save np array eigenvals
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs',vecs) #save np array eigenvecs
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis') #create visualization output directory
    if args.mask_med_wall==True: #save unmasked vecs for visualization purposes
        unmasked_vecs=uts.unmask_medial_wall_vecs(vecs,u['mask'])
    else:
        unmasked_vecs = vecs
    #visualization
    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,unmasked_vecs) #save out harmonics.vtk
        print(f'[CHAP] Saved harmonics for {sub} {ses}')
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,unmasked_vecs) #save out harmonics.vtk
        print(f'[CHAP] Saved harmonics for {sub}')
    check_func(args,sub,ses,u,vecs,vals)

def check_func(args,sub,ses,u,vecs,vals):
    if args.skip_func == False:
        if 'is_func' in u[f'{sub}_info'][ses]: #func stuff
            inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}/func') #func output folder
            if u[f'{sub}_info'][ses]['is_func'] == 'cift': #bids method
                bids.bids_spectra_prep(args,sub,ses,u,vecs,vals)
            else:  #functional stuff, HCP method (it saves is_func elsewhere)
                hcp_prep.hcp_spectra_prep(args,sub,ses,u,vecs,vals)    
    print(f'[CHAP] Finished session: {ses}')

def func_spectra(args, sub, ses, timeseries, task, bids_stuff, vecs, vals): #for each timeseries
    task_dir = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{task}'
    #check for prev data
    if os.path.exists(f'{task_dir}/powerspectra/{bids_stuff}_mean_power_spectrum.npy'):
        print(f'[CHAP] Spectra computed previously for {bids_stuff}. If you want to run again, delete the old stuff, chap')
    else: #calculate spectra
        inout.if_not_exist_make(f'{task_dir}')
        for spec in ['powerspectra', 'energyspectra','reconspectra']:
            inout.if_not_exist_make(f'{task_dir}/{spec}')
            inout.if_not_exist_make(f'{task_dir}/criticality')
            inout.if_not_exist_make(f'{task_dir}/criticality/{spec}')
        #power spectra
        print(f'[CHAP] Computing mean, dynamic, and normalized power spectra for {bids_stuff}...')
        mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs) #average power over the whole scan (average of dynamic for each harmonic)
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_mean_power_spectrum', mean_power_spectrum)
        dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals) #power at each TR
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_dynamic_power_spectrum', dynamic_power_spectrum)
        normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
        np.save(f'{task_dir}/powerspectra/{bids_stuff}_normalized_power_spectrum', normalized_power_spectrum)
        print(f'[CHAP] Saved power spectra')
        #energy spectra
        print(f'[CHAP] Computing mean and dynamic energy spectra for {bids_stuff}...')
        mean_energy_spectrum = cs.mean_energy_spectrum(timeseries, vecs, vals) #average energy over the whole scan (average of dynamic for each harmonic)
        np.save(f'{task_dir}/energyspectra/{bids_stuff}_mean_energy_spectrum', mean_energy_spectrum)
        dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals) #energy at each TR
        np.save(f'{task_dir}/energyspectra/{bids_stuff}_dynamic_energy_spectrum', dynamic_energy_spectrum)
        print(f'[CHAP] Saved energy spectra')
        #reconstruction spectrum
        dynamic_reconstruction_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals) #takes on negative values
        np.save(f'{task_dir}/reconspectra/{bids_stuff}_dynamic_reconstruction_spectrum', dynamic_reconstruction_spectrum)
        print(f'[CHAP] Saved dynamic reconstruction spectrum for {bids_stuff} scan')
        #criticality
        if args.criticality == True:
            inout.if_not_exist_make(f'{task_dir}/criticality/')
            power_criticality = cs.criticality(dynamic_power_spectrum, 'power')
            power_criticality.to_csv(f'{task_dir}/criticality/power/power_criticality.csv', index=False)
            energy_criticality = cs.criticality(dynamic_energy_spectrum, 'energy')
            energy_criticality.to_csv(f'{task_dir}/criticality/power/energy_criticality.csv', index=False)
            #TODO: reconstruction spectrum has no 'mean_recon_spectra' - should we still be calcluating criticality for this spectrum
            recon_criticality = cs.criticality(dynamic_reconstruction_spectrum, 'reconstruction')
            recon_criticality.to_csv(f'{task_dir}/criticality/power/reconstruction_criticality.csv', index=False)

      
    
  
    
  
    
  
                
                
                