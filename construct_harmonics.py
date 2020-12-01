#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 16:01:27 2020
used in CHAP entrypoint_script

    sc - array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    ec - array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]
    tol - search radius of nearest neighbor search for matching endpoints to surface vertices
    NNnum - number of nearest neighboring surface vertices to assign to each endpoint
"""
import decomp as dcp
import input_output as inout
import subprocess
import numpy as np
import os
import matrix_methods as mm
import compute_spectra as cs
from scipy import sparse

#this fxn isn't finished, fix user_info session stuff
def prep_harmonics_bids(args, sub, file, user_info, multises, ses=""):
    tck_name = file.split('/')[-1][:-4]
    print('[CHAP] Saving streamline endpoints and converting to vtk...')
    subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
    user_info[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk'
    for file in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}'):
        if '_endpoints.tck' in file:
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
    print('[CHAP] Finished MRtrix commands')
    #construct surface coordinates, surface endpoints
    user_info[f'{sub}_info'][ses]['surfs']['lh'] = f'{args.surf_dir}/sub-{sub}/surf/lh.white.corresponded.vtk'
    user_info[f'{sub}_info'][ses]['surfs']['rh'] = f'{args.surf_dir}/sub-{sub}/surf/rh.white.corresponded.vtk'
      

def construct_harmonics_calculate_spectra(args, sub, ses, user_info, multises): 
    if 'vtk' in user_info[f'{sub}_info'][ses]['surfs']['lh']:
        sc,si=inout.read_vtk_surface_both_hem(user_info[f'{sub}_info'][ses]['surfs']['lh'], user_info[f'{sub}_info'][ses]['surfs']['rh'])
    else:
        sc,si=inout.read_gifti_surface_both_hem(user_info[f'{sub}_info'][ses]['surfs']['lh'], user_info[f'{sub}_info'][ses]['surfs']['rh'])
    print('[CHAP] Saved surface coordinates and surface indices')
    ec=inout.read_streamline_endpoints(user_info[f'{sub}_info'][ses]['endpoints']) #read endpoint locations into numpy array
    print('[CHAP] Saved endpoint coordinates')
    print('[CHAP] Constructing surface matrix...')
    surf_mat=mm.construct_surface_matrix(sc,si) #construct surface matrix from sc and si
    print(f'shapes: ec is {ec.shape}, sc is {sc.shape}, si is {si.shape}')
    print('[CHAP] Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc, ec, tol=3, NNnum = args.nnum) #construct struc conn matrix from ec and sc 
    sparse.save_npz(f'{args.output_dir}/chap/sub-{sub}/{ses}/struc_conn_mat', struc_conn_mat) #save structural connectivity matrix
    print('[CHAP] Saved structural connectivity matrix')
    print('[CHAP] Computing harmonics...')
    vals,vecs=dcp.lapDecomp(struc_conn_mat, args.evecs) #laplacian decomposition
    os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis') #visualization output directory
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vals',vals)
    np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs',vecs)
    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,vecs) 
        print(f'[CHAP] Saved harmonics for {sub} {ses}')
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,vecs)
        print(f'[CHAP] Saved harmonics for {sub}')
    if args.fprep_dir: #if functional images are specified
        if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/func'):
            os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/func') #func output folder
        for vol in user_info[f'{sub}_info']['func']:
            if ses in vol:
                task = inout.get_task(vol) #get taskname
                full_path_lh = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/surfmapped_vol_lh.func.gii'
                full_path_rh = f'{args.output_dir}/chap/sub-{sub}/{ses}/func/surfmapped_vol_rh.func.gii'           
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
                bids_stuff = inout.get_bids_stuff(full_path_lh)
                #read functional timeseries of surface mapped volume
                timeseries = cs.read_functional_timeseries(full_path_lh, full_path_rh)
                #power spectra
                if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra'):
                    os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra')
                mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_mean_power_spectrum', mean_power_spectrum)
                dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_dynamic_power_spectrum', dynamic_power_spectrum)
                normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/powerspectra/{bids_stuff}_normalized_power_spectrum', normalized_power_spectrum)
                print('[CHAP] Computed power spectra')
                #energy spectra
                if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra'):
                    os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra')
                mean_energy_spectrum = cs.mean_energy_spectrum(timeseries, vecs, vals)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra/{bids_stuff}_mean_energy_spectrum', mean_energy_spectrum)
                dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/energyspectra/{bids_stuff}_dynamic_energy_spectrum', dynamic_energy_spectrum)
                print('[CHAP] Computed energy spectra')
                #reconstruction spectrum
                if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/reconspectra'):
                    os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/reconspectra')
                dynamic_reconstruction_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals)
                np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/reconspectra/{bids_stuff}_dynamic_reconstruction_spectrum', dynamic_reconstruction_spectrum)
                print('[CHAP] Computed reconstruction spectra')
    print(f'[CHAP] Finished {ses}')
                
                

                
                
                
                
                
                