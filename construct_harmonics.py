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

def construct_harmonics_calculate_spectra(args, sub, output_dir, file, multises, ses=""):
    tck_name = file.split('/')[-1][:-4]
    print('[CHAP] Saving streamline endpoints and converting to vtk...')
    subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
    for file in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}'):
        if '_endpoints.tck' in file:
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
    print('[CHAP] Finished MRtrix commands')
    #construct surface coordinates, surface endpoints
    lh_surf_path = f'{args.surf_dir}/sub-{sub}/surf/lh.white.corresponded.vtk'
    rh_surf_path = f'{args.surf_dir}/sub-{sub}/surf/rh.white.corresponded.vtk'
    if lh_surf_path.endswith('.vtk'):
        sc,si=inout.read_vtk_surface_both_hem(lh_surf_path, rh_surf_path)
        print('[CHAP] Saved surface coordinates and surface indices')
    else:
        sc,si=inout.read_gifti_surface_both_hem(lh_surf_path, rh_surf_path)
        print('[CHAP] Saved surface coordinates and surface indices')
    streamline_path = f'{output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk'
    ec=inout.read_streamline_endpoints(streamline_path) #read endpoint locations into numpy array
    print('[CHAP] Saved endpoint coordinates')
    print('[CHAP] Constructing surface matrix...')
    surf_mat=mm.construct_surface_matrix(sc,si) #construct surface matrix from sc and si
    ihc_mat=mm.construct_inter_hemi_matrix(sc,tol=3)
    print('[CHAP] Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc, ec, tol=3, NNnum = args.nnum) #construct struc conn matrix from ec and sc 
    sparse.save_npz(f'{output_dir}/chap/sub-{sub}/{ses}/struc_conn_mat', struc_conn_mat) #save structural connectivity matrix
    print('[CHAP] Saved structural connectivity matrix')
    print('[CHAP] Computing harmonics...')
    vals,vecs=dcp.lapDecomp(struc_conn_mat, args.evecs) #laplacian decomposition
    os.mkdir(f'{output_dir}/chap/sub-{sub}/{ses}/vis') #visualization output directory
    np.save(f'{output_dir}/chap/sub-{sub}/{ses}/vals',vals)
    np.save(f'{output_dir}/chap/sub-{sub}/{ses}/vecs',vecs)
    if multises:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,vecs) 
    else:
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_harmonics.vtk',sc,si,vecs)
    print('[CHAP] Saved harmonics for {ses}')
    #Compute spectra as specified
    #TODO: add correct filepaths once volume-to-surface mapping is complete
    if args.fprep_dir:
        full_path_lh = "placeholder_lh.gii"
        full_path_rh = "placeholder_rh.gii"
        timeseries = cs.read_functional_timeseries(full_path_lh, full_path_rh)
        os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'powerspectra')
        mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs)
        dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals)
        normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'powerspectra', mean_power_spectrum)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'powerspectra', dynamic_power_spectrum)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'owerspectra', normalized_power_spectrum)
        os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'energyspectra')
        dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals)
        normalized_energy_spectrum = cs.normalized_energy_spectrum(timeseries, vecs)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'powerspectra', mean_energy_spectrum)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'powerspectra', dynamic_energy_spectrum)
        os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'reconspectra')
        recon_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals)
        np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'reconspectra', dynamic_reconstruction_spectrum)
