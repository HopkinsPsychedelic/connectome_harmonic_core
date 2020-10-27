#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 16:01:27 2020

@author: quintinfrerichs
"""
import numpy as np
import matrix_methods as mm
import compute_spectra as cs

def construct_harmonics_calculate_spectra(output_dir, file, ses=""):
    tck_name = file.split('/')[-1][:-4]
    os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'endpoints')
    subprocess.check_call("./mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/'+ses+'dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/'+ses+'endpoints'), shell=True)
    #construct surface coordinates, surface endpoints
    lh_surf_path =  user_info[f'{sub}_info']['lh_surf'][0]
    rh_surf_path = user_info[f'{sub}_info']['rh_surf'][0]
    if lh_surf_path.endswith('.vtk'):
        sc,si=inout.read_vtk_surface_both_hem(lh_surf_path, rh_surf_path)
    else:
        sc,si=inout.read_gifti_surface_both_hem(lh_surf_path, rh_surf_path)
    streamline_path = f'{output_dir}/chap/sub-{sub}/'+ses+'endpoints/{tck_name}.vtk'
    ec=inout.read_streamline_endpoints(streamline_path)
    print('Constructing surface matrix...')
    surf_mat=mm.construct_surface_matrix(sc,si)
    ihc_mat=mm.construct_inter_hemi_matrix(sc,tol=3)
    print('Constructing structural connectivity matrix...')
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc,ec,tol=3,NNnum=45)
    print('Computing harmonics...')
    vals,vecs=dcp.lapDecomp(struc_conn_mat,args.number)
    os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'vecsvals')
    os.mkdir(f'{output_dir}/chap/sub-{sub}/'+ses+'vis')
    np.save(f'{output_dir}/chap/sub-{sub}/'+ses+'vecsvals/',[vals,vecs])
    inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/'+ses+'vis/',sc,si,vecs)
    #Compute spectra as specified
    #TODO: add correct filepaths once volume-to-surface mapping is complete
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