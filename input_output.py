#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:16:15 2020

@author: patricktaylor
"""
import numpy as np
import vtk
import meshio
#from tvtk.api import tvtk, write_data
import nibabel as nib
import test_retest_fxns as t_rt
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler,normalize
from nilearn import plotting
import matplotlib.pyplot as plt
'''
def save_surface(filename,points,edges,feature=None):
    mesh = tvtk.PolyData(points=points, polys=edges)
    if feature is not None:
        mesh.point_data.scalars=feature
    write_data(mesh, filename)
    return
'''
'''
mask=generate_mask_from_parc()
m= conn
masked_connectivity=utility.mask_connectivity_matrix(m,mask)

vals,mvecs=lap_decomp(masked_connectivity,n)

vecs=util.unmask_medial_wall_vecs(masked_connectivity,mask)
'''


#ts=read_functional_timeseries()

#maskedts=uts.mask_timeseries(ts, mask)

def generate_mask_from_parc_hcp(lhparc,rhparc):
    #105923/MNINonLinear/fsaverage_LR32k/105923.L.aparc.32k_fs_LR.label.gii'105923/MNINonLinear/fsaverage_LR32k/105923.R.aparc.32k_fs_LR.label.gii'
    parc=read_gifti_feature_both_hem(lhparc,rhparc).astype('int32')
    inds1=np.where(parc==-1)[0]
    mask=np.zeros(len(parc))
    mask[inds1]=1
    return mask 


def save_eigenvector(filename,points,edges,vecs):
    Cells = {"triangle": edges}
    V={}
    for i in range (0,len(vecs[0,:])):
        V.update({"ev%d" % (i):vecs[:,i]})
    mesh=meshio.Mesh(points,Cells,V)
    meshio.write(filename,mesh)
    return

def save_eigenvector_to_hems(filename,points,edges,vecs):
    half=int(len(points)/2)
    lhsc=points[:half]
    rhsc=points[half:]
    lhsi=edges[:int(len(edges)/2)]
    rhsi=edges[int(len(edges)/2):]-half
    lhvec=vecs[:half]
    rhvec=vecs[half:]
    save_eigenvector(filename %'lh',lhsc,lhsi,lhvec)
    save_eigenvector(filename %'rh',rhsc,rhsi,rhvec)
    return


def read_vtk_feature(filename,featurename):
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    scalars=data.GetPointData()
    Scalars=scalars.GetArray(featurename)
    feature=np.zeros((data.GetNumberOfPoints(),))
    for i in range (data.GetNumberOfPoints()):
        feature[i]=Scalars.GetValue(i)
    return feature

def read_vtk_feature_both_hem(lfile,rfile,featurename):
    l=read_vtk_feature(lfile,featurename)
    r=read_vtk_feature(rfile,featurename)
    feature=np.hstack((l,r))
    return feature


def read_vtk_surface(filename):
    #reads a vtk surface mesh and returns the coordinates of vertices (nvert,3), and the connections definiing the mesh (ncon,3) as numpy arrays
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    CellArray = data.GetPolys()
    Polygons = CellArray.GetData()
    edges=np.zeros((CellArray.GetNumberOfCells(),3))
    for i in range (0, CellArray.GetNumberOfCells()):
        edges[i,:]=[Polygons.GetValue(j) for j in range (i*4+1,i*4+4)]
    points=np.zeros((data.GetNumberOfPoints(),3))
    for i in range(data.GetNumberOfPoints()):
            points[i,:] = data.GetPoint(i)
    return points, edges

def read_vtk_surface_both_hem(Lfile,Rfile):
    #takes full path+file name of two hemispheres of a surface, loads them, and combines the coordinates and connections from both into 2 full surface arrays
    lhc,lhi=read_vtk_surface(Lfile)
    rhc,rhi=read_vtk_surface(Rfile)
    coords,si=combine_hemis(lhc,rhc,lhi,rhi)
    return coords,si

def read_gifti_surface(filename,hcp=False):
    if hcp:
        data=nib.load(filename)
        points=data.darrays[0].data
        edges=data.darrays[1].data
        return points,edges
    else:
        data=nib.load(filename)
        points=data.darrays[1].data
        edges=data.darrays[0].data
        return points,edges

def read_gifti_surface_both_hem(Lfile,Rfile,hcp=False):
    lhc,lhi=read_gifti_surface(Lfile,hcp)
    rhc,rhi=read_gifti_surface(Rfile,hcp)
    points,edges=combine_hemis(lhc,rhc,lhi,rhi)
    return points,edges

def combine_hemis(lhc,rhc,lhi,rhi):
    #concatenates surface coordinates of two hemispheres and creates connectivity array for full surface
    coords=np.vstack((lhc,rhc))
    si=np.vstack((lhi,rhi+len(lhc)))
    return coords, si

def read_functional_timeseries(lhfunc,rhfunc,bcp=True):
    l = nib.load(lhfunc).darrays
    r = nib.load(rhfunc).darrays
    timeseries = np.zeros((2*len(l[0].data), len(r)))
    if bcp:
        timeseries=np.concatenate((np.array(l[0].data),np.array(r[0].data)))
        return timeseries
    else:
        for i in range(len(l)):
            lt = np.array(l[i].data)
            rt = np.array(r[i].data)
            tp = np.concatenate((lt, rt))
            timeseries[:, i] = tp
        return timeseries



def read_streamline_endpoints(filename):
    #reads endpoint locations of vtk file containing only the endpoints of a tractogram. returns numpy array of size (nEndpoints,3).
    #endpoint 0 and endpoint 1 correspond to the same fiber. endpoint 2, endpoint 3 correspond to the same fiber... etc
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    points=np.zeros((data.GetNumberOfPoints(),3))
    for i in range(data.GetNumberOfPoints()):
            points[i,:] = data.GetPoint(i)
    return points


'''
def split_vtk_feature_to_hems(path,fname):
    sc,si=read_vtk_surface(path+fname)
    feature=read_vtk_feature(path+fname,'scalars')
    lsc=sc[:int(len(sc)/2)]
    rsc=sc[int(len(sc)/2):]
    lsi=si[:int(len(si)/2)]
    rsi=si[int(len(si)/2):]-int(len(sc)/2)
    save_surface(path+'lh_'+fname,lsc,lsi,feature[:int(len(sc)/2)])
    save_surface(path+'rh_'+fname,rsc,rsi,feature[int(len(sc)/2):])
    return
'''

def gifti_to_scalar(L,R):
    l=L.darrays
    r=R.darrays
    La=np.array([l[0].data]).T
    Ra=np.array([r[0].data]).T
    scal=np.vstack((La,Ra))
    return scal

def read_gifti_feature_both_hem(lfname,rfname):
    L=nib.load(lfname)
    R=nib.load(rfname)
    featurevec=gifti_to_scalar(L,R)
    return featurevec

def generate_mask_from_parc(lhparc,rhparc):
    parc=read_vtk_feature_both_hem(lhparc,rhparc).astype('int32')
    inds1=np.where(parc==1639705)[0]
    #inds3=np.where(parc==1639704)[0]
    inds2=np.where(parc==3294840)[0]
    #inds4=np.where(parc==3294839)[0]
    mask=np.zeros(len(parc))
    mask[inds1]=1
    mask[inds2]=1
    #mask[inds3]=1
    #mask[inds4]=1
    return mask

def get_run(fname):
    runstart = fname.find('run') + 4
    x = fname[runstart:]
    run = x.split('_')[0]
    return str(run)

def get_acq(fname):
    acqstart = fname.find('acq') + 4
    x = fname[acqstart:]
    acq = x.split('_')[0]
    return acq

def get_task(fname):
    taskstart = fname.find('task') + 5
    x = fname[taskstart:]
    task = x.split('_')[0]
    return task
    
def get_bids_stuff(lh_full_path):
    stuff_start = lh_full_path.find('vol_') + 4
    x = lh_full_path[stuff_start:]
    bids_stuff = x.split('lh')[0][:-1]
    return(bids_stuff)

def if_not_exist_make(path):
    import os
    if not os.path.exists(path):
        os.mkdir(path)

def normalize_ts(ts):
    ts=ts.T
    nts=(ts - np.mean(ts, axis=0)) / np.std(ts, axis=0)
    return nts.T

def combine_pe(ts_lr, ts_rl):
    ts_lr_n = normalize_ts(ts_lr)
    ts_rl_n = normalize_ts(ts_rl)
    return np.hstack((ts_lr_n, ts_rl_n))

def normalize_ev(ev):
    ev = ev.reshape(1,-1)
    ev = ev.T
    scaler = MinMaxScaler(feature_range=(-1,1))
    ev = scaler.fit_transform(ev)
    return ev.T

def scale_ev(ev):
    return np.interp(ev, (ev.min(), ev.max()), (-1, +1))

def network_verts(network, parcel_csv, dtseries):
    network_parcels = np.where(parcel_csv['Community']==network)[0]
    network_parcels = [p+1 for p in network_parcels]
    dtseries = [int(p) for p in dtseries]
    network_verts = []
    for p in dtseries:
        if p in network_parcels:
            network_verts.append(1)
        else:
            network_verts.append(0)
    return np.array(network_verts)
#usage: if only one session, input 'ses' as 'single', if multises, use session number
def visualize_harmonics(sub, ses, si, sc, vecs, output_dir):   
    for i in range(vecs[0].shape[0]):
        fig,ax = plt.subplots(2, 2, subplot_kw={'projection': '3d'}, dpi=1000)
        fig.suptitle(f'ev_{i} {sub} {ses}')
        plot = plotting.plot_surf_stat_map([sc,si],vecs[:,i],view='dorsal',cmap='RdBu',output_file=f'{output_dir}ev_{i}_ses_{ses}.png', colorbar=False,vmax=.005, axes=ax[0][0], figure=fig)
        plot = plotting.plot_surf_stat_map([sc,si],vecs[:,i],view='medial',cmap='RdBu',output_file=f'{output_dir}ev_{i}_ses_{ses}.png', colorbar=False,vmax=.005, axes=ax[0][1], figure=fig)
        plot = plotting.plot_surf_stat_map([sc[:int(sc.shape[0]/2)],si[:int(si.shape[0]/2)]],vecs[:int(vecs.shape[0]/2),i],view='medial',cmap='RdBu',output_file=f'{output_dir}ev_{i}_ses_{ses}.png', colorbar=False,vmax=.005, axes=ax[1][0], figure=fig)
        plot = plotting.plot_surf_stat_map([sc,si],vecs[:,i],view='lateral',cmap='RdBu',output_file=f'{output_dir}ev_{i}_ses_{ses}.png', colorbar=False,vmax=.005, axes=ax[1][1], figure=fig)
        plt.close(fig)


'''
parcel_csv = pd.read_csv('/Users/bwinston/Downloads/Parcels/Parcels.csv')
dtseries = np.array(np.loadtxt('/Users/bwinston/Downloads/Gordon_Parcels_LR.dtseries.txt'))
dtseries = np.expand_dims(dtseries,1)
masked_vecs = np.load('/Users/bwinston/Downloads/chap_out_test/sub-105923/ses-retest/vecs.npy')
masked_vecs = np.delete(masked_vecs,0,axis=1)
net_verts = {}
for network in list(set(parcel_csv['Community'])):
    net_verts[network] = {} 
    net_verts[network]['verts'] = network_verts(network, parcel_csv, dtseries)
    net_verts[network]['corrs'] = []
    for i in range(0,99):
       net_verts[network]['corrs'].append(abs(pearsonr(net_verts[network]['verts'],masked_vecs[:,i])[0])) 
    
    
    

def read_gifti_surface(filename):
    data=nib.load(filename)
    points=data.darrays[1].data
    edges=data.darrays[0].data
    return points,edges

def read_gifti_surface_both_hem(Lfile,Rfile):
    lhc,lhi=read_gifti_surface(Lfile)
    rhc,rhi=read_gifti_surface(Rfile)
    points,edges=combine_hemis(lhc,rhc,lhi,rhi)
    return points,edges   
    
def combine_hemis(lhc,rhc,lhi,rhi):
    #concatenates surface coordinates of two hemispheres and creates connectivity array for full surface
    coords=np.vstack((lhc,rhc))
    si=np.vstack((lhi,rhi+len(rhc)))
    return coords, si
'''
    
    
    
    
