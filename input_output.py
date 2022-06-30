#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:16:15 2020

@author: patricktaylor
"""
import numpy as np
#import vtk
import meshio
#from tvtk.api import tvtk, write_data
import nibabel as nib
#import test_retest_fxns as t_rt
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler,normalize
import utility_functions as uts
import os
from glob import glob
import statistics as stats
import pickle
from scipy import sparse
import pandas as pd
import test_retest_fxns as t_rt
import datetime
import icc as ic
import json

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
    fname = fname.split('/')[-1]
    taskstart = fname.find('task') + 5
    x = fname[taskstart:]
    task = x.split('_')[0]
    return task

def add_bids_thing_to_fname(bids_thing,vol,full_path_lh,full_path_rh):
    if bids_thing in vol:
        func = 'get_{bids_thing}'
        thing = eval(func+'(vol)')
        full_path_lh = full_path_lh[:-11] + f'{bids_thing}-{thing}_' + full_path_lh[-11:]
        full_path_rh = full_path_rh[:-11] + f'{bids_thing}-{thing}_' + full_path_rh[-11:]
        return full_path_lh,full_path_rh
    
def get_bids_stuff(dts):
    dts = dts.split('/')[-1] #get just filename not path
    stuff_end = dts.find('desc') - 1 #get everything before desc
    return(dts[:stuff_end])

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


def net_verts():
    parcel_csv = pd.read_csv('/data2/Brian/connectome_harmonics/Parcels/Parcels.csv')
    dtseries = np.array(np.loadtxt('/data2/Brian/connectome_harmonics/Gordon_Parcels_LR.dtseries.txt'))
    dtseries = np.expand_dims(dtseries,1)
    masked_vecs = np.load('/data/hcp_test_retest_pp/derivatives/chap/sub-114823/ses-test/vecs.npy')
    masked_vecs = np.delete(masked_vecs,0,axis=1)
    net_verts = {}
    for network in list(set(parcel_csv['Community'])):
        net_verts[network] = {} 
        net_verts[network]['verts'] = network_verts(network, parcel_csv, dtseries)
        net_verts[network]['unmasked_verts'] = uts.unmask_medial_wall(net_verts[network]['verts'],np.load('/data2/Brian/connectome_harmonics/mask.npy'))
    return net_verts
'''
#NET VERTS mac
parcel_csv = pd.read_csv('/Users/bwinston/Downloads/Parcels/Parcels.csv')
dtseries = np.array(np.loadtxt('/Users/bwinston/Downloads/Gordon_Parcels_LR.dtseries.txt'))
dtseries = np.expand_dims(dtseries,1)
masked_vecs = np.load('/Users/bwinston/Downloads/chap_out_test/sub-105923/ses-retest/vecs.npy')
masked_vecs = np.delete(masked_vecs,0,axis=1)
net_verts = {}
for network in list(set(parcel_csv['Community'])):
    net_verts[network] = {} 
    net_verts[network]['verts'] = network_verts(network, parcel_csv, dtseries)
    net_verts[network]['unmasked_verts'] = uts.unmask_medial_wall(net_verts[network]['verts'],np.load('/Users/bwinston/Documents/connectome_harmonics/hcp_mask.npy'))
'''

'''
627549 missing resting state
859671 retest structural harmonics?
'''
def get_subs(chap_dir,functional=False, rest = False, t_rt=False):
   subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
   subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
   if t_rt==False:
       for sub in ['test_avg', 'retest_avg', 'total_avg','114823','115320','139839','172332','192439','185542','185442','859671']: #add bad subs here 859671 just there bc other pipeline
            if os.path.exists(f'{chap_dir}/sub-{sub}'):
                subs.remove(sub)
   else:
       for sub in ['test_avg', 'retest_avg', 'total_avg','859671','187547','662551','433839']: #add bad subs here
            if os.path.exists(f'{chap_dir}/sub-{sub}'):
                subs.remove(sub)        
   if functional == True:
        subs.remove('341834')
        subs.remove('627549')
   elif rest == True:
       subs = [sub for sub in subs if sub not in ['187547','341834','859671','627549','177746']]
       #subs = [sub for sub in subs if sub not in ['200109','146129','783462','187547','287248','158035','111312','660951','341834','859671','135528']]
   #subs = ['105923','103818','111312']
   return subs

def mofl(list_of_lists):
    return np.mean(np.array(list_of_lists),axis=0)

'''
def across_avg(subs,av,dic,fxn,data,mofl=True): #dic doesn't have to be overall dict #data is just name
    dic[f'across_subj_all_{data}'] = []
    for sub in subs:
        dic[sub][f'c_sub_all_{data}'] = []
        for c_sub in subs:
            if c_sub != sub:
                if sub not in dic[c_sub]:
                dic[sub][c_sub] = {}
                for ses in ['test','retest']: #take mofl out of below if isn't mofl #gd needs mofl
                    if mofl:
                        dic[sub][c_sub][ses] = fxn(mofl(dic[sub][ses][f'{data}']),mofl(dic[c_sub][ses][f'{data}']))
                    else:
                        dic[sub][c_sub][ses] = fxn(av[sub][ses]['vecs'],av[c_sub][ses]['vecs'],)
                dic[sub][f'c_sub_all_{data}'].append((dic[sub][c_sub]['test'] + dic[sub][c_sub]['retest'])/2)
        dic[f'across_subj_all_{data}'].append(stats.mean(dic[sub][f'c_sub_all_{data}']))
        dic[f'across_subj_all_{data}'] = list(set(dic[f'across_subj_all_{data}']))
    dic[f'across_subj_avg_{data}'] = stats.mean(dic[sub][f'c_sub_all_{data}'])
'''
def across_avg(subs,dic,fxn,data,mofl=True,seslist = ['test','retest']): #dic doesn't have to be overall dict #data is just name
    dic[f'across_subj_all_{data}'] = []
    for sub in subs:
        dic[sub][f'c_sub_all_{data}'] = []
        for c_sub in subs:
            if c_sub != sub:
                if sub not in dic[c_sub]:
                    dic[sub][c_sub] = {}
                    for ses in seslist: #take mofl out of below if isn't mofl #gd needs mofl
                        if mofl:
                            dic[sub][c_sub][ses] = fxn(mofl(dic[sub][ses][f'{data}']),mofl(dic[c_sub][ses][f'{data}']))
                        else:
                            dic[sub][c_sub][ses] = fxn(dic[sub][ses][f'{data}'],dic[c_sub][ses][f'{data}'])
                    dic[sub][f'c_sub_all_{data}'].append((dic[sub][c_sub]['test'] + dic[sub][c_sub]['retest'])/2)
        dic[f'across_subj_all_{data}'].append(dic[sub][f'c_sub_all_{data}'])
    dic[f'across_subj_avg_{data}'] = stats.mean(sum(dic[f'across_subj_all_{data}'],[]))

                
def abs_pearson(x,y,fisher=False, abso = True):
    if fisher==True and abso==True:
        return np.arctanh(abs(pearsonr(x,y)[0]))
    elif fisher==False and abso==True:
        return abs(pearsonr(x,y)[0])
    elif fisher==True and abso==False:
        return np.arctanh(pearsonr(x,y)[0])
    elif fisher==False and abso==False:
        return pearsonr(x,y)[0]

def load_pkl(file): #file should end in .pkl
    file_to_read = open(file,"rb")
    loaded_dic = pickle.load(file_to_read)
    return loaded_dic

'''
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
    
def dts_to_func_gii(dts,out):
    os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {dts} COLUMN -metric CORTEX_LEFT {out}_hem-l.func.gii')
    os.system(f'bash /home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -cifti-separate {dts} COLUMN -metric CORTEX_RIGHT {out}_hem-r.func.gii') 
    
def get_rdists(rdist,how_many):
    rdist_dic = {}
    for x in range(how_many):
        rdist_dic[x] = rdist[:,x]
    return rdist_dic

def save_pickle(a_dict,fname):
    with open(fname,'wb') as file:
        pickle.dump(a_dict,file,protocol=pickle.HIGHEST_PROTOCOL)
        
def load_pickle(fname):
    with open(fname,'rb') as file:
        return pickle.load(file)
    
def save_degree_vector(sub):
    mat = sparse.load_npz(f'/data/HCP_Raw/derivatives/chap/sub-{sub}/struc_conn_mat.npz')
    mat = mat.toarray()
    jeff = []
    for vert in range(len(mat)):
        jeff.append(sum(mat[vert]))
    sums = np.array(jeff)
    sums_unmasked = uts.unmask_medial_wall_vecs(sums,'/usr/local/connectome_harmonic_core/connectome_harmonic_core/hcp_mask.npy')
    return sums

def get_sc(sub): #gets HCP_Raw surface info
    sc = {}
    sc['sc'],sc['si'] = read_gifti_surface_both_hem(f'/data/HCP_Raw/derivatives/ciftify/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.L.white.32k_fs_LR.surf.gii',f'/data/HCP_Raw/derivatives/ciftify/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.R.white.32k_fs_LR.surf.gii', True)
    sc['lhc'],sc['lhi'] = read_gifti_surface(f'/data/HCP_Raw/derivatives/ciftify/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.L.white.32k_fs_LR.surf.gii',hcp=True)
    sc['rhc'],sc['rhi'] = read_gifti_surface(f'/data/HCP_Raw/derivatives/ciftify/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.R.white.32k_fs_LR.surf.gii',hcp=True)
    return sc 

def icc_vecs(vec1,vec2, icc_module=True):
    begin_time = datetime.datetime.now()
    vec2 = t_rt.check_polarity(vec1,vec2)
    if icc_module == True:
        measures = pd.DataFrame({'vec1':vec1, 'vec2':vec2})
        return ic.icc(measures,model='twoway',type='agreement',unit='single')[0]
    target_59 = list(range(59412))
    target = target_59 + list(range(59412))
    rater = ['vec1'] * 59412
    rater = rater + (['vec2'] * 59412)
    rating = np.hstack((vec1,vec2))
    df = pd.DataFrame({'vertex':target, 'vecs':rater, 'value':rating})
    icc_time = datetime.datetime.now()
    icc = intraclass_corr(data=df, targets = 'vertex', raters = 'vecs', ratings = 'value')
    icc.set_index('Type')
    print(f'icc itself took {datetime.datetime.now() - icc_time} h:m:s')
    print(f'total time {datetime.datetime.now() - begin_time} h:m:s')
    return icc['ICC'][0]

def rem_json_field(fname,field):                                                               
    obj  = json.load(open(fname))
    try:
        del(obj[field])
        open(fname, "w").write(
            json.dumps(obj, sort_keys=True, indent=4, separators=(',', ': ')))
    except:
        print('no volume timing')

