"""
@author: patricktaylor
"""

import numpy as np
import nibabel as nib


def dynamic_energy_spectrum(timeseries,vecs,vals):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros((len(vecs[0,:]),len(timeseries[0,:])))
    for k in range(len(spectrum)):
        v=vecs[:,k]
        for tp in range(len(zeromean[0,:])):
            spectrum[k][tp]=(np.abs(np.dot(v,zeromean[:,tp]))**2)*vals[k]**2
    return spectrum 

def dynamic_power_spectrum(timeseries,vecs,vals):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros((len(vecs[0,:]),len(timeseries[0,:])))
    for k in range(len(spectrum)):
        #print(k)
        v=vecs[:,k]
        for tp in range(len(zeromean[0,:])):
            spectrum[k][tp]=np.abs(np.dot(v,zeromean[:,tp]))
    return spectrum

def dynamic_reconstruction_spectrum(timeseries,vecs,vals):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros((len(vecs[0,:]),len(timeseries[0,:])))
    for k in range(len(spectrum)):
        #print(k)
        v=vecs[:,k]
        for tp in range(len(zeromean[0,:])):
            spectrum[k][tp]=np.dot(v,zeromean[:,tp])
    return spectrum

def mean_energy_spectrum(timeseries,vecs,vals):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros(len(vecs[0,:]))
    for k in range(len(spectrum)):
      v=vecs[:,k]
      for tp in range(len(zeromean[0,:])):
          spectrum[k]+=np.abs(np.dot(v,zeromean[:,tp]))**2/len(zeromean[0,:])*vals[k]**2
    return spectrum 
def mean_power_spectrum(timeseries,vecs):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros(len(vecs[0,:]))
    for k in range(len(spectrum)):
        #print(k)
        v=vecs[:,k]
        for tp in range(len(zeromean[0,:])):
            spectrum[k]+=np.abs(np.dot(v,zeromean[:,tp]))/len(zeromean[0,:])
    return spectrum 

def normalized_power_spectrum(timeseries,vecs):
    zeromean = np.zeros(np.shape(timeseries))
    for i in range(len(timeseries)):
        zeromean[i, :] = (timeseries[i, :] - np.mean(timeseries[i, :]))
    spectrum=np.zeros(len(vecs[0,:]))
    for k in range(len(spectrum)):
        #print(k)
        v=vecs[:,k]
        for tp in range(len(zeromean[0,:])):
            spectrum[k]+=np.abs(np.dot(v,zeromean[:,tp]))/len(zeromean[0,:])/np.linalg.norm(v)/np.linalg.norm(zeromean[:,tp])
    return spectrum

def read_functional_timeseries(lhfunc,rhfunc):
    l = nib.load(lhfunc).darrays
    r = nib.load(rhfunc).darrays
    timeseries = np.zeros((2*len(l[0].data), len(r)))
    for i in range(len(l)):
        lt = np.array(l[i].data)
        rt = np.array(r[i].data)
        tp = np.concatenate((lt, rt))
        timeseries[:, i] = tp
    return timeseries   

    
    
    