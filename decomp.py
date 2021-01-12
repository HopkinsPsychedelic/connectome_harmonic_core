#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 19:36:03 2020

@author: patricktaylor
"""
import scipy.sparse.linalg
from scipy.sparse import csgraph
import time
import numpy as np

from sklearn.decomposition import PCA

def get_group_pca_comp(evlist,num):
    tempmat = np.zeros((np.shape(evlist[0])[1]-1,num*len(evlist)))
    for i in range (len(evlist)):
        tempmat[:,i*num:(i+1)*num]=evlist[i][:,1:num] #from 0-10, 10-20, 20-30   
    pca = PCA(n_components=num)
    pca.fit(tempmat.T)
    components = pca.components_
    components = components.T
    return components 

def get_group_pca_comp_b(evlist,num):
    tempmat = np.zeros((np.shape(evlist[0])[0],num*len(evlist)))
    for i in range (len(evlist)):
        tempmat[:,i*num:(i+1)*num]=evlist[i][:,0:num] #from 0-10, 10-20, 20-30   
    pca = PCA(n_components=num)
    pca.fit(tempmat.T)
    components = pca.components_
    components = components.T
    return components 

def lapDecomp(Asparse,num):
    #Asparse- csr format sparse adjacency matrix
    #num - number of laplacian eigenmodes to compute
    #returns tuple of eigenvalues (vals) and eigenvectors (vecs) as arrays of size (dim(Asparse),1) and  size (dim(Asparse),num)
    start = time.time()
    LA=csgraph.laplacian(Asparse,normed=True)
    L=LA.astype('float64')
    vals,vecs=scipy.sparse.linalg.eigsh(L,num,which='SM')
    end = time.time()
    print('decomposition time=', (end-start),'seconds, which is', (end - start)/60, 'minutes')
    return vals,vecs

def cos_norm(v1,v2):
    
    dot = np.dot(v1,v2)
    
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    
    return np.abs(dot/(n1*n2))

def reconstruct_vec(vtarg,vbasis):
    vrecon=np.zeros(np.shape(vtarg))
    for i in range(len(vbasis[0,:])):
        c = np.dot(vtarg,vbasis[:,i])
        
        vrecon+=c*vbasis[:,i]
        
    return vrecon

def reconstruct_basis(vecs,basisvecs):
    reconvecs=np.zeros(np.shape(vecs))
    for i in range (len(vecs[0,:])):
        reconvecs[:,i]=reconstruct_vec(vecs[:,i], basisvecs)
    
    return reconvecs


def get_projection_coefs(vtarg,vbasis,metric='cos'):
    coefs=np.zeros(len(vbasis[0,:]))
    
    for i in range(len(vbasis[0,:])):
        if metric=='cos':
            coefs[i]=cos_norm(vtarg,vbasis[:,i])
        if metric=='dot':
            coefs[i]=np.dot(vtarg,vbasis[:,i])
    
    return coefs

def recon_npercent_best(vtarg,vbasis,percent):
    coscoefs=get_projection_coefs(vtarg, vbasis)
    dotcoefs=get_projection_coefs(vtarg, vbasis,metric='dot')
    
    numcoefs=int(len(coscoefs)*(percent/100))
    
    sortinds=np.argsort(coscoefs)[::-1]
    sorted_coscoefs=coscoefs[sortinds]
    sorted_dotcoefs=dotcoefs[sortinds]
    sorted_vbasis=vbasis[:,sortinds]
    
    resultvec=np.zeros(len(vtarg))
    
    for i in range(numcoefs):
        resultvec+=sorted_dotcoefs[i]*sorted_vbasis[:,i]
    
    return resultvec

def get_recon_as_function_of_number_of_vecs_used(vtarg,vbasis):
    percents=np.arange(1,100)
    reconquals=[]
    for p in percents:
        recvec=recon_npercent_best(vtarg,vbasis,p)
        q = cos_norm(recvec,vtarg)
        reconquals.append(q)
    return reconquals

def get_number_of_vecs_needed(vtarg,vbasis,tol=.01):
    
    reconquals=get_recon_as_function_of_number_of_vecs_used(vtarg, vbasis)
    
    finalval=reconquals[-1]
    for i in range(len(reconquals)-1):
        if (finalval-reconquals[i])<=tol:
            return i

    
def get_av_num_vecs_needed(vecstarg,vecsbasis):
    num=0
    for i in range (len(vecstarg[0,:])):
        n=get_number_of_vecs_needed(vecstarg[:,i], vecsbasis)
        print(n)
        num+=n
    
    return num/len(vecstarg[0,:])


def subspace_distance_eff(v1,v2):
    start=time.time()
    v1=v1.astype('float32')
    v2=v2.astype('float32')
    p1=np.matmul(v1,v1.T)
    print('projmat1 done')
    p2=np.matmul(v2,v2.T)
    print('projmat2 done')
    dif=p1-p2
    print('D computed')
    u,singular_values,right_evecs=scipy.sparse.linalg.svds(dif,k=1,which='LM',return_singular_vectors=True)
    print('done')
    end=time.time()
    print(f'{end-start} seconds taken for dimension {np.shape(v1)}')
    return singular_values[0]