#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:27:26 2020

@author: patricktaylor
"""
from scipy import sparse
import numpy as np
import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked
import time
from scipy.sparse import csgraph
import gdist
import sklearn

def construct_struc_conn(sc,ec,tol=2):
    ind,dist=ut.neighbors(sc,ec,1)
    indstart=ind[::2][:,0]
    indend=ind[1::2][:,0]
    diststart=dist[::2][:,0]
    distend=dist[1::2][:,0]
    
    tolmask= (diststart>tol) | (distend>tol)
    
    M=sparse.lil_matrix((len(sc),len(sc)))
    
    for i in range (len(indend)):
        
        if not tolmask[i]:
            
            M[indstart[i],indend[i]]+=1
    
    return M.tocsr()

def construct_incidence_matrices(sc,ec,tol=2):
    start=time.time()
    ind,dist=ut.neighbors(sc,ec,1)
    indstart=ind[::2][:,0]
    indend=ind[1::2][:,0]
    diststart=dist[::2][:,0]
    distend=dist[1::2][:,0]
    
    tolmask= (diststart>tol) | (distend>tol)
    
    startinc=sparse.lil_matrix( (len(sc), (len(tolmask) -tolmask.sum()) ), dtype=np.float32)
    endinc=sparse.lil_matrix( (len(sc), (len(tolmask) -tolmask.sum()) ), dtype=np.float32)
    
    good_inds=[i for i in range (len(tolmask)) if not tolmask[i]]
    x=np.arange(len(good_inds))
    startinc[indstart[good_inds],x]=1
    endinc[indend[good_inds],x]=1
    
    #for j,i in enumerate(good_inds):
        
    #    startinc[indstart[i],j]=1
    #    endinc[indend[i],j]=1
    end=time.time()
    print(f'{end-start} seconds taken for incidence matrix construction')
    return startinc.tocsr(), endinc.tocsr()
    
def construct_smoothing_matrix(sc,si,mask,sigma=2,epsilon=0.05):
    start=time.time()
    h=int(len(sc)/2)
    hsi=int(len(si)/2)
    lg=construct_smoothing_matrix_one_hem(sc[:h], si[:hsi],sigma=sigma,epsilon=epsilon)
    rg=construct_smoothing_matrix_one_hem(sc[h:], si[hsi:]-h ,sigma=sigma,epsilon=epsilon)
    
    sm=sparse.vstack((sparse.hstack((lg,sparse.csr_matrix((lg.shape[0], rg.shape[1]), dtype=lg.dtype))).tocsr(),sparse.hstack((sparse.csr_matrix((rg.shape[0], lg.shape[1]), dtype=lg.dtype),rg)).tocsr()))
    sm=uts.mask_connectivity_matrix(sm,mask)
    sm=sklearn.preprocessing.normalize(sm, norm='l1')
    end=time.time()
    print (end-start, 'seconds taken for smoothing matrix construction')
    return sm
    
def construct_smoothing_matrix_one_hem(sc,si,sigma=2,epsilon=0.05):
    maxd=sigma * (-2 * np.log(epsilon)) ** (1 / 2)
    dists=gdist.local_gdist_matrix(sc.astype(np.float64), si.astype(np.int32),maxd )
    dists[dists > maxd] = 0
    dists = dists.minimum(dists.T)
    dists.eliminate_zeros()
    dists = dists.tolil()
    dists.setdiag(0)
    dists = dists.tocsr()
    
    g = -(dists.power(2) / (2 * (sigma ** 2)))
    np.exp(g.data, out=g.data)
    g += sparse.eye(g.shape[0], dtype=g.dtype).tocsr()

    return g
    
def smooth_incidence_matrices(start, end, coefs,binarize=False,return_unsmoothed=False):
    
    if binarize:
        M=start.dot(end.T)
        if return_unsmoothed:
            mat=M
        M[M>1]=1
        
        M=coefs.T.dot(M.dot(coefs))
    else:
        smooth_start = start.T.dot(coefs).T
        smooth_end = end.T.dot(coefs).T
        M = smooth_start.dot(smooth_end.T)
    if return_unsmoothed:
        return M+M.T,mat+mat.T
    else:
        return M+M.T

def construct_smoothed_connectivity_matrix(sc,si,ec,mask,tol=2,sigma=3,epsilon=0.2,binarize=False,return_unsmoothed=False):
    start=time.time()
    
    starti,endi=construct_incidence_matrices(uts.mask_medial_wall_vecs(sc,mask), ec, tol)
    print('incidence matrices computed')
    
    smoothing_coefs=construct_smoothing_matrix(sc,si,mask,sigma,epsilon)
        
    print('smoothing coefficients computed')
    if return_unsmoothed:
        smoothA,A=smooth_incidence_matrices(starti,endi,smoothing_coefs,binarize=binarize,return_unsmoothed=return_unsmoothed)
        end=time.time()
        print(f'{end-start} seconds taken')
        return smoothA,A
    else:
        A=smooth_incidence_matrices(starti,endi,smoothing_coefs,binarize=binarize)
    
        end=time.time()
        print(f'{end-start} seconds taken')
        return A


def construct_surface_matrix(SC,SI):
    """
     SC- surface coordinates (used for determining size of surface matrix only, can be empty array of length=len(SC))
     SI - array of vertex connections. Each row of SI contains indices of vertices in SC that form a triangle in the mesh
     """
    M=sparse.lil_matrix((len(SC),len(SC))) #initialize
    #set each M_ij, where i,j are vertices connected in the mesh, equal to 1.
    M[SI[:,0],SI[:,1]]=1
    M[SI[:,1],SI[:,0]]=1
    M[SI[:,1],SI[:,2]]=1
    M[SI[:,2],SI[:,1]]=1
    M[SI[:,2],SI[:,0]]=1
    M[SI[:,0],SI[:,2]]=1
    return M.tocsr()


def construct_inter_hemi_matrix(SC,tol=4):
    """
    creates interhemisphere connection matrix IHC for a given surface mesh with coordinates SC by connecting all vertices on the midline fissure
    to their nearest neighbor on the opposite hemisphere. IHC has dimension (len(SC),len(SC)).
    """
    half=int(len(SC)/2)
    li,ld=ut.neighbors(SC[:half],SC[half:],1)
    ri,rd=ut.neighbors(SC[half:],SC[:half],1)
    IHC=sparse.lil_matrix((half*2,half*2))
    #R=sparse.lil_matrix((half,half))
    for i in range (half):
        if ld[i]<tol:
            IHC[i+half,li[i]]=1
            IHC[li[i],i+half]=1
    for i in range (half):
        if rd[i]<tol:
            IHC[ri[i]+half,i]=1
            IHC[i,ri[i]+half]=1
    IHC=IHC.tocsr()
    return IHC



def construct_structural_connectivity_matrix(SC,EC,tol=3,NNnum=45,binarize=True):
    '''
    SC- array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    EC- array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]. also, EC[0] and EC[1] are endpoints of the same streamline etc.
    tol- maximum search radius of nearest neighbor search for matching endpoints to surface vertices (in mm)
    NNnum- number of nearest neighboring surface vertices to assign to each endpoint
    '''
    ind,dist=ut.neighbors(SC,EC,1) #computes 1 nearest neighbor of each ec in sc space and returns numpy arrays size (len(ec),1) of 1) the index NN surface vertex, and 2) distance from endpoint to vertex
    bad=[] 
    c=np.arange(len(dist)) 
    even=c[::2] #0,2,4,6, etc.
    odd=c[1::2] #1,3,5,7, etc.
    for i in range (int(len(dist)/2)): #from 0 to length of ec/2
        if (dist[even[i]]>=tol or dist[odd[i]]>=tol): #if distance of ec to sc at either endpoint is greater than tol, throw out streamline
            bad.append(even[i])
            bad.append(odd[i])
    newEC=np.delete(EC,bad,axis=0) #throw out bad fibers
    s2eInd, s2eDist=ut.neighbors(SC,newEC,1) #compute new ind and dist from "filtered" streamlines
    Rind,Rdist=ut.neighbors(newEC,SC,NNnum) #find NNum of endpoints for each surface vertex as a minimum
    OtherEndInd=np.zeros(np.shape(Rind))
    for i in range(len(Rind)):
        for j in range (NNnum):
            if Rind[i][j]%2==0 :
                OtherEndInd[i][j]=int(s2eInd[Rind[i][j]+1])
            else:
                OtherEndInd[i][j]=int(s2eInd[Rind[i][j]-1])

    M=sparse.lil_matrix((len(SC),len(SC)))
    x=np.arange(len(SC))
    for i in range (NNnum):
        AccSurfInd=np.column_stack((x,OtherEndInd[:,i]))
        U,C=np.unique(AccSurfInd,axis=0,return_counts=True)
        M[U[:,0],U[:,1]]+=C
        M[U[:,1],U[:,0]]+=C 

    x=np.arange(len(SC))
    M[x,x]=0 #remove self connections?
    print(M.nnz/2) #nonzero elements/2 (connections)
    M=M.tocsr()
    if binarize:
        M[M>=1]=1
    return M.tocsr(),Rdist

def diffusion_matrix(A,t=0):
    Lap,D_sqrt_vec=csgraph.laplacian(A,normed=True,return_diag=True)
    D_inv_sqrt_vec=1/D_sqrt_vec
    D_inv_sqrt=sparse.lil_matrix(np.shape(Lap))
    x=np.arange(np.shape(Lap)[0])
    D_inv_sqrt[x,x]=D_inv_sqrt_vec
    L_sqrt=(D_inv_sqrt.dot(A)).dot(D_inv_sqrt)
    Lap,D_alpha=csgraph.laplacian(L_sqrt,return_diag=True)
    D_alpha_mat=sparse.lil_matrix(np.shape(Lap))
    D_alpha_mat[x,x]=1/D_alpha
    P=D_alpha_mat.dot(L_sqrt)
    print('diffusion matrix computed')
    if t==0:
        return P
    else:
        return P**t

def cosine_similarity_matrix(M):
    gen = pairwise_distances_chunked(M,metric='cosine',n_jobs=-1)
    cosine_similarity=sparse.lil_matrix(np.shape(M))
    start=0
    for item in gen:
        end=start+len(item)
        g=1-item
        #g=np.where(g<0.1,0,g)
        print(np.shape(g))
        cosine_similarity[start:end,:]=g
        start=end
    #cosine_similarity=1-cosine_distance
    x=np.arange(np.shape(M)[0])
    cosine_similarity[x,x]=0
    return cosine_similarity
