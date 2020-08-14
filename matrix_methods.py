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




def construct_structural_connectivity_matrix(SC,EC,tol=3,NNnum=45):
    '''
    SC- array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    EC- array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]
    tol- search radius of nearest neighbor search for matching endpoints to surface vertices
    NNnum- number of nearest neighboring surface vertices to assign to each endpoint
    '''
    ind,dist=ut.neighbors(SC,EC,1)
    bad=[]
    c=np.arange(len(dist))
    even=c[::2]
    odd=c[1::2]
    for i in range (int(len(dist)/2)):
        if (dist[even[i]]>=tol or dist[odd[i]]>=tol):
            bad.append(even[i])
            bad.append(odd[i])
    newEC=np.delete(EC,bad,axis=0)
    s2eInd, s2eDist=ut.neighbors(SC,newEC,1)
    Rind,Rdist=ut.neighbors(newEC,SC,NNnum)
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
    print(M.nnz)

    x=np.arange(len(SC))
    M[x,x]=0
    print(M.nnz)
    M=M.tocsr()
    return M.tocsr()


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
