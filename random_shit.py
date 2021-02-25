#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 21:00:19 2021

@author: bwinston
"""
import numpy as np
from scipy.stats import entropy
from math import log, e
import pandas as pd

def entropy1(labels, base=None):
  value,counts = np.unique(labels, return_counts=True)
  return entropy(counts, base=base)

def entropy2(labels, base=None):
  """ Computes entropy of label distribution. """

  n_labels = len(labels)

  if n_labels <= 1:
    return 0

  value,counts = np.unique(labels, return_counts=True)
  probs = counts / n_labels
  n_classes = np.count_nonzero(probs)

  if n_classes <= 1:
    return 0

  ent = 0.

  # Compute entropy
  base = e if base is None else base
  for i in probs:
    ent -= i * log(i, base)

  return ent

def entropy3(labels, base=None):
  vc = pd.Series(labels).value_counts(normalize=True, sort=False)
  base = e if base is None else base
  return -(vc * np.log(vc)/np.log(base)).sum()

def entropy4(labels, base=None):
  value,counts = np.unique(labels, return_counts=True)
  norm_counts = counts / counts.sum()
  base = e if base is None else base
  return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()

labels = [1,3,5,2,3,5,3,2,1,3,4,5]

print(entropy1(labels))
print(entropy2(labels))
print(entropy3(labels))
print(entropy4(labels))

entropy1(jeff, base=2)

import matplotlib.pyplot as plt

fig,ax = plt.subplots(2,2,subplot_kw={'projection': '3d'})
plotting.plot_surf_stat_map([sc,si],ivp['105923']['test']['unmasked_vecs'][:,4],view='dorsal',cmap='RdBu',colorbar=True,output_file=None,title='4',vmax=.005,axes=ax[0][0],figure=fig)
plotting.plot_surf_stat_map([sc,si],ivp['105923']['test']['unmasked_vecs'][:,5],view='dorsal',cmap='RdBu',colorbar=True,output_file=None,title='5',vmax=.005,axes=ax[1][0],figure=fig)
plotting.plot_surf_stat_map([sc,si],ivp['105923']['test']['unmasked_vecs'][:,6],view='dorsal',cmap='RdBu',colorbar=True,output_file=None,title='6',vmax=.005,axes=ax[0][1],figure=fig)
plotting.plot_surf_stat_map([sc,si],ivp['105923']['test']['unmasked_vecs'][:,7],view='dorsal',cmap='RdBu',colorbar=True,output_file=None,title='7',vmax=.005,axes=ax[1][1],figure=fig)




