"""
@author: patricktaylor
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.ticker import MaxNLocator
from statistics import stdev
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from powerlaw import plot_pdf, Fit, pdf
import powerlaw
import plfit
import math
from sklearn.metrics import mean_squared_error
import pandas as pd


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

def plot_spectrum(spectrum, spectrum_type):
    #TODO: replace with args.evecs
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Wavenumber, Ψ')
    plt.ylabel(spectrum_type)
    plt.bar(range(20), spectrum)
    plt.show()  

def criticality(mean_spectrum, dynamic_spectrum, spectrum_type, comparison_distributions=['exponential', 'truncated_power_law', 'lognormal']):
    #alpha is fitted parameter for the powerlaw distribution, sigma the standard error
    #R is the loglikelihood ratio between the distributions,  
    max_spectrum = [max(row) for row in list(dynamic_spectrum)]
    std_spectrum = [stdev(row) for row in list(dynamic_spectrum)]
    criticality_df_row = pd.DataFrame()
    for dist in comparison_distributions:
        fit_mean = powerlaw.Fit(mean_spectrum)
        fit_max = powerlaw.Fit(max_spectrum)
        fit_std = powerlaw.Fit(std_spectrum)
        R_mean, p_mean = fit_mean.distribution_compare('power_law', dist)
        R_max, p_max = fit_max.distribution_compare('power_law', dist)
        R_std, p_std = fit_std.distribution_compare('power_law', dist)
        criticality_df_row[f'mps_powerlaw_vs_{dist}'] = [{'R':R_mean, 'p':p_mean, 'alpha':fit_mean.power_law.alpha, 'sigma':fit_mean.power_law.sigma}]
        criticality_df_row[f'maxps_powerlaw_vs_{dist}'] = [{'R':R_max, 'p':p_max, 'alpha':fit_max.power_law.alpha, 'sigma':fit_max.power_law.sigma}]
        criticality_df_row[f'stdps_powerlaw_vs_{dist}'] = [{'R':R_std, 'p':p_std, 'alpha':fit_std.power_law.alpha, 'sigma':fit_std.power_law.sigma}]
        print(f'Mean {spectrum_type} spectrum: ')
        print(f'Powerlaw vs. {dist} R: '+str(R_mean)+" p: "+str(p_mean))
        print(f'Max {spectrum_type} spectrum: ')
        print(f'Powerlaw vs. {dist} R: '+str(R_max)+" p: "+str(p_max))
        print(f'STD {spectrum_type} spectrum: ')
        print(f'Powerlaw vs. {dist} R: '+str(R_std)+" p: "+str(p_std))
    return criticality_df_row

def plot_rmse_criticality(mean_spectrum, dynamic_spectrum, spectrum_type):
        #TODO: replace with args.evecs
    wavenumbers = list(np.arange(20+1))[1:]
    mean_spectrum = list(mean_spectrum)
    max_spectrum = [max(row) for row in list(dynamic_spectrum)]
    std_spectrum = [stdev(row) for row in list(dynamic_spectrum)]
    log_wavenumbers = list(map(lambda x: math.log(x, 10), wavenumbers))
    log_mean_spectrum = list(map(lambda x: math.log(x,10), mean_spectrum))
    log_std_spectrum = list(map(lambda x: math.log(x,10), max_spectrum))
    log_max_spectrum = list(map(lambda x: math.log(x,10), std_spectrum))
    fig, axs=plt.subplots(3,1)
    axs[0].scatter(log_wavenumbers, log_mean_spectrum)
    
    axs[0].set_ylabel("Mean "+spectrum_type)
    coefficient_best_fit_mean = np.polyfit(log_wavenumbers, log_mean_spectrum, 1)
    pred_best_fit_mean =  list(map(lambda x: x*coefficient_best_fit_mean[0]+coefficient_best_fit_mean[1], log_wavenumbers))
    #[[i, pred_best_fit[i]] for i in range(20)],
    line_mean = mlines.Line2D(log_wavenumbers, pred_best_fit_mean, color='red')
    axs[0].add_line(line_mean)
    print("RSME Mean powerlaw fit: "+str(math.sqrt(mean_squared_error(mean_spectrum, pred_best_fit_mean))))
    axs[1].scatter(log_wavenumbers, log_std_spectrum)
    axs[1].set_ylabel("Std "+spectrum_type)
    coefficient_best_fit_std = np.polyfit(log_wavenumbers, log_std_spectrum, 1)
    pred_best_fit_std =  list(map(lambda x: x*coefficient_best_fit_std[0]+coefficient_best_fit_std[1], log_wavenumbers))
    #[[i, pred_best_fit[i]] for i in range(20)],
    line_std = mlines.Line2D(log_wavenumbers, pred_best_fit_std, color='red')
    axs[1].add_line(line_std)
    print("RSME STD powerlaw fit: "+str(math.sqrt(mean_squared_error(std_spectrum, pred_best_fit_mean))))
    axs[2].scatter(log_wavenumbers, log_max_spectrum)
    axs[2].set_xlabel('Log Wavenumber, Ψ')
    axs[2].set_ylabel("Max "+spectrum_type)
    coefficient_best_fit_max = np.polyfit(log_wavenumbers, log_max_spectrum, 1)
    pred_best_fit_max =  list(map(lambda x: x*coefficient_best_fit_max[0]+coefficient_best_fit_max[1], log_wavenumbers))
    #[[i, pred_best_fit[i]] for i in range(20)],
    line = mlines.Line2D(log_wavenumbers, pred_best_fit_max, color='red')
    axs[2].add_line(line)
    print("RSME Max powerlaw fit: "+str(math.sqrt(mean_squared_error(max_spectrum, pred_best_fit_max))))   