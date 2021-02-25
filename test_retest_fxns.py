#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:38:29 2021

@author: bwinsto2
"""
import numpy as np
import os
from glob import glob
import matrix_methods as mm
import input_output as inout
import decomp as dcp
from scipy import sparse
import numpy as np
import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked, pairwise
import time
from scipy.sparse import csgraph
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.stats import entropy
from scipy.spatial import distance
import statistics as stats  
import pandas as pd   
from nilearn import plotting
    
def get_key(my_dict, val):
    for key, value in my_dict.items():
         if val == value:
             return key

#looking at test retest reliability STILL ORDERING EFFECTS IN THIS
def test_retest_rel_cos(chap_dir, n_evecs):
    n_evecs = n_evecs-1
    global cd 
    cd = {} #for chap_data
    all_bcorrs, bcorr_plot = [], []
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        cd[sub] = {}
        cd[sub]['bcorrs'] = []
        for ses in ['test','retest']:
           cd[sub][ses] = {}
           cd[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy') 
           cd[sub][ses]['vecs'] = np.delete(cd[sub][ses]['vecs'], 0, axis=1)
        global hp
        hp, hp['holy'] = {}, {} #dict for within these fxns
        hp['corr_orig'] = np.empty((n_evecs,n_evecs))
        for evec_test in range(0,n_evecs): 
            for evec_retest in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                hp['corr_orig'][evec_retest,evec_test] = chap_cosine_sim(cd[sub]['test']['vecs'][:,evec_test], cd[sub]['retest']['vecs'][:,evec_retest]) #comparing column with column (ev with ev)
        hp['corr_all'] = hp['corr_orig'].swapaxes(0,1) #prepare to turn into dicts
        hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])} #turn into dicts
        find_bcorrs(sub, hp, 0, n_evecs) #run find bcorrs function, get best pairs
        for ev in range(n_evecs):
            cd[sub]['bcorrs'].append(cd[sub]['pairs'][ev]['bcorr']) #all ideal corrs (w/ no repeats)
    for sub in subs:        
        all_bcorrs.append(cd[sub]['bcorrs']) #list of lists of bcorrs
    cd['bcorr_avg'] = np.average(np.array(all_bcorrs), axis=0) #average of each spot in bcorrs
    #for avg in range(1,n_evecs):
        #bcorr_plot.append(stats.mean(cd['bcorr_avg'][1:avg])) 
    #plt.plot(bcorr_plot)

def test_retest_rel_r(chap_dir, n_evecs):
    n_evecs = n_evecs-1
    global cd 
    cd = {} #for chap_data
    all_bcorrs = []
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        cd[sub] = {}
        cd[sub]['bcorrs'] = []
        for ses in ['test','retest']:
           cd[sub][ses] = {}
           cd[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy') 
           cd[sub][ses]['vecs'] = np.delete(cd[sub][ses]['vecs'], 0, axis=1)
        global hp
        hp, hp['holy'] = {}, {} #dict for within these fxns
        hp['corr_orig'] = np.empty((n_evecs,n_evecs))
        for evec_test in range(0,n_evecs): 
            for evec_retest in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                hp['corr_orig'][evec_retest,evec_test] = np.arctanh(abs(pearsonr(cd[sub]['test']['vecs'][:,evec_test], cd[sub]['retest']['vecs'][:,evec_retest])[0])) #comparing column with column (ev with ev)
        hp['corr_all'] = hp['corr_orig'].swapaxes(0,1) #prepare to turn into dicts
        hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])} #turn into dicts
        find_bcorrs(sub, hp, 0, n_evecs) #run find bcorrs function, get best pairs
        for ev in range(n_evecs):
            cd[sub]['bcorrs'].append(cd[sub]['pairs'][ev]['bcorr']) #all ideal corrs (w/ no repeats)
    for sub in subs:        
        all_bcorrs.append(cd[sub]['bcorrs']) #list of lists of bcorrs
    cd['bcorr_avg'] = np.average(np.array(all_bcorrs), axis=0) #average of each spot in bcorrs
    #for avg in range(1,n_evecs):
        #bcorr_plot.append(stats.mean(cd['bcorr_avg'][1:avg])) 
    #plt.plot(bcorr_plot)

def find_bcorrs(sub, hp, run, n_evecs): 
    while len(hp['corr_all']) > 0:
        hp[run] = {}
        hp[run]['maxes'] = {}
        for ev in hp['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[run]['maxes'][ev] = {} #test session ev's dict
            hp[run]['maxes'][ev]['bcorr'] = max(hp['corr_all'][ev].values()) #best correlation w/ any retest ev
            hp[run]['maxes'][ev]['ret_ind'] = get_key(hp['corr_all'][ev], max(hp['corr_all'][ev].values())) #which retest ev was above?
        hp[run]['ret_used'], hp[run]['win_ind'], hp[run]['good_boys'] = [], [], []
        for ev in hp['corr_all']:
            hp[run]['ret_used'].append(hp[run]['maxes'][ev]['ret_ind']) #tally retest indices used above
        ret_dups = set([x for x in hp[run]['ret_used'] if hp[run]['ret_used'].count(x) > 1]) #retest evs that were used multiple times
        if ret_dups == 0:
            break
        hp[run]['ret_used'] = set(hp[run]['ret_used']) #take out duplicates from ret_used
        for ret_dup in ret_dups: #for each retest ev in demand by multiple
            competition = [] #competition is within the retest ev
            for ev in hp[run]['maxes']: 
                if ret_dup == hp[run]['maxes'][ev]['ret_ind']: #if a test evec is involved in a dispute
                    competition.append(hp[run]['maxes'][ev]['bcorr']) #add its correlation to retest ev's competition 
            competition = max(competition) 
            for ev in hp[run]['maxes']:
                if hp[run]['maxes'][ev]['bcorr'] == competition: #if its correlation was a winner, leave as is
                    hp[run]['win_ind'].append(ev) #add test ev that won to win_ind
        for ev in hp['corr_all']:
            if hp[run]['maxes'][ev]['ret_ind'] not in ret_dups: #if test ev not involved in any disputes, add to good_boys
               hp[run]['good_boys'].append(ev) 
        for ev in list(hp['corr_all']):
            if ev in (set(hp[run]['good_boys']) | set(hp[run]['win_ind'])):
               hp['holy'][ev] = hp[run]['maxes'][ev] #add good boys and winners to holy
               del hp['corr_all'][ev] #delete them from corr_all
        for ev in hp['corr_all']:
            for r_ev in hp[run]['ret_used']:
                del hp['corr_all'][ev][r_ev] #idk what this is doing
        run = run + 1
        find_bcorrs(sub, hp, run, len(hp['corr_all'])) #rerun function on smaller corr_all (leftovers)
    cd[sub]['pairs'] = hp['holy']

'''
average connectome/across subjects shit:
'''
#generate average-connectome harmonics
def avg_harms(chap_dir):
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    global avg
    avg = {}
    avg['connectomes_test'], avg['connectomes_retest'], avg['connectomes_all']= 0,0,0
    for sub in subs:
       avg['connectomes_test'] = avg['connectomes_test'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-test/connectome.npz') #sum all test connectomes
       avg['connectomes_retest'] = avg['connectomes_retest'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-retest/connectome.npz')
       avg['connectomes_all'] = avg['connectomes_all'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-test/connectome.npz') + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-retest/connectome.npz')
    avg['connectomes_test'] = avg['connectomes_test'] / len(subs) #take average
    avg['connectomes_retest'] = avg['connectomes_retest'] / len(subs)
    avg['connectomes_all'] = avg['connectomes_all'] / (len(subs)*2)
    inout.if_not_exist_make(f'{chap_dir}/sub-test_avg')
    inout.if_not_exist_make(f'{chap_dir}/sub-retest_avg')
    inout.if_not_exist_make(f'{chap_dir}/sub-total_avg')
    avg_test_vals,avg_test_vecs = dcp.lapDecomp(avg['connectomes_test'], 200) 
    avg_retest_vals,avg_retest_vecs = dcp.lapDecomp(avg['connectomes_retest'], 200) 
    avg_total_vals,avg_total_vecs = dcp.lapDecomp(avg['connectomes_all'], 200) 
    np.save(f'{chap_dir}/sub-test_avg/vecs', avg_test_vecs)
    np.save(f'{chap_dir}/sub-retest_avg/vecs', avg_retest_vecs)
    np.save(f'{chap_dir}/sub-total_avg/vecs', avg_total_vecs)

#average connectome harmonics as components compare to individuals
def ind_vs_avg(chap_dir, n_evecs):
    global iva
    iva = {} 
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        subs.remove(sub)
    for sub in subs:
        iva[sub] = {}    
        for ses in ['test','retest']:
           iva[sub][ses],iva[f'{ses}_avg'] = {}, {}
           iva[sub][ses]['bcorrs'] = []
           iva[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
           iva[f'{ses}_avg']['vecs'] = np.load(f'{chap_dir}/sub-{ses}_avg/vecs.npy')
        global hp
        hp = {}
        for ses in ['test','retest']:
            hp[ses], hp[ses]['holy'] = {}, {}  
            hp[ses]['corr_orig'] = np.empty((n_evecs,n_evecs)) #init comparison matrix 
            for evec_avg in range(0,n_evecs): 
                for evec_ses in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                    hp[ses]['corr_orig'][evec_ses,evec_avg] = abs(pearsonr(iva[f'{ses}_avg']['vecs'][:,evec_avg], iva[sub][ses]['vecs'][:,evec_ses])[0])
            hp[ses]['corr_all'] = hp[ses]['corr_orig'].swapaxes(0,1)
            hp[ses]['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp[ses]['corr_all'])}
            find_bcorrs_iva(sub, ses, hp, 0, hp[ses]['corr_all'], n_evecs)
    for sub in subs:
        for ses in ['test','retest']:
            for ev in range(n_evecs):
                iva[sub][ses]['bcorrs'].append(iva[sub][ses]['pairs'][ev]['bcorr'])
    iva['all_test_bcorrs'],iva['all_retest_bcorrs'] = [],[]
    for sub in subs:
        for ses in ['test','retest']:
            iva[f'all_{ses}_bcorrs'].append(iva[sub][ses]['bcorrs'])
    iva['bcorr_test_avg'] = np.average(np.array(iva['all_test_bcorrs']), axis = 0)
    iva['bcorr_retest_avg'] = np.average(np.array(iva['all_retest_bcorrs']), axis = 0)
        
#finds best correlations btwn avg and ind
def find_bcorrs_iva(sub, ses, hp, run, corr_mat, n_evecs): 
    while len(hp[ses]['corr_all']) > 0:
        hp[ses][run] = {}
        hp[ses][run]['maxes'] = {}
        for ev in hp[ses]['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[ses][run]['maxes'][ev] = {}
            hp[ses][run]['maxes'][ev]['bcorr'] = max(hp[ses]['corr_all'][ev].values())
            hp[ses][run]['maxes'][ev][f'{ses}_ind'] = get_key(hp[ses]['corr_all'][ev], max(hp[ses]['corr_all'][ev].values()))
        hp[ses][run][f'{ses}_used'], hp[ses][run]['win_ind'], hp[ses][run]['good_boys'] = [], [], []
        for ev in hp[ses]['corr_all']:
            hp[ses][run][f'{ses}_used'].append(hp[ses][run]['maxes'][ev][f'{ses}_ind'])
        ses_dups = set([x for x in hp[ses][run][f'{ses}_used'] if hp[ses][run][f'{ses}_used'].count(x) > 1])
        if ses_dups == 0:
            break
        hp[ses][run][f'{ses}_used'] = set(hp[ses][run][f'{ses}_used'])
        for ses_dup in ses_dups:
            competition = []
            for ev in hp[ses][run]['maxes']:
                if ses_dup == hp[ses][run]['maxes'][ev][f'{ses}_ind']:
                    competition.append(hp[ses][run]['maxes'][ev]['bcorr'])
            competition = max(competition)
            for ev in hp[ses][run]['maxes']:
                if hp[ses][run]['maxes'][ev]['bcorr'] == competition:
                    hp[ses][run]['win_ind'].append(ev)
        for ev in hp[ses]['corr_all']:
            if hp[ses][run]['maxes'][ev][f'{ses}_ind'] not in ses_dups:
               hp[ses][run]['good_boys'].append(ev) 
        for ev in list(hp[ses]['corr_all']):
            if ev in (set(hp[ses][run]['good_boys']) | set(hp[ses][run]['win_ind'])):
               hp[ses]['holy'][ev] = hp[ses][run]['maxes'][ev]
               del hp[ses]['corr_all'][ev]
        for ev in hp[ses]['corr_all']:
            for s_ev in hp[ses][run][f'{ses}_used']:
                del hp[ses]['corr_all'][ev][s_ev]
        run = run + 1
        find_bcorrs_iva(sub, ses, hp, run, hp[ses]['corr_all'], len(hp[ses]['corr_all']))
    iva[sub][ses]['pairs'] = hp[ses]['holy']

'''
Take everyone from test session, generate 100 pca harmonics. for each subject, find the best test session pairs with 
those 100 pca harmonics. 

Then, find the average correlation btwn each pca harmonic and its best pair across all subjects--this shows 
the most reliable components across subjects (maybe). Pick 20.

Then, for each subject, find test-retest reliability of those particular harmonics
(the ones that paired w/ the 20 best pcas) with their best pair from the retest session.
'''
def ind_vs_pca(chap_dir, n_evecs, n_comp, mask, sc, si):
    n_evecs = n_evecs-1
    global ivp
    ivp, ivp['test'], ivp['retest'] = {}, {}, {}
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    ivp['evlist'] = []
    for sub in subs:
        ivp[sub] = {}    
        for ses in ['test','retest']:
           ivp[sub][ses],ivp[f'{ses}_avg'] = {}, {}
           ivp[sub][ses]['bcorrs'],ivp[sub][ses]['inds'] = [],[]
           ivp[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
           ivp[sub][ses]['vecs'] = np.delete(ivp[sub][ses]['vecs'], 0, axis=1)
           ivp[sub][ses]['unmasked_vecs'] = np.empty([64984,n_evecs])
           mask = np.load('/Users/bwinston/Documents/connectome_harmonics/hcp_mask.npy')
           for ev in range(n_evecs):
               ivp[sub][ses]['unmasked_vecs'][:,ev]=ut.unmask_medial_wall(ivp[sub][ses]['vecs'][:,ev],mask)
           ivp['evlist'].append(ivp[sub][ses]['vecs'][:,0:n_evecs])
    ivp['pca_harms'] = dcp.get_group_pca_comp_brian(ivp['evlist'], n_comp, n_evecs)[0]
    ivp['variance'] = dcp.get_group_pca_comp_brian(ivp['evlist'], n_comp, n_evecs)[1]
    ivp['unmasked_pca_harms'] = np.empty([64984,n_comp])
    for ev in range(n_comp):
        ivp['unmasked_pca_harms'][:,ev] = ut.unmask_medial_wall(ivp['pca_harms'][:,ev],mask)
    for ev in range(n_evecs):
        ivp['unmasked_pca_harms']
    for sub in subs:
        global hp
        hp = {}
        for ses in ['test','retest']:
            hp[ses], hp[ses]['holy'] = {}, {}  
            hp[ses]['corr_orig'] = np.empty((n_evecs,n_comp)) #init comparison matrix 
            for evec_pca in range(0,n_comp): 
                for evec_ses in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                    hp[ses]['corr_orig'][evec_ses,evec_pca] = abs(pearsonr(ivp['pca_harms'][:,evec_pca], ivp[sub][ses]['vecs'][:,evec_ses])[0])
            hp[ses]['corr_all'] = hp[ses]['corr_orig'].swapaxes(0,1)
            hp[ses]['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp[ses]['corr_all'])}           
            find_bcorrs_ivp(sub, ses, hp, 0, n_evecs)  
    for sub in subs:
        for ses in ['test','retest']:
            for ev in range(n_comp):
                ivp[sub][ses]['bcorrs'].append(ivp[sub][ses]['pairs'][ev]['bcorr'])
                ivp[sub][ses]['inds'].append(ivp[sub][ses]['pairs'][ev][f'{ses}_ind'])
    ivp['all_test_bcorrs'],ivp['all_retest_bcorrs'] = [],[]
    for sub in subs:
        for ses in ['test','retest']:
            ivp[f'all_{ses}_bcorrs'].append(ivp[sub][ses]['bcorrs']) #list of lists
    ivp['bcorr_test_avg'] = np.average(np.array(ivp['all_test_bcorrs']), axis = 0)
    ivp['bcorr_retest_avg'] = np.average(np.array(ivp['all_retest_bcorrs']), axis = 0)
    ivp['pca_to_match_avg'] = (ivp['bcorr_retest_avg']+ivp['bcorr_retest_avg'])/2
    ivp['PCs'] = {}
    for pc in range(n_comp):
        ivp['PCs'][pc] = {}
        display = plotting.plot_surf_stat_map([sc,si],ivp['unmasked_pca_harms'][:,pc],view='dorsal',cmap='RdBu',output_file=None,colorbar=True,title=f'PC{pc}',vmax=.005)
        plotting.show()
        plt.close(display)
        fig,ax = plt.subplots(len(subs),2,subplot_kw={'projection': '3d'})
        fig.suptitle(f'PC {pc}')
        row = 0
        for sub in subs:
            ivp['PCs'][pc][sub] = {}
            for ses in ['test','retest']:
                ivp['PCs'][pc][sub][ses] = ivp[sub][ses]['pairs'][pc][f'{ses}_ind'] #which vec to plot
                ind = ivp['PCs'][pc][sub][ses]
                col=0 if ses=='test' else 1    
                if pearsonr(ivp['unmasked_pca_harms'][:,pc],ivp[sub][ses]['unmasked_vecs'][:,ind])[0] < 0:
                    vec = np.negative(ivp[sub][ses]['unmasked_vecs'][:,ind])
                else:
                    vec = ivp[sub][ses]['unmasked_vecs'][:,ind]
                display = plotting.plot_surf([sc,si],vec,view='dorsal',cmap='RdBu',output_file=None,colorbar=False,vmax=.005,vmin=-.005,axes=ax[row][col],figure=fig)
                ax[row][col].set_title(f'{sub} {ses} harmonic {ind}', fontsize=10)
            row = row+1
        plotting.show()
        plt.close(fig)
    
    bcorr_t_sort = -np.sort(-ivp['bcorr_test_avg'])[::1]
    bcorr_r_sort = -np.sort(-ivp['bcorr_retest_avg'])[::1]
    bcorr_t_sort = bcorr_t_sort[:20]
    bcorr_r_sort = bcorr_r_sort[:20]
    ivp['pca_t_harms'], ivp['pca_r_harms'] = [], []
    for c in bcorr_t_sort:
        ivp['pca_t_harms'].append(np.where(ivp['bcorr_test_avg'] == c)[0][0]) #best pca harmonics
    for c in bcorr_r_sort:
        ivp['pca_r_harms'].append(np.where(ivp['bcorr_retest_avg'] == c)[0][0]) 
    for sub in subs:
        ivp[sub]['enchilada'] = {} #enchilada is test-retest but with the test ones we care about bc they were similar to PCAs
        ivp[sub]['enzymatique'] = []
        ivp[sub]['ti_2_check'] = []
        for pca_harm in ivp['pca_t_harms']:
            ivp[sub]['ti_2_check'].append(ivp[sub]['test']['pairs'][pca_harm]['test_ind'])
    test_retest_rel(chap_dir, n_evecs)
    all_enzymatiques = []
    for sub in subs:
        for ti in ivp[sub]['ti_2_check']:
            ivp[sub]['enchilada'][ti] = cd[sub]['pairs'][ti]
            ivp[sub]['enzymatique'].append(ivp[sub]['enchilada'][ti]['bcorr'])
        all_enzymatiques.append(ivp[sub]['enzymatique'])
    ivp['test_retest_reliabilty_avg'] = np.average(np.array(all_enzymatiques), axis = 0)

    
def find_bcorrs_ivp(sub, ses, hp, run, n_evecs): 
    while len(hp[ses]['corr_all']) > 0:
        hp[ses][run] = {}
        hp[ses][run]['maxes'] = {}
        for ev in hp[ses]['corr_all']: #for each pca evec, which {ses} evec is the best match?
            hp[ses][run]['maxes'][ev] = {}
            hp[ses][run]['maxes'][ev]['bcorr'] = max(hp[ses]['corr_all'][ev].values())
            hp[ses][run]['maxes'][ev][f'{ses}_ind'] = get_key(hp[ses]['corr_all'][ev], max(hp[ses]['corr_all'][ev].values()))
        hp[ses][run][f'{ses}_used'], hp[ses][run]['win_ind'], hp[ses][run]['good_boys'] = [], [], []
        for ev in hp[ses]['corr_all']:
            hp[ses][run][f'{ses}_used'].append(hp[ses][run]['maxes'][ev][f'{ses}_ind'])
        ses_dups = set([x for x in hp[ses][run][f'{ses}_used'] if hp[ses][run][f'{ses}_used'].count(x) > 1])
        if ses_dups == 0:
            break
        hp[ses][run][f'{ses}_used'] = set(hp[ses][run][f'{ses}_used'])
        for ses_dup in ses_dups:
            competition = []
            for ev in hp[ses][run]['maxes']:
                if ses_dup == hp[ses][run]['maxes'][ev][f'{ses}_ind']:
                    competition.append(hp[ses][run]['maxes'][ev]['bcorr'])
            competition = max(competition)
            for ev in hp[ses][run]['maxes']:
                if hp[ses][run]['maxes'][ev]['bcorr'] == competition:
                    hp[ses][run]['win_ind'].append(ev)
        for ev in hp[ses]['corr_all']:
            if hp[ses][run]['maxes'][ev][f'{ses}_ind'] not in ses_dups:
               hp[ses][run]['good_boys'].append(ev) 
        for ev in list(hp[ses]['corr_all']):
            if ev in (set(hp[ses][run]['good_boys']) | set(hp[ses][run]['win_ind'])):
               hp[ses]['holy'][ev] = hp[ses][run]['maxes'][ev]
               del hp[ses]['corr_all'][ev]
        for ev in hp[ses]['corr_all']:
            for s_ev in hp[ses][run][f'{ses}_used']:
                del hp[ses]['corr_all'][ev][s_ev]
        run = run + 1
        find_bcorrs_ivp(sub, ses, hp, run, len(hp[ses]['corr_all']))
    ivp[sub][ses]['pairs'] = hp[ses]['holy']

'''2v shit'''

def test_retest_rel_2v(vec_1, vec_2, n_evecs):
    n_evecs = n_evecs-1
    global cd 
    cd = {} #for chap_data
    cd['bcorrs'] = []
    if type(vec_1) != np.ndarray:
        cd['vec_1'] = np.load(vec_1)
    else:
        cd['vec_1'] = vec_1
    if type(vec_2) != np.ndarray:   
        cd['vec_2'] = np.load(vec_2)
    else:
        cd['vec_2'] = vec_2
    for i in ['1', '2']:
        cd[f'vec_{i}'] = np.delete(cd[f'vec_{i}'], 0, axis=1)
    global hp
    hp, hp['holy'] = {}, {} #dict for within these fxns
    hp['corr_orig'] = np.empty((n_evecs,n_evecs))
    for evec_1 in range(0,n_evecs): 
        for evec_2 in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
            hp['corr_orig'][evec_2,evec_1] = np.arctanh(abs(pearsonr(cd['vec_1'][:,evec_1], cd['vec_2'][:,evec_2])[0])) #comparing column with column (ev with ev)
    hp['corr_all'] = hp['corr_orig'].swapaxes(0,1) #prepare to turn into dicts
    hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])} #turn into dicts
    find_bcorrs_2v(hp, 0, n_evecs) #run find bcorrs function, get best pairs
    for ev in range(n_evecs):
        cd['bcorrs'].append(cd['pairs'][ev]['bcorr']) #all ideal corrs (w/ no repeats)
    return cd['bcorrs']

def find_bcorrs_2v(hp, run, n_evecs): 
    while len(hp['corr_all']) > 0:
        hp[run] = {}
        hp[run]['maxes'] = {}
        for ev in hp['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[run]['maxes'][ev] = {} #test session ev's dict
            hp[run]['maxes'][ev]['bcorr'] = max(hp['corr_all'][ev].values()) #best correlation w/ any retest ev
            hp[run]['maxes'][ev]['ret_ind'] = get_key(hp['corr_all'][ev], max(hp['corr_all'][ev].values())) #which retest ev was above?
        hp[run]['ret_used'], hp[run]['win_ind'], hp[run]['good_boys'] = [], [], []
        for ev in hp['corr_all']:
            hp[run]['ret_used'].append(hp[run]['maxes'][ev]['ret_ind']) #tally retest indices used above
        ret_dups = set([x for x in hp[run]['ret_used'] if hp[run]['ret_used'].count(x) > 1]) #retest evs that were used multiple times
        if ret_dups == 0:
            break
        hp[run]['ret_used'] = set(hp[run]['ret_used']) #take out duplicates from ret_used
        for ret_dup in ret_dups: #for each retest ev in demand by multiple
            competition = [] #competition is within the retest ev
            for ev in hp[run]['maxes']: 
                if ret_dup == hp[run]['maxes'][ev]['ret_ind']: #if a test evec is involved in a dispute
                    competition.append(hp[run]['maxes'][ev]['bcorr']) #add its correlation to retest ev's competition 
            competition = max(competition) 
            for ev in hp[run]['maxes']:
                if hp[run]['maxes'][ev]['bcorr'] == competition: #if its correlation was a winner, leave as is
                    hp[run]['win_ind'].append(ev) #add test ev that won to win_ind
        for ev in hp['corr_all']:
            if hp[run]['maxes'][ev]['ret_ind'] not in ret_dups: #if test ev not involved in any disputes, add to good_boys
               hp[run]['good_boys'].append(ev) 
        for ev in list(hp['corr_all']):
            if ev in (set(hp[run]['good_boys']) | set(hp[run]['win_ind'])):
               hp['holy'][ev] = hp[run]['maxes'][ev] #add good boys and winners to holy
               del hp['corr_all'][ev] #delete them from corr_all
        for ev in hp['corr_all']:
            for r_ev in hp[run]['ret_used']:
                del hp['corr_all'][ev][r_ev] #idk what this is doing
        run = run + 1
        find_bcorrs_2v(hp, run, len(hp['corr_all'])) #rerun function on smaller corr_all (leftovers)
    cd['pairs'] = hp['holy']     
    
##METRIC 1
#chap_out = '/Users/bwinston/Downloads/chap_out_test'        

def struc_metric_1(chap_dir, n_evecs):
    global sm 
    sm = {} #for chap_data
    sm['within_subj_all'], sm['across_subj_all'] = [],[]
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        sm[sub] = {}
        for ses in ['test','retest']:
           sm[sub][ses] = {}
           sm[sub][ses]['vecs'] = f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy'
        sm[sub][sub] = stats.mean(test_retest_rel_2v(sm[sub]['test']['vecs'], sm[sub]['retest']['vecs'], n_evecs))
        sm['within_subj_all'].append(sm[sub][sub])
        sm[sub]['c_sub_all'] = [] #where empty averages will go
    for sub in subs:
        for c_sub in subs:
            if c_sub != sub:
                sm[sub][c_sub] = {}
                for ses in ['test','retest']:
                    sm[sub][c_sub][ses] = stats.mean(test_retest_rel_2v(sm[sub][ses]['vecs'], sm[c_sub][ses]['vecs'], n_evecs))
                sm[sub][c_sub]['avg'] = (sm[sub][c_sub]['test'] + sm[sub][c_sub]['retest']) / 2
                sm[sub]['c_sub_all'].append(sm[sub][c_sub]['avg'])
        sm['across_subj_all'].append(stats.mean(sm[sub]['c_sub_all']))
    sm['within_subj_avg'] = stats.mean(sm['within_subj_all']) 
    sm['across_subj_avg'] = stats.mean(sm['across_subj_all'])
    return sm['within_subj_avg'], sm['across_subj_avg']
                   
#let's check where it stops separating!        
def struc_metric_1_sep(chap_dir, n_evecs):
    print('ehllo')
    global sm 
    sm = {} #for chap_data
    sm['within_subj_all'], sm['across_subj_all'] = [],[]
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        sm[sub] = {}
        for ses in ['test','retest']:
           sm[sub][ses] = {}
           sm[sub][ses]['vecs'] = f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy'
        sm[sub][sub] = test_retest_rel_2v(sm[sub]['test']['vecs'], sm[sub]['retest']['vecs'], n_evecs)
        sm['within_subj_all'].append(sm[sub][sub])
        sm[sub]['c_sub_all'] = [] #where empty averages will go
    sm['within_subj_avg'] = np.average(np.array(sm['within_subj_all']), axis=0)
    for sub in subs:
        sm[sub]['across_all'], sm[sub]['across_avg'] = [], []
        for c_sub in subs:
            if c_sub != sub:
                sm[sub][c_sub] = {}
                for ses in ['test','retest']:
                    sm[sub]['across_all'].append(test_retest_rel_2v(sm[sub][ses]['vecs'], sm[c_sub][ses]['vecs'], n_evecs))
        sm[sub]['across_avg'] = np.average(np.array(sm[sub]['across_all']), axis=0)
        sm['across_subj_all'].append(sm[sub]['across_avg'])
        print(f'finished {sub}')
    sm['across_subj_avg'] = np.average(np.array(sm['across_subj_all']), axis=0)
    plt.plot(sm['within_subj_avg'])
    plt.plot(sm['across_subj_avg'])
    return sm['within_subj_avg'], sm['across_subj_avg']        
         
##Metric 2 (sparsity shit)
def struc_metric_2(chap_dir):
    global sp 
    print('hey')
    sp = {} #for chap_data
    sp['within_subj_all'], sp['across_subj_all'] = [],[]
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    c_subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for c_sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{c_sub}'):
            c_subs.remove(c_sub)
    for sub in c_subs:
        sp[sub] = {}
        for ses in ['test','retest']:
           sp[sub][ses] = {}
           sp[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
        sp[sub][sub] = dcp.get_av_num_vecs_needed(sp[sub]['test']['vecs'], sp[sub]['retest']['vecs'])
        sp['within_subj_all'].append(sp[sub][sub])
        sp[sub]['c_sub_all'] = [] #where empty averages will go
    print('within stuff done')

def struc_metric_2_across(chap_dir):   
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        for c_sub in subs:
            if c_sub != sub:
                sp[sub][c_sub] = {}
                for ses in ['test','retest']:
                    sp[sub][c_sub][ses] = dcp.get_av_num_vecs_needed(sp[sub][ses]['vecs'], sp[c_sub][ses]['vecs'])
                sp[sub][c_sub]['avg'] = (sp[sub][c_sub]['test'] + sp[sub][c_sub]['retest']) / 2
                sp[sub]['c_sub_all'].append(sp[sub][c_sub]['avg'])
        sp['across_subj_all'].append(stats.mean(sp[sub]['c_sub_all']))
        print(f'sub {sub} done')
    sp['across_subj_avg'] = stats.mean(sp['across_subj_all'])
    return sp['across_subj_avg'], sp['across_subj_all']

#metric_2_forreal_across = struc_metric_2_across('/data/hcp_test_retest_pp/derivatives/chap')      
 
'''
spectra/functional shit
'''

#frequency

def zero_crossings(ts):
    return ((ts[:-1] * ts[1:]) < 0).sum()        
        
def ev_freq(chap_dir, n_evecs):
    global s
    s = {}
    s['freq_avg'] = []
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in subs:
        s[sub] = {}
        for ses in ['test','retest']:
            s[sub][ses] = {}
            for scan in ['REST1', 'REST2']:
                s[sub][ses][scan] = {}
                s[sub][ses][scan]['power'], s[sub][ses][scan]['energy'] = {}, {}
                s[sub][ses][scan]['recon'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/reconspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_dynamic_reconstruction_spectrum.npy')
                s[sub][ses][scan]['recon'] = np.delete(s[sub][ses][scan]['recon'], 0, axis=0)
                for spec in ['power', 'energy']:
                    for t in ['mean', 'dynamic']:
                        s[sub][ses][scan][spec][t] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/{spec}spectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_{t}_{spec}_spectrum.npy')
                        s[sub][ses][scan][spec][t] = np.delete(s[sub][ses][scan][spec][t], 0, axis=0)
                s[sub][ses][scan]['power']['normalized'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/powerspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_normalized_power_spectrum.npy')
                s[sub][ses][scan]['power']['normalized'] = np.delete(s[sub][ses][scan]['power']['normalized'], 0, axis=0)
    for sub in subs:
        s[sub]['freq'] = []
        for ev in range(n_evecs-1):
            s[sub][ev] = []
            for ses in ['test','retest']:
                for scan in ['REST1', 'REST2']:
                    s[sub][ev].append(zero_crossings(s[sub][ses][scan]['recon'][ev]))
            s[sub]['freq'].append(stats.mean(s[sub][ev]))
        s['freq_avg'].append(s[sub]['freq'])
    s['freq_avg'] = np.average(np.array(s['freq_avg']), axis=0)
    plt.plot(s['freq_avg'])
    
def freq_comp(chap_dir, n_evecs):
    global s
    s = {}
    s['freq_avg'] = []
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in subs:
        s[sub] = {}
        for ses in ['test','retest']:
            s[sub][ses] = {}
            for scan in ['REST1', 'REST2']:
                s[sub][ses][scan] = {}
                s[sub][ses][scan]['power'], s[sub][ses][scan]['energy'] = {}, {}
                s[sub][ses][scan]['recon'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/reconspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_dynamic_reconstruction_spectrum.npy')
                s[sub][ses][scan]['recon'] = np.delete(s[sub][ses][scan]['recon'], 0, axis=0)
                for spec in ['power', 'energy']:
                    for t in ['mean', 'dynamic']:
                        s[sub][ses][scan][spec][t] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/{spec}spectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_{t}_{spec}_spectrum.npy')
                        s[sub][ses][scan][spec][t] = np.delete(s[sub][ses][scan][spec][t], 0, axis=0)
                s[sub][ses][scan]['power']['normalized'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/powerspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_normalized_power_spectrum.npy')
                s[sub][ses][scan]['power']['normalized'] = np.delete(s[sub][ses][scan]['power']['normalized'], 0, axis=0)  
    s['within_subj_all'] = []
    s['across_subj_all'] = []
    for sub in subs:
        for ses in ['test','retest']:
            s[sub][ses]['freq'] = []
            for ev in range(n_evecs-1):
                s[sub][ses][ev] = []
                for scan in ['REST1', 'REST2']:
                    s[sub][ses][ev].append(zero_crossings(s[sub][ses][scan]['recon'][ev]))
                s[sub][ses][ev] = stats.mean(s[sub][ses][ev])
                s[sub][ses]['freq'].append(s[sub][ses][ev])
        s[sub][sub] = np.arctanh(pearsonr(s[sub]['test']['freq'], s[sub]['retest']['freq'])[0])
        s['within_subj_all'].append(s[sub][sub])
    for sub in subs:
        s[sub]['c_sub_all'] = []
        for c_sub in subs:
            if c_sub != sub:
                s[sub][c_sub] = {}
                for ses in ['test','retest']:
                    s[sub][c_sub][ses] = pearsonr(s[sub][ses]['freq'], s[c_sub][ses]['freq'])[0]
                s[sub][c_sub]['avg'] = (s[sub][c_sub]['test'] + s[sub][c_sub]['retest']) / 2
                s[sub]['c_sub_all'].append(s[sub][c_sub]['avg'])
        s['across_subj_all'].append(stats.mean(s[sub]['c_sub_all']))
    s['within_subj_avg'] = stats.mean(s['within_subj_all']) 
    s['across_subj_avg'] = stats.mean(s['across_subj_all'])
    return s['within_subj_avg'], s['across_subj_avg']
           
'''shan entropy'''
def shan_entropy_kite(dist):
    #https://www.kite.com/python/answers/how-to-calculate-shannon-entropy-in-python
    pd_series = pd.Series(dist)
    counts = pd_series.value_counts()
    shan_entropy = entropy(counts)   
    return shan_entropy

def shan_entropy(dist):
    bins = list(np.linspace(-65, 65, 40))
    p = np.digitize(dist, bins=bins)
    return shan_entropy_kite(p)

def ev_shan_ent(chap_dir, n_evecs):
    global se
    se = {}
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in subs:
        se[sub] = {}
        for ses in ['test','retest']:
            se[sub][ses] = {}
            for scan in ['REST1', 'REST2']:
                se[sub][ses][scan] = {}
                se[sub][ses][scan]['power'], se[sub][ses][scan]['energy'] = {}, {}
                se[sub][ses][scan]['recon'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/reconspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_dynamic_reconstruction_spectrum.npy')
                se[sub][ses][scan]['recon'] = np.delete(se[sub][ses][scan]['recon'], 0, axis=0)
                for spec in ['power', 'energy']:
                    for t in ['mean', 'dynamic']:
                        se[sub][ses][scan][spec][t] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/{spec}spectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_{t}_{spec}_spectrum.npy')
                        se[sub][ses][scan][spec][t] = np.delete(se[sub][ses][scan][spec][t], 0, axis=0)
                se[sub][ses][scan]['power']['normalized'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/func/{scan}/powerspectra/sub-{sub}_ses-{ses}_task-{scan.lower()}_normalized_power_spectrum.npy')
                se[sub][ses][scan]['power']['normalized'] = np.delete(se[sub][ses][scan]['power']['normalized'], 0, axis=0)
    se['across_sub_avg'] = []
    se['within_sub_avg'] = []
    for sub in subs:
        for ses in ['test','retest']:
            se[sub][ses]['ent'] = []
            for scan in ['REST1', 'REST2']:
                se[sub][ses][scan]['ent'] = []
                for ev in range(n_evecs-1):
                    se[sub][ses][scan]['ent'].append(shan_entropy(se[sub][ses][scan]['recon'][ev])) 
                se[sub][ses]['ent'].append(se[sub][ses][scan]['ent'])
            se[sub][ses]['ent'] = np.average(np.array(se[sub][ses]['ent']), axis=0)
        se[sub]['within_sub'] = np.arctanh(pearsonr(se[sub]['test']['ent'], se[sub]['retest']['ent'])[0])
        se['within_sub_avg'].append(se[sub]['within_sub'])
    for sub in subs:
        se[sub]['ent_corrs'] = []
        for c_sub in subs:
            if c_sub != sub:
                se[sub][c_sub] = {}
                for ses in ['test','retest']:
                    se[sub][c_sub][ses] = np.arctanh(pearsonr(se[sub][ses]['ent'], se[c_sub][ses]['ent'])[0])
                se[sub]['ent_corrs'].append((se[sub][c_sub]['test'] + se[sub][c_sub]['retest'])/2)
        se[sub]['across_subs'] = stats.mean(se[sub]['ent_corrs'])
        se['across_sub_avg'].append(se[sub]['across_subs'])
    se['across_sub_avg'] = stats.mean(se['across_sub_avg'])
    se['within_sub_avg'] = stats.mean(se['within_sub_avg'])
    return se['within_sub_avg'], se['across_sub_avg']

'''RSN STUFF'''
#takes output of ivp, net_verts
def vecs_vs_rsn(ivp,net_verts, chap_dir):
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for network in list(set(net_verts)):
        ivp[network], ivp[network]['pcs'] = {}, {}
        for pc in range(len(ivp['pca_harms'][0])):
            ivp[network]['pcs'][pc] = []
            for ses in ['test','retest']:
                for sub in subs:
                   ind = ivp[sub][ses]['pairs'][pc][f'{ses}_ind']
                   ivp[network]['pcs'][pc].append(abs(pearsonr(ivp[sub][ses]['vecs'][:,ind],net_verts[network]['verts'])[0]))
            ivp[network]['pcs'][pc] = stats.mean(ivp[network]['pcs'][pc])
        ivp[network][f'pcs_{network}'] = []
        for pc in range(len(ivp['pca_harms'][0])):
            ivp[network][f'pcs_{network}'].append(ivp[network]['pcs'][pc])
        plt.plot(ivp[network][f'pcs_{network}'])
        plt.xlabel(network)
        plt.show()

#vecs_vs_rsn(ivp,net_verts,'/Users/bwinston/Downloads/chap_out_test')


'''                    
                        
                    sm[sub][c_sub][ses] = stats.mean(test_retest_rel_2v(sm[sub][ses]['vecs'], sm[c_sub][ses]['vecs'], n_evecs))
                sm[sub][c_sub]['avg'] = (sm[sub][c_sub]['test'] + sm[sub][c_sub]['retest']) / 2
                sm[sub]['c_sub_all'].append(sm[sub][c_sub]['avg'])
        sm['across_subj_all'].append(stats.mean(sm[sub]['c_sub_all']))
    sm['within_subj_avg'] = stats.mean(sm['within_subj_all']) 
    sm['across_subj_avg'] = stats.mean(sm['across_subj_all'])
    return sm['within_subj_avg'], sm['across_subj_avg']

'''