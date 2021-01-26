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
import matplotlib.pylab as plt
from scipy.stats.stats import pearsonr
from scipy.spatial import distance
import statistics as stats     
    
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
                hp['corr_orig'][evec_retest,evec_test] = abs(pearsonr(cd[sub]['test']['vecs'][:,evec_test], cd[sub]['retest']['vecs'][:,evec_retest])[0]) #comparing column with column (ev with ev)
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
def ind_vs_pca(chap_dir, n_evecs, n_pca):
    n_evecs = n_evecs-1
    global ivp
    ivp, ivp['test'], ivp['retest'] = {}, {}, {}
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    ivp['test']['evlist'], ivp['retest']['evlist'] = [], []
    for sub in subs:
        ivp[sub] = {}    
        for ses in ['test','retest']:
           ivp[sub][ses],ivp[f'{ses}_avg'] = {}, {}
           ivp[sub][ses]['bcorrs'] = []
           ivp[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
           ivp[sub][ses]['vecs'] = np.delete(ivp[sub][ses]['vecs'], 0, axis=1)
           ivp[ses]['evlist'].append(ivp[sub][ses]['vecs'][:,0:n_evecs])
    for ses in ['test','retest']:
        ivp[ses]['pca_harms'] = dcp.get_group_pca_comp_b(ivp[ses]['evlist'], 100)
    for sub in subs:
        global hp
        hp = {}
        for ses in ['test','retest']:
            hp[ses], hp[ses]['holy'] = {}, {}  
            hp[ses]['corr_orig'] = np.empty((n_evecs,n_pca)) #init comparison matrix 
            for evec_pca in range(0,n_pca): 
                for evec_ses in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                    hp[ses]['corr_orig'][evec_ses,evec_pca] = abs(pearsonr(ivp[ses]['pca_harms'][:,evec_pca], ivp[sub][ses]['vecs'][:,evec_ses])[0])
            hp[ses]['corr_all'] = hp[ses]['corr_orig'].swapaxes(0,1)
            hp[ses]['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp[ses]['corr_all'])}           
            find_bcorrs_ivp(sub, ses, hp, 0, hp[ses]['corr_all'], n_evecs)  
    for sub in subs:
        for ses in ['test','retest']:
            for ev in range(n_pca):
                ivp[sub][ses]['bcorrs'].append(ivp[sub][ses]['pairs'][ev]['bcorr'])
    ivp['all_test_bcorrs'],ivp['all_retest_bcorrs'] = [],[]
    for sub in subs:
        for ses in ['test','retest']:
            ivp[f'all_{ses}_bcorrs'].append(ivp[sub][ses]['bcorrs']) #list of lists
    ivp['bcorr_test_avg'] = np.average(np.array(ivp['all_test_bcorrs']), axis = 0)
    ivp['bcorr_retest_avg'] = np.average(np.array(ivp['all_retest_bcorrs']), axis = 0)
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

    
def find_bcorrs_ivp(sub, ses, hp, run, corr_mat, n_evecs): 
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
        find_bcorrs_ivp(sub, ses, hp, run, hp[ses]['corr_all'], len(hp[ses]['corr_all']))
    ivp[sub][ses]['pairs'] = hp[ses]['holy']

'''2v shit'''

def test_retest_rel_2v(vec_1, vec_2, n_evecs):
    n_evecs = n_evecs-1
    global cd 
    cd = {} #for chap_data
    cd['bcorrs'] = []
    cd['vec_1'] = np.load(vec_1)
    cd['vec_2'] = np.load(vec_2)
    for i in ['1', '2']:   
        cd[f'vec_{i}'] = np.delete(cd[f'vec_{i}'], 0, axis=1)
    global hp
    hp, hp['holy'] = {}, {} #dict for within these fxns
    hp['corr_orig'] = np.empty((n_evecs,n_evecs))
    for evec_1 in range(0,n_evecs): 
        for evec_2 in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
            hp['corr_orig'][evec_2,evec_1] = abs(pearsonr(cd[f'vec_1'][:,evec_1], cd[f'vec_2'][:,evec_2])[0]) #comparing column with column (ev with ev)
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
    sm['across_subj_avg'] = np.average(np.array(sm['across_subj_all']), axis=0)
    plt.plot(sm['within_subj_avg'])
    plt.plot(sm['across_subj_avg'])
    return sm['within_subj_avg'], sm['across_subj_avg']        
         
##Metric 2 (sparsity shit)
def struc_metric_2(chap_dir):
    global sp 
    sp = {} #for chap_data
    sp['within_subj_all'], sp['across_subj_all'] = [],[]
    subject_dirs = glob(os.path.join(chap_dir, "sub-*")) #get subs
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        if os.path.exists(f'{chap_dir}/sub-{sub}'):
            subs.remove(sub)
    for sub in subs:
        sp[sub] = {}
        for ses in ['test','retest']:
           sp[sub][ses] = {}
           sp[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
        sp[sub][sub] = dcp.get_av_num_vecs_needed(sp[sub]['test']['vecs'], sp[sub]['retest']['vecs'])
        sp['within_subj_all'].append(sp[sub][sub])
        sp[sub]['c_sub_all'] = [] #where empty averages will go
    for sub in subs:
        for c_sub in subs:
            if c_sub != sub:
                sp[sub][c_sub] = {}
                for ses in ['test','retest']:
                    sp[sub][c_sub][ses] = dcp.get_av_num_vecs(sp[sub][ses]['vecs'], sp[c_sub][ses]['vecs'])
                sp[sub][c_sub]['avg'] = (sp[sub][c_sub]['test'] + sp[sub][c_sub]['retest']) / 2
                sp[sub]['c_sub_all'].append(sp[sub][c_sub]['avg'])
        sp['across_subj_all'].append(stats.mean(sp[sub]['c_sub_all']))
    sp['within_subj_avg'] = stats.mean(sp['within_subj_all']) 
    sp['across_subj_avg'] = stats.mean(sp['across_subj_all'])
    return sp['within_subj_avg'], sp['across_subj_avg']
        
        
        
        
        

