
#!/bin/env python

import numpy as np 
import pandas as pd
import pystan
import scipy.io as scipy
import _pickle as cPickle
#from scipy.io 
import stan_utility
import pandas as pd
import os
import pickle


def inv_logit(arr):
    '''Elementwise inverse logit (logistic) function.'''
    return 1 / (1 + np.exp(-arr))

def phi_approx(arr):
    '''Elementwise fast approximation of the cumulative unit normal. 
    For details, see Bowling et al. (2009). "A logistic approximation 
    to the cumulative normal distribution."'''
    return inv_logit(0.07056 * arr ** 3 + 1.5976 * arr)

def run_model():
	# Create standata (load from data generated by dataForQ_bySession_all.m)

	data = loadmat('/Volumes/witten/Paper_code/human_bandit/06-Feb-2022/data_toStan.mat')

	NS = data['NS']
	r  = data['r']
	c  = data['c']
	NT = data['NT']
	NT_all = data['NT_all']

	r=r.astype(int)
	c=c.astype(int)

	standata = {'NS':NS[0][0],'NT':NT[0][0],'NT_all':NT_all[:,0], 
           'r':r, 'c':c}


	# Compile the model
	#sm = CmdStanModel(stan_file="/qlearning_stay_side.stan")
	sm = pystan.StanModel(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/qmodel/','qlearning_stay_side.stan'))



	# Fit model
	fit = sm.sampling(data=standata, iter=1000, warmup=250, chains = 4, control=dict(adapt_delta=0.99, max_treedepth=15))


	print(fit)

	# Save data
	summary = fit.summary()
	summary = pd.DataFrame(summary['summary'], columns=summary['summary_colnames'], index=summary['summary_rownames'])
	summary.to_csv(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/','summary_int_pc.csv'))

	extract = fit.extract()
	for k, v in extract.items(): extract[k] = v
	with open(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/' 'StanFit_int_pc.pickle'), 'wb') as fn: cPickle.dump(extract, fn)

	samplerParams = fit.get_sampler_params()
	samplerParams = pd.DataFrame(samplerParams)
	summary.to_csv(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'samplerParams.csv'))

	with open(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'fit.pkl'), 'wb') as f:
		pickle.dump({'model' : sm, 'fit' : fit}, f, protocol=-1)

	with open(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/' 'fit2.pkl'), 'wb') as f:
		pickle.dump([sm, fit], f, protocol=-1)

	#    Extract data and save to mat
	alphas = phi_approx(extract['alphas'])
	betas = extract['betas']
	sides = extract['sides']
	stays = extract['stays']

	scipy.savemat(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'alphas.mat'),{'vect':alphas})
	scipy.savemat(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'betas.mat'),{'vect':betas})
	scipy.savemat(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'sides.mat'),{'vect':sides})
	scipy.savemat(os.path.join('/jukebox/witten/Paper_Code/Cox_2021/processed_data/dataForQ/human_bandit/all/', 'stays.mat'),{'vect':stays})


	print(pystan.check_hmc_diagnostics(fit))
