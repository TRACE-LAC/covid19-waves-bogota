# -*- coding: utf-8 -*-
"""
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions
@author: davidsantiagoquevedo
@author: ntorresd
"""
import yaml
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

strat='wave'

ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'epidemiological_distributions')
SCRIPTS_PATH = config['PATHS']['SCRIPTS_PATH'].format(dir = 'epidemiological_distributions')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'epidemiological_distributions')

import sys
sys.path.append(UTILS_PATH)
import utilities_epi_dist as ut
SEED = config['MODELS']['SEED']
ITER = config['MODELS']['ITER']
CHAINS = config['MODELS']['CHAINS']
MIN_VAL = config['MODELS']['MIN_VAL']
MAX_VAL = config['MODELS']['MAX_VAL']

# Model variables
strat_name = 'wave'
n_strats = 4
params = ['beta']

# 1. Prepare the data
all_dfs, columns = ut.prepare_confirmed_cases_data()

##############################################################################
# print n samples and range of data
for df in all_dfs:
    col = str(df.columns[4])
    print(col, len(df[col].index), df[col].min(), '-', df[col].max())
    
##############################################################################
# load the samples (models fits)

dist_posteriors_gamma = {}
dist_posteriors_wei = {}
dist_posteriors_lognorm = {}
dist_posteriors_exp = {}
dist_posteriors_gln = {}

for i in range(len(columns)):
    col = columns[i]
    dist_posteriors_gamma.update({col: pd.read_csv(OUT_PATH + col +'-samples-gamma.csv')})
    dist_posteriors_lognorm.update({col: pd.read_csv(OUT_PATH + col +'-samples-log_normal.csv')})
    dist_posteriors_wei.update({col: pd.read_csv(OUT_PATH + col +'-samples-weibull.csv')})
    dist_posteriors_exp.update({col: pd.read_csv(OUT_PATH + col +'-samples-exponential.csv')})  
    try:
        dist_posteriors_gln.update({col: pd.read_csv(OUT_PATH + col + '-samples-gln.csv')})
    except:
        print(col, 'did not find samples for gln')

##############################################################################
# Laplace approx & Bayes Factors to find the best fit

def fit_posterior(distr, col):
    if distr == 'gamma':
        df = dist_posteriors_gamma[col].iloc[:,:-2]
        samples = df.values
    elif distr == 'lognormal':
        df = dist_posteriors_lognorm[col].iloc[:,:-2]
        samples = df.values
    elif distr == 'weibull':
        df = dist_posteriors_wei[col].iloc[:,:-2]
        samples = df.values
    elif distr == 'exponential':
        df = dist_posteriors_exp[col].iloc[:,:-1]
        samples = df.values
    elif distr == 'gln':
        df = dist_posteriors_gln[col].iloc[:,:-3]
        samples = df.values
    else:
        print('Invalid distribution')
        return {}
    cov = np.cov(np.array(samples), rowvar = False)
    theta_mean = samples.mean(axis = 0)
    return {'mu': theta_mean, 'cov': cov, 'distr': distr}

columns_BF = columns

posterior_gamma = {}
posterior_wei = {}
posterior_lognorm = {}
posterior_exp = {}
posterior_gln = {}

for col in columns_BF:
    posterior_gamma.update({col: fit_posterior('gamma', col)})
    posterior_lognorm.update({col: fit_posterior('lognormal', col)})
    posterior_wei.update({col: fit_posterior('weibull', col)})
    posterior_exp.update({col: fit_posterior('exponential', col)})
    posterior_gln.update({col: fit_posterior('gln', col)})
    

def Logf(posterior, col):
    col_id = columns.index(col)
    data = all_dfs[col_id]
    distr = posterior[col]['distr']
    theta = posterior[col]['mu']
    x = data[col].values + 0.5 
    strat_vals = data[strat+'_id'].values # int index for the state
    k = len(data[strat+'_id'].unique())
    
    log_prior = 0
    log_likelihood = 0
    
    if distr == 'gamma':
        sigma_alpha = theta[-2]
        sigma_beta = theta[-1]
        alpha_dist = dist_posteriors_gamma[col].alpha.values.mean()
        beta_dist =  dist_posteriors_gamma[col].beta.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_alpha)
        log_prior += stats.norm(0,1).logpdf(sigma_beta)
        for i in range(k):
            log_prior += stats.norm(alpha_dist,sigma_alpha).logpdf(theta[i]) # alphas
            log_prior += stats.norm(beta_dist,sigma_beta).logpdf(theta[i+k]) # betas
        for i in range(len(x)):
            alpha = theta[strat_vals[i]-1]
            beta = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.gamma(a=alpha, scale = 1.0/beta).logpdf(x[i])

    elif distr == 'lognormal':
        sigma_mu = theta[-2]
        sigma_sigma = theta[-1]
        mu_dist = dist_posteriors_lognorm[col].mu.values.mean()
        sigma_dist = dist_posteriors_lognorm[col].sigma.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_mu)
        log_prior += stats.norm(0,1).logpdf(sigma_sigma)
        for i in range(k):
            log_prior += stats.norm(mu_dist,sigma_mu).logpdf(theta[i]) # mu
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
        for i in range(len(x)):
            mu = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.lognorm(s=sigma, scale = np.exp(mu)).logpdf(x[i])  
            
    elif distr == 'weibull':
        sigma_alpha = theta[-2]
        sigma_sigma = theta[-1]
        alpha_dist = dist_posteriors_wei[col].alpha.values.mean()
        sigma_dist = dist_posteriors_wei[col].sigma.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_alpha)
        log_prior += stats.norm(0,1).logpdf(sigma_sigma)
        for i in range(k):
            log_prior += stats.norm(alpha_dist,sigma_alpha).logpdf(theta[i]) # alphas
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
        for i in range(len(x)):
            alpha = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            log_likelihood += stats.weibull_min(c=alpha, scale = sigma).logpdf(x[i])
    
    elif distr == 'exponential':
        sigma_beta = theta[-2]
        beta_dist = dist_posteriors_exp[col].beta.values.mean()
        log_prior += stats.norm(0,1).logpdf(sigma_beta)
        for i in range(k):
            log_prior += stats.norm(beta_dist,sigma_beta).logpdf(theta[i]) # betas
        for i in range(len(x)):
            beta = theta[strat_vals[i]-1]
            log_likelihood += stats.expon(scale = 1/beta).logpdf(x[i])
            
    elif distr == 'gln':
        sigma_mu = theta[-3]
        sigma_sigma = theta[-2]
        sigma_g = theta[-1]
        mu_dist = dist_posteriors_gln[col].mu.values.mean()
        sigma_dist = dist_posteriors_gln[col].sigma.values.mean()
        g_dist = dist_posteriors_gln[col].g.values.mean()
        log_prior += stats.norm(2,0.5).logpdf(sigma_mu)
        log_prior += stats.norm(0.5,0.5).logpdf(sigma_sigma)
        log_prior += stats.norm(1.5,0.5).logpdf(sigma_g)

        for i in range(k):
            log_prior += stats.norm(mu_dist,sigma_mu).logpdf(theta[i]) # mu
            log_prior += stats.norm(sigma_dist,sigma_sigma).logpdf(theta[i+k]) # sigmas
            log_prior += stats.norm(g_dist,sigma_g).logpdf(theta[i+2*k]) # g
        for i in range(len(x)):
            mu = theta[strat_vals[i]-1]
            sigma = theta[strat_vals[i]-1 + k]
            g = theta[strat_vals[i]-1 + 2*k]
            log_likelihood += ut.gln_lpdf(x[i], mu, sigma, g)
    else:
        print('Invalid distribution')
        return np.nan

    return log_likelihood + log_prior    

# careful, this is slow
for col in columns_BF:
    print(col)
    posterior_gamma[col].update({'Logf': Logf(posterior_gamma, col)})
    posterior_gamma[col].update({'LogLaplace': ut.LogLaplaceCovariance(posterior_gamma, col)})
    
    posterior_lognorm[col].update({'Logf': Logf(posterior_lognorm, col)})
    posterior_lognorm[col].update({'LogLaplace': ut.LogLaplaceCovariance(posterior_lognorm, col)})
    
    
    posterior_wei[col].update({'Logf': Logf(posterior_wei, col)})
    posterior_wei[col].update({'LogLaplace': ut.LogLaplaceCovariance(posterior_wei, col)})
    
    posterior_exp[col].update({'Logf': Logf(posterior_exp, col)})
    posterior_exp[col].update({'LogLaplace': ut.LogLaplaceCovariance(posterior_exp, col)})
    
    posterior_gln[col].update({'Logf': Logf(posterior_gln, col)})
    posterior_gln[col].update({'LogLaplace': ut.LogLaplaceCovariance(posterior_gln, col)})
    
    
columns_models = ['Gamma', 'Lognormal', 'Weibull', 'Exponential', 'Gen Lognormal']
bayesfactorsBR = pd.DataFrame(columns = columns_models, index = columns_models)
BayesFactorsdist = {}
evidence_dict = {'Gamma': posterior_gamma, 'Lognormal': posterior_lognorm,
                 'Exponential': posterior_exp, 'Weibull': posterior_wei,
                 'Gen Lognormal': posterior_gln
                } 
                  
all_BF = {}
#row - col (negative: in favor of col; positive: in favor of row)
for col in columns_BF:
    for c1 in columns_models:
        for c2 in columns_models:
            bayesfactorsBR.loc[c1,c2] = 2*(evidence_dict[c1][col]['LogLaplace'] - evidence_dict[c2][col]['LogLaplace'])
    all_BF.update({col: bayesfactorsBR})
    print(col)
    print(bayesfactorsBR, '\n')
    df_bf=pd.DataFrame(bayesfactorsBR)
    df_bf.to_csv(f'{OUT_PATH}/bf_{col}.csv')