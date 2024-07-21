# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""

import yaml
import pandas as pd
import numpy as np
from scipy import stats
import scipy as scipy
import matplotlib.pyplot as plt
import mpl_axes_aligner as mpla
from matplotlib.ticker import MaxNLocator
import seaborn as sns

config = yaml.load(open("config.yml", "r"))["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'epidemiological_distributions')

import sys
sys.path.append(UTILS_PATH)
import utilities_epi_dist as ut

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

##############################################################################
# 1. Prepare the data
all_dfs, columns = ut.prepare_confirmed_cases_data()

##############################################################################
# load the samples (models fits for every epidemiological distribution)
dist_posteriors = ut.load_samples()

# get best models
df_best_models = ut.best_model()
df_best_models = df_best_models.transpose()
df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")

########################################################################
# Plot CDFs
def plot_cdf(df, epi_dist, max_val, ax, cdf_null_hyp, n, n_subset = None, subset = 'wave'):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    if n_subset:
        data = df.loc[(df[subset] == n_subset) & (df[epi_dist]<=max_val)][epi_dist]
        ext = f'[{n_subset}]'
    else:
        data = df.loc[df[epi_dist]<=max_val][epi_dist]   
        ext = ''
    data = data + 0.5
    data_sorted = np.sort(np.array(data.drop_duplicates()))
    #get cdf from obs data
    count, bins_count = np.histogram(data, bins=len(data_sorted))  
    # finding the PDF of the histogram using count values
    pdf_obs = count / sum(count)
    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf_obs = np.cumsum(pdf_obs)    
    axt = ax.twinx()
    axt.plot(data_sorted, cdf_obs, linestyle='--', marker='',color='black')
    axt.plot(data_sorted, cdf_null_hyp,  linestyle='-',color=colors[n])
    """
    #get max different (K-S statistic)
    arr_dif_abs = np.abs(cdf_null_hyp-np.array(cdf_obs))
    dn_ks = max(arr_dif_abs)
    # Plot distances between cdf_null_hyp and observed data
    for x1,x2, y1, y2 in zip(data_sorted,data_sorted, cdf_obs, cdf_null_hyp):
        axt.plot([x1, x2], [y1, y2], color='green', alpha=0.4)
    """
    axt.spines.right.set_visible(True) #It was set as False by default in the .mpstyle file
    return axt
    
########################################################################
# Plot function (PDFs and CDFs)
def plot_dist(n_df, epi_dist, max_val, ax, n_subset = None, subset = 'wave',
              dist_list = ['gamma', 'lognormal', 'weibull', 'exponential'],
              bin_unit = 1, title=None):
    df = all_dfs[n_df]
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    best = df_best_models[df_best_models[epi_dist] == 0].index[0]
    if n_subset:
        data = df.loc[(df[subset ]== n_subset) & (df[epi_dist] <= max_val)][epi_dist]
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[subset] == n_subset) & (df[epi_dist] <= max_val)][epi_dist]
        ext = f'[{n_subset}]'
    else:
        data = df.loc[df[epi_dist]<=max_val][epi_dist]   
        if max_val > data.max(): 
            max_val = data.max()
            data = df.loc[(df[subset] == n_subset) & (df[epi_dist]<=max_val)][epi_dist]
        ext = ''
    data = data + 0.5
    data_sorted = np.sort(np.array(data.drop_duplicates()))
    params = dist_posteriors[epi_dist]      
    bins = int(max_val/bin_unit)
    ax.hist(data, density=True, bins = bins, alpha = 0.4, edgecolor = 'white')#, label=col)
    xx = np.linspace(min(data),max(data),1000)
    for n, dist in enumerate(dist_list):
        if dist == 'Gamma':
            a = params[dist]['alpha'+ext]
            b = params[dist]['beta'+ext]
            fit_stan = stats.gamma.pdf(xx, a = a, scale = 1/b, loc = 0)
            if dist == best:
                cdf_null_hyp = [stats.gamma.cdf(x, a = a, scale = 1/b, loc = 0) for x in data_sorted]
                axt = plot_cdf(df, epi_dist = epi_dist, max_val = max_val,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Lognormal':
            s = params[dist]['sigma'+ext]
            mu = params[dist]['mu'+ext]
            fit_stan = stats.lognorm.pdf(xx, s = s, scale = np.exp(mu))
            if dist == best:
                cdf_null_hyp = [stats.lognorm.cdf(x, s = s, scale = np.exp(mu)) for x in data_sorted]
                axt = plot_cdf(df, epi_dist = epi_dist, max_val = max_val,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Weibull':
            a = params[dist]['alpha'+ext]
            s = params[dist]['sigma'+ext]
            fit_stan = stats.weibull_min.pdf(xx, c = a, scale = s)
            if dist == best:
                cdf_null_hyp = [stats.weibull_min.cdf(x, c = a, scale = s) for x in data_sorted]
                axt = plot_cdf(df, epi_dist = epi_dist, max_val = max_val,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n=n, n_subset = n_subset)
        elif dist == 'Exponential':
            b = params[dist]['beta'+ext]
            fit_stan = stats.expon.pdf(xx, loc = 0, scale = 1/b)
            if dist == best:
                cdf_null_hyp = [stats.expon.cdf(x, loc = 0 , scale = 1/b) for x in data_sorted]
                axt = plot_cdf(df, epi_dist = epi_dist, max_val = max_val,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        elif dist == 'Gen Lognormal':
            mu = params[dist]['mu'+ext]
            s = params[dist]['sigma'+ext]
            g = params[dist]['g'+ext]
            fit_stan = ut.gln_pdf(xx, mu = mu, sigma = s, g = g)
            if dist == best:
                cdf_null_hyp = [ut.gln_cdf(x, mu = mu, sigma = s, g = g) for x in data_sorted]
                axt = plot_cdf(df, epi_dist = epi_dist, max_val = max_val,
                         ax = ax, cdf_null_hyp = cdf_null_hyp, n = n, n_subset = n_subset)
        #Fit
        ax.plot(xx, fit_stan, label = dist, color = colors[n])        
    ax.set_xlim([0,max_val])
    ax.set_xlabel(epi_dist.replace('_',' '))
    ax.set_ylabel('Probability')
    ax.set_title(title)
    mpla.align.yaxes(ax, 0, axt, 0, 0.03)
    return ax

########################################################################
########################################################################
########################################################################
# Plot best model
def get_best_error(dist):
    df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")
    df_dist = df_best_models_sum[df_best_models_sum.dist == dist]
    df_dist = df_dist.set_index('stat')
    name = dist.replace('_', ' ')
    mean = np.array([df_dist['wave_1']['mean'], df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['mean'], df_dist['wave_4']['mean']])
    err = abs(np.array([[df_dist['wave_1']['mean'] - df_dist['wave_1']['q025'], 
                     df_dist['wave_2']['mean'] - df_dist['wave_2']['q025'], 
                     df_dist['wave_3']['mean'] - df_dist['wave_3']['q025'], 
                     df_dist['wave_4']['mean'] - df_dist['wave_4']['q025']],
                    [df_dist['wave_1']['q975'] - df_dist['wave_1']['mean'], 
                     df_dist['wave_2']['q975'] - df_dist['wave_2']['mean'], 
                     df_dist['wave_3']['q975'] - df_dist['wave_3']['mean'], 
                     df_dist['wave_4']['q975'] - df_dist['wave_4']['mean']]]))
    return mean, err

def plot_best_model_line(dist, ax, n, title):
    df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")
    df_dist = df_best_models_sum[df_best_models_sum.dist == dist]
    df_dist = df_dist.set_index('stat')
    name = dist.replace('_', ' ')
    
    mean, err = get_best_error(dist)
    
    ax.errorbar(['1', '2', '3', '4'], mean, err, 
                ls = '-', marker = 'o', 
                color = colors[n],
                label = name)
    ax.set_xlabel('Wave')
    ax.set_ylabel('Days')
    ax.set_title(title)
    return ax

def plot_best_model_bar(dist, ax, n, name_y, title):
    df_dist = df_best_models_sum[df_best_models_sum.dist == dist]
    df_dist = df_dist.set_index('stat')
    name = dist.replace('_', ' ')
    
    mean, err = get_best_error(dist)
    
    ax.bar(['1', '2', '3', '4'], mean, yerr = err,
    color = colors[n], label = name)
    ax.set_xlabel('Wave')
    ax.set_ylabel(name_y)
    ax.set_title(title)
    
def plot_violin(var, name_y, title, ax):
    if var == 'icu_stay':
        n = 0
    if var == 'hosp_stay':
        n = 1
    if var == 'onset_icu':
        n = 2
    if var == 'onset_hosp':
        n = 3
    if var == 'onset_death':
        n = 4
    df = all_dfs[n]
    df = df[(df[var] > 0) & (df[var] <= 80)]
    df_best_models_sum = pd.read_csv(OUT_PATH + "best_fit_summary.csv")
    df_dist = df_best_models_sum[df_best_models_sum.dist == var]
    df_dist = df_dist.set_index('stat')
    df['wave'] = df['wave'].astype(int)
    sns.violinplot(data=df, x="wave", y=var, dodge=False, ax= ax, color = colors[n])
    ax.set_xlabel('Wave')
    ax.set_ylabel(name_y)
    ax.set_title(title)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

def plot_best_model_bar_all(dist, ax, w, n, wt=0.1):
    mean, err = get_best_error(dist)
    ax.bar([1+w, 2+w, 3+w, 4+w], mean, yerr = err, 
           width = abs(wt), 
           color = colors[n])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_yscale('log')
    ax.set_ylabel('Average value of delay time (Days)')
    ax.set_xlabel('Wave')
    
def plot_best_model_bar_outcome(wave, ax, w, wt=0.1):
    dist = ['onset_hosp', 'onset_icu', 'onset_death', 'hosp_stay', 'icu_stay']
    xlabels= ['Onset to hospitalisation', 
            'Onset to ICU entrance',
            'Onset to death',
            'Hospital stay',
            'ICU stay']
    mean = []
    err = [[],[]]
    for d in dist:
        m, e = get_best_error(d)
        mean.append(m[wave-1])
        err[0].append(e[0][wave-1])
        err[1].append(e[1][wave-1])
    mean = np.array(mean)
    err = np.array(err)
    
    n=wave-1
    ax.bar([1+w, 2+w, 3+w, 4+w, 5+w], mean, yerr = err, 
           width = abs(wt), 
           color = colors[n])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.set_yscale('log')
    ax.set_ylabel('Average value of delay time (Days)')
    ax.set_xlabel('Epidemiological distribution')
    ax.set_xticks([1,2,3,4,5], xlabels)