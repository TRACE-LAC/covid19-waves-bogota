# -*- coding: utf-8 -*-
"""
Created on Thr Jul 31 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""

import yaml
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from met_brewer import met_brew #  Package installation : https://github.com/BlakeRMills/MetBrewer

config = yaml.load(open("config.yml", "r"))["default"]

# Paths
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')

# Plot style
plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# Read data 
strat = 'wave'
# percentages
df_percentages = pd.read_csv(OUT_PATH + 'percentages.csv')
df_percentages[strat] = df_percentages[strat].astype(int)
# proportions
df_proportions_all = pd.read_csv(OUT_PATH+'proportions_all.csv')
df_proportions_60p = pd.read_csv(OUT_PATH+'proportions_60p.csv')
# rates 
df_ratios = pd.read_csv(OUT_PATH+'ratios.csv')
df_ratios = df_ratios[~df_ratios['age_group'].isin(['all'])].sort_values(by = ['wave', 'age_group'])

# Auxiliar plot function
def plot_xyvar(df, ax, n_strat, varx='age_group', vary='percentage'):
    data = df.loc[df[strat]==n_strat]
    ax.plot(data[varx], data[vary], label=strat+str(n_strat), ls = '-', marker = ".", lw = 1)

# Wave cases percentage distribution by age group
def plot_percentage(ax):
    vary = 'percentage'
    strat_list = df_percentages[strat].unique()
    for n_strat in strat_list:
        plot_xyvar(df_percentages[[strat, 'age_group', 'HOSP-%']], ax=ax[0], n_strat=n_strat, vary='HOSP-%')
        plot_xyvar(df_percentages[[strat, 'age_group', 'ICU-%']], ax=ax[1], n_strat=n_strat, vary='ICU-%')        
        plot_xyvar(df_percentages[[strat, 'age_group', 'DEATH-%']], ax=ax[2], n_strat=n_strat, vary='DEATH-%')

# Wave counts by age group
def plot_counts(ax):
    strat_list = df_ratios[strat].unique()
    for axi in ax:
        axi.tick_params(axis='x', labelrotation=90)
        axi.set_ylabel('cases')
        axi.set_xlabel('age group')
    for wave in strat_list:
        plot_xyvar(df_ratios[[strat, 'age_group','hosp']], ax=ax[0], n_strat=wave, vary='hosp')
        plot_xyvar(df_ratios[[strat, 'age_group', 'icu']], ax=ax[1], n_strat=wave, vary='icu')
        plot_xyvar(df_ratios[[strat, 'age_group', 'deaths']], ax=ax[2], n_strat=wave, vary='deaths')

# Stacked histogram of cases by wave and age group
def plot_stacked_histogram(df, ax, strat=strat, group_var='age_group', vary='counts', pallete='VanGogh1'): 
    strat_list = df[strat].unique()
    group_list = df[group_var].unique()
    color_list = met_brew(pallete, n=len(group_list), brew_type='continuous')
    
    y_offset = np.zeros(len(strat_list))
    counter=0
    for group in group_list:
        mask = (df[group_var]==group)
        counts_list =  df.loc[mask][vary].values
        ax.bar(strat_list, counts_list, label=group, bottom=y_offset, color=color_list[counter]) 
        y_offset = y_offset + counts_list 
        counter+=1

def plot_counts_histograms(ax, strat=strat, group_var='age_group', vary='counts'):

    plot_stacked_histogram(df = df_ratios[[strat, 'age_group','hosp']], ax=ax[0], strat=strat, group_var=group_var, vary='hosp')
    plot_stacked_histogram(df = df_ratios[[strat, 'age_group','icu']], ax=ax[1], strat=strat, group_var=group_var, vary='icu')
    plot_stacked_histogram(df = df_ratios[[strat, 'age_group','deaths']], ax=ax[2], strat=strat, group_var=group_var, vary='deaths')


# Proportion histogram with binomial confidence interval
# Auxiliar plot function
def plot_proportions_histogram(df, ax, outcome_var, prop_lower, prop_upper, strat=strat, group='all',side='center', width=0.3):
    proportions = df[outcome_var].values
    yerr_lower = proportions - df[prop_lower].values
    yerr_upper = df[prop_upper].values - proportions
    yerr = [yerr_lower, yerr_upper]
    if side=='center':
        ax.bar(df[strat], proportions, yerr=yerr, width=width, label=group)
    elif side=='left':
        ax.bar(df[strat]-width/2, proportions, yerr=yerr, width=width, label=group)
    elif side=='right':
        ax.bar(df[strat]+width/2, proportions, yerr=yerr, width=width, label=group)
    
def plot_proportions_hist(ax):
    plot_proportions_histogram(df_proportions_all, ax[0], 
                        outcome_var='hosp', prop_lower='hosp_lower', prop_upper='hosp_upper', side='left', group='all')
    plot_proportions_histogram(df_proportions_60p, ax[0], 
                        outcome_var='hosp', prop_lower='hosp_lower', prop_upper='hosp_upper', side='right', group='60+')

    plot_proportions_histogram(df_proportions_all, ax[1], 
                        outcome_var='icu', prop_lower='icu_lower', prop_upper='icu_upper', side='left', group='all')
    plot_proportions_histogram(df_proportions_60p, ax[1], 
                        outcome_var='icu', prop_lower='icu_lower', prop_upper='icu_upper', side='right', group='60+')

    plot_proportions_histogram(df_proportions_all, ax[2], 
                        outcome_var='death', prop_lower='death_lower', prop_upper='death_upper', side='left', group='all')
    plot_proportions_histogram(df_proportions_60p, ax[2], 
                        outcome_var='death', prop_lower='death_lower', prop_upper='death_upper', side='right', group='60+')

# Rates
# TODO: Use plot_var_err instead
def plot_ratios(ax, var, var_name):
    for wave in range(1,5):
        df_temp = df_ratios[df_ratios['wave'] == wave]
        ratio = df_temp[var] 
        ratio_lower = ratio - df_temp[var+'_lower']
        ratio_upper = df_temp[var+'_upper'] - ratio
        yerr = np.array(list(zip(ratio_lower, ratio_upper))).T
        ax.errorbar(df_temp['age_group'], df_temp[var], yerr = yerr,
                    fmt='o', color = colors[wave-1], label = 'wave '+str(wave),
                    capsize = 5, ls = '-', marker = ".", lw = 1)
    ax.set_xlabel('age group')
    ax.set_ylabel(var_name)

# Plot variable with error bar
def plot_var_err(df, axi, var):
    for wave in range(1, 5):
        df_temp = df[df['wave'] == wave]
        data = df_temp[var]
        data_lower = data - df_temp[var + '_lower']
        data_upper = df_temp[var + '_upper'] - data
        yerr = np.array(list(zip(data_lower, data_upper))).T
        axi.errorbar(df_temp['age_group'], df_temp[var], yerr = yerr,
                     fmt='o', color = colors[wave-1], label = 'wave '+str(wave),
                     capsize = 5, ls = '-', marker = ".", lw = 1)
        
# Plot percentages with erros
def plot_percentage_err(ax):
    plot_var_err(df_percentages, axi=ax[0], var = 'HOSP-%')
    plot_var_err(df_percentages, axi=ax[1], var = 'ICU-%')
    plot_var_err(df_percentages, axi=ax[2], var = 'DEATH-%')
