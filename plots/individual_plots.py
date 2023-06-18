# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""     
import warnings
warnings.filterwarnings('ignore')

import yaml
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from met_brewer import met_brew

config = yaml.load(open("config.yml", "r"))["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
sys.path.append(SCRIPTS_PATH)

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

#################################### ####################################
############################## Genomics #################################
#################################### ####################################
DATE_GENOMICS = config['UPDATE_DATES']['GENOMICS']
import results_genomics as results_genomics

# Multinomial
fig, ax = plt.subplots()
results_genomics.plot_multinomial(ax, ['2021-12', '2022-30'])
fig.savefig(FIG_PATH + 'variants_multinomial.png')

# Prevalence
fig, ax = plt.subplots(figsize = (7.5,5))
results_genomics.plot_prevalence(ax)
ax.set_xlabel('')
ax.set_ylabel('Prevalence')
ax.tick_params(axis = 'x', rotation = 90)
fig.savefig(FIG_PATH + 'variants_prevalence_' + DATE_GENOMICS + '.png')

#################################### ####################################
################################# RT ####################################
#################################### ####################################
import results_rt as results_rt

# rt
fig, ax = plt.subplots()
results_rt.plot_rt(ax)
ax.set_xlabel('')
ax.set_ylabel('R(t)') 
ax.set_ylim(top=2.3)
ax.tick_params(axis='x', rotation=90)
fig.savefig(FIG_PATH + 'rt.png')

#################################### ####################################
############################## Overview #################################
#################################### ####################################
import overview as overview

# Population pyramid
fig, ax = plt.subplots(figsize=(7.5, 5))
overview.plot_pyramid(ax)
ax.set_xlim(left=-260000)
ax.set_xlabel('Cases')
ax.set_ylabel('Age group')
ax.legend()     
fig.savefig(FIG_PATH + 'population_pyramid.png')

# Cases and cum. death
fig, ax = plt.subplots(figsize=(7.5, 5))
overview.plot_cases_death_cum(ax)
fig.savefig(FIG_PATH + 'cases_death_cum.png')

#################################### ####################################
############################### Waves ##################################
#################################### ####################################
import results_waves as results_waves
fig, ax = plt.subplots(figsize = (15,5))
results_waves.plot_waves(ax)
fig.savefig(FIG_PATH + 'waves.png')

#################################### ####################################
#########################  Epi distributions ############################
#################################### ####################################
import results_epidemiological_distributions as results_epidist
dist_list = ['Gamma', 'Lognormal', 'Weibull', 'Exponential', 'Gen Lognormal']
max_val_plot = 60
fig, ax = plt.subplots(5, 4, figsize = (10, 11))

#ONSET HOSP
var = 'onset_hosp'
n_df = 3 #df_onset_hosp
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='a.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='b.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='c.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[0][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='d.')

#ONSET ICU
var = 'onset_icu'
n_df = 2 #df_onset_icu
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='e.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='f.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='g.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[1][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='h.')

#ONSET DEATH
var = 'onset_death'
n_df = 4 #df_onset_death
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='i.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='j.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='k.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[2][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='l.')

#HOSP stay
var = 'hosp_stay'
n_df = 1 #df_hosp_stay
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[3][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='m.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[3][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='n.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[3][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='o.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[3][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='p.')

#ICU stay
var = 'icu_stay'
n_df = 0 #df_icu_stay
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[4][0], n_subset = 1, dist_list = dist_list,
          bin_unit = 1, title='q.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[4][1], n_subset = 2, dist_list = dist_list,
          bin_unit = 1, title='r.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[4][2], n_subset = 3, dist_list = dist_list,
          bin_unit = 1, title='s.')
results_epidist.plot_dist(n_df, epi_dist = var, max_val = max_val_plot, 
          ax = ax[4][3], n_subset = 4, dist_list = dist_list,
          bin_unit = 1, title='t.')

handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.99, -0.03), ncol = len(dist_list))
fig.savefig(FIG_PATH + 'distributions.png')

#################################### ####################################
 #########################  Severe outcomes ############################
#################################### ####################################
import results_severe_outcomes as results_severe_outcomes

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
results_severe_outcomes.plot_percentage_err(ax)
ax[0].set_ylabel('Hospitalization percentage by age-group')
ax[1].set_ylabel('ICU percentage by age-group')
ax[2].set_ylabel('Deaths percentage by age-group')
for axi in ax:
    axi.tick_params(axis='x', labelrotation=90)
    axi.set_xlabel('Age group')
ax[0].set_title('a.')
ax[1].set_title('b.')
ax[2].set_title('c.')
handles, labels = ax[2].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor = (0.79, -0.02), ncol = 6)
fig.savefig(FIG_PATH+'hosp_icu_death_percentages.png')

fig, ax = plt.subplots(1, 3, figsize=(12, 4))
results_severe_outcomes.plot_counts(ax)
handles, labels = ax[0].get_legend_handles_labels()
ax[0].set_title('a.')
ax[1].set_title('b.')
ax[2].set_title('c.')
fig.legend(handles, labels, bbox_to_anchor = (0.8, -0.03), ncol = len(labels)) 
fig.savefig(FIG_PATH+'hosp_icu_death_counts.png')

fig, ax = plt.subplots(1, 3, figsize=(12, 4))
for axi in ax:
    axi.set_xlabel('Wave')
ax[0].set_ylabel('Hospitalization counts by age')
ax[1].set_ylabel('ICU counts by age')
ax[2].set_ylabel('Death counts by age')
ax[0].set_title('a.')
ax[1].set_title('b.')
ax[2].set_title('c.')
results_severe_outcomes.plot_counts_histograms(ax)
fig.savefig(FIG_PATH+'hosp_icu_death_counts_hist.png') 

# CFR
fig, ax = plt.subplots()
var = 'CFR'
var_name = 'CFR'
results_severe_outcomes.plot_ratios(ax, var, var_name)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left',numpoints=1)
fig.savefig(FIG_PATH + f'{var}.png')

# HCR
fig, ax = plt.subplots()
var = 'HCR'
var_name = 'HCR'
results_severe_outcomes.plot_ratios(ax, var, var_name)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left',numpoints=1)
fig.savefig(FIG_PATH + f'{var}.png')

# HCR - I
fig, ax = plt.subplots()
var = 'ICU-CR'
var_name = 'ICU-CR'
results_severe_outcomes.plot_ratios(ax, var, var_name)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left',numpoints=1)
fig.savefig(FIG_PATH + f'{var}.png')

# HFR - I
fig, ax = plt.subplots()
var = 'HFR'
var_name = 'HFR'
results_severe_outcomes.plot_ratios(ax, var, var_name)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left',numpoints=1)
fig.savefig(FIG_PATH + f'{var}.png')

# HFR - I
fig, ax = plt.subplots()
var = 'ICU-FR'
var_name = 'ICU-FR'
results_severe_outcomes.plot_ratios(ax, var, var_name)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left',numpoints=1)
fig.savefig(FIG_PATH + f'{var}.png')