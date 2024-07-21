# -*- coding: utf-8 -*-
"""
Created on Th December 1 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""  

import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

config = yaml.load(open("config.yml", "r"))["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
sys.path.append(SCRIPTS_PATH)

import results_severe_outcomes as results_severe_outcomes

fig, ax = plt.subplots(2, 3, figsize=(12, 8), sharex = True)

results_severe_outcomes.plot_ratios(ax=ax[0][1], var='HCR', var_name='HCR')
results_severe_outcomes.plot_ratios(ax=ax[0][2], var='ICU-CR', var_name=' ICU-CR')

results_severe_outcomes.plot_ratios(ax=ax[1][0], var='CFR', var_name='CFR')
results_severe_outcomes.plot_ratios(ax=ax[1][1], var='HFR', var_name='HFR')
results_severe_outcomes.plot_ratios(ax=ax[1][2], var='ICU-FR', var_name='ICU-FR')

handles, labels = ax[1][1].get_legend_handles_labels()
ax[0][0].legend(handles, labels, loc='center')
ax[0][0].axis('off')

ax[0][1].tick_params(axis='x', rotation=90)
ax[0][2].tick_params(axis='x', rotation=90)
ax[1][0].tick_params(axis='x', rotation=90)
ax[1][1].tick_params(axis='x', rotation=90)
ax[1][2].tick_params(axis='x', rotation=90)

ax[0][1].set_title('a.')
ax[0][2].set_title('b.')
ax[1][0].set_title('c.')
ax[1][1].set_title('d.')
ax[1][2].set_title('e.')

ax[0][1].set_xlabel('')
ax[0][2].set_xlabel('')
ax[1][0].set_xlabel('Age group')
ax[1][1].set_xlabel('Age group')
ax[1][2].set_xlabel('Age group')

ax[0][1].set_ylabel('Hospitalisation Case Ratio')
ax[0][2].set_ylabel('ICU Case Ratio')
ax[1][0].set_ylabel('Case Fatality Ratio')
ax[1][1].set_ylabel('Hospitalisation Fatality Ratio')
ax[1][2].set_ylabel('ICU Fatality Ratio')

fig.savefig(FIG_PATH + 'figure_3.jpg')