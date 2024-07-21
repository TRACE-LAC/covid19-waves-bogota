# -*- coding: utf-8 -*-
"""
Created on Wed 1 Dic 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""
import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
sys.path.append(SCRIPTS_PATH)

import results_epidemiological_distributions as red

########### ########### ###########
###########   Option 0  ###########
########### ########### ###########
fig, ax = plt.subplots(figsize = (15,5))
red.plot_best_model_bar_outcome(1, ax, w = - 0.15)
red.plot_best_model_bar_outcome(2, ax, w = - 0.05)
red.plot_best_model_bar_outcome(3, ax, w =  0.05)
red.plot_best_model_bar_outcome(4, ax, w =  0.15)

leglabels= ['Wave 1',
            'Wave 2',
            'Wave 3',
            'Wave 4',
           ]
ax.legend(leglabels, bbox_to_anchor=(0.75, -0.15), ncol=5)
fig.savefig(FIG_PATH + 'figure_5_v0.jpg')
########### ########### ###########
###########   Option 1  ###########
########### ########### ###########
fig, ax = plt.subplots(figsize = (15,5))

dist = 'onset_hosp'
w=-0.2
red.plot_best_model_bar_all(dist, ax, w, 1)

dist = 'onset_icu'
w=-0.1
red.plot_best_model_bar_all(dist, ax, w, 2)

dist = 'onset_death'
w= 0
red.plot_best_model_bar_all(dist, ax, w, 3)

dist = 'hosp_stay'
w=0.1
red.plot_best_model_bar_all(dist, ax, w, 4)
dist = 'icu_stay'
w=0.2
red.plot_best_model_bar_all(dist, ax, w, 5)

## Legends
labels= ['Onset to hospitalisation', 
           'Onset to ICU entrance',
           'Onset to death',
           'Hospital stay',
           'ICU stay']
ax.legend(labels, bbox_to_anchor=(0.95, -0.15), ncol=5)

fig.tight_layout()
fig.show()
fig.savefig(FIG_PATH + 'figure_5_v1.jpg')

########### ########### ###########
###########   Option 2  ###########
########### ########### ###########
fig, ax = plt.subplots(5,2, figsize = (10, 14), sharex = True)
dist = 'onset_hosp'
red.plot_best_model_bar(dist, ax[0][0], 3, 'Onset to hospitalization', 'a.')
red.plot_violin(dist, 'Onset to hospitalization', 'f.', ax[0][1])

dist = 'onset_icu'
red.plot_best_model_bar(dist, ax[1][0], 2, 'Onset to ICU entrance', 'b.')
red.plot_violin(dist, 'Onset to ICU entrance', 'g.', ax[1][1])

dist = 'onset_death'
red.plot_best_model_bar(dist, ax[2][0], 4, 'Onset to death', 'c.')
red.plot_violin(dist, 'Onset to death', 'h.', ax[2][1])

dist = 'hosp_stay'
red.plot_best_model_bar(dist, ax[3][0], 1, 'Hospital stay', 'd.')
red.plot_violin(dist, 'Hospital stay', 'i.', ax[3][1])

dist = 'icu_stay'
red.plot_best_model_bar(dist, ax[4][0], 0, 'ICU stay', 'e.')
red.plot_violin(dist, 'ICU stay', 'j.', ax[4][1])

ax[0][0].set_xlabel('')
ax[0][1].set_xlabel('')
ax[0][1].set_xlabel('')
ax[1][0].set_xlabel('')
ax[1][1].set_xlabel('')
ax[2][0].set_xlabel('')
ax[2][1].set_xlabel('')
ax[3][0].set_xlabel('')
ax[3][1].set_xlabel('')

ax[0][1].set_ylabel('')
ax[1][1].set_ylabel('')
ax[2][1].set_ylabel('')
ax[3][1].set_ylabel('')
ax[4][1].set_ylabel('')

fig.savefig(FIG_PATH + 'figure_5_v2.jpg')