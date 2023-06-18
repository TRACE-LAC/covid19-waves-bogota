# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""  
import warnings
warnings.filterwarnings('ignore')

import sys
import yaml
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

SCRIPTS_PATH = config['PATHS']['PLOT_PATH']
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'plots')
sys.path.append(SCRIPTS_PATH)

import results_rt as results_rt
import results_genomics as results_genomics
import results_waves as results_waves
import overview as overview

fig, ax = plt.subplots(2,2, figsize = (14,8))

axi = ax[0][0]
overview.plot_cases_death_cum(ax = axi)
results_waves.draw_waves(axi)
axi.set_title('a.')

axi = ax[0][1]
overview.plot_pyramid(ax = axi)
axi.set_xlim(left=-270000)
axi.set_xlabel('Cases')
axi.set_ylabel('Age group')
axi.set_title('b.')
axi.legend(loc = 'upper left')  

axi = ax[1][0]
results_rt.plot_rt(ax = axi)
results_waves.draw_waves(axi)
axi.set_xlabel('')
axi.set_ylabel('R(t)') 
axi.set_ylim(top = 2.4)
axi.tick_params(axis = 'x', rotation = 90)
axi.set_title('c.')

axi = ax[1][1]
results_genomics.plot_prevalence(ax = axi)
axi.tick_params(axis = 'x', rotation = 90)
axi.set_xlabel('')
axi.set_ylabel('Prevalence')
axi.set_title('d.', loc = 'left')

fig.savefig(FIG_PATH + 'figure_1.png')