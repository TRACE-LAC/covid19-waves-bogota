# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 2022
@author: davidsantiagoquevedo
@author: ntorresd
""" 

import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

#Paths 
UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']
DATA_PATH = config['PATHS']['DATA_PATH']
# WE NEED TO UNIFY THE UTILITIES
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'severe_outcomes')

# Plot style
plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# Import useful functions
sys.path.append(UTILS_PATH)
import utilities_severity as ut

# Read data
df_confirmed_bogota = pd.read_csv(DATA_PATH + 'confirmed_cases_' + UPDATE + '.csv')
df_confirmed_bogota = df_confirmed_bogota.astype({'age':int})

df_confirmed_bogota['onset'] = pd.to_datetime(df_confirmed_bogota['onset'], errors='coerce')
df_confirmed_bogota['death'] = pd.to_datetime(df_confirmed_bogota['death'], errors='coerce')

# Population pyramid
age_group_dec_dic={
    '0-9':[0,9],
    '10-19':[10,19],
    '20-29':[20,29], 
    '30-39':[30,39], 
    '40-49':[40,49], 
    '50-59':[50,59],
    '60-69':[60,69],
    '70-79':[70,79], 
    '80+':[80,9999]
}
df_confirmed_bogota = ut.age_group_dec(df_confirmed_bogota,'age','age_unit', age_group_dec_dic)
df_counts = df_confirmed_bogota.groupby(by=['age_group', 'sex']).size().reset_index(name='counts')

df_female = df_counts[df_counts['sex']=='F']
df_male = df_counts[df_counts['sex']=='M']
df_prop = pd.crosstab(index=df_confirmed_bogota['age_group'],
                            columns=df_confirmed_bogota['sex'],
                            normalize="index")*100                
def plot_pyramid(ax):
        hist_f = ax.barh(df_female['age_group'], -df_female['counts'], align='center', label='female')
        hist_m = ax.barh(df_male['age_group'], df_male['counts'], align='center', label='male')

        ax.bar_label(hist_m, labels=['{:.1f}%'.format(prop) for prop in df_prop['M'].round(1)], label_type='edge')
        ax.bar_label(hist_f, labels=['{:.1f}%'.format(prop) for prop in df_prop['F'].round(1)], label_type='edge')

        ticks =  ax.get_xticks()
        ax.set_xticklabels([int(abs(tick)) for tick in ticks])

# Incidences 
df = df_confirmed_bogota
df = ut.age_group_60(df,'age','age_unit')
df.dropna(subset=['age_group'], inplace=True)

## Conditions
mask_60p = (df['age_group']=='60+')
mask_death = (df['death'].notnull()) & ~(df_confirmed_bogota['condition'] =='Fallecido (No aplica No causa Directa)')

df_cases_all = ut.counts(df, var='onset',columns=['date', 'cases'])
df_cases_60p = ut.counts(df[mask_60p], var='onset', columns=['date', 'cases'])

df_deaths_all = ut.counts(df[mask_death], var='death',columns=['date', 'deaths'])
df_deaths_all['deaths_cum'] = df_deaths_all['deaths'].cumsum().values

df_deaths_60p = ut.counts(df[mask_death & mask_60p], var='death',columns=['date', 'deaths'])
df_deaths_60p['deaths_cum'] = df_deaths_60p['deaths'].cumsum().values

def plot_cases_death_cum(ax):
        ln1 = ax.plot(df_cases_all['date'], df_cases_all['cases'], label='cases all')
        ln2 = ax.plot(df_cases_60p['date'], df_cases_60p['cases'], label='cases 60+')

        # fig, ax2 =plt.subplots(figsize=(7.5, 5))
        ax_twin = ax.twinx()

        ln3 = ax_twin.plot(df_deaths_all['date'], df_deaths_all['deaths_cum'], label='cumulative deaths all', linestyle='dashed')
        ln4 = ax_twin.plot(df_deaths_60p['date'], df_deaths_60p['deaths_cum'], label='cumulative deaths 60+', linestyle='dashed')

        lns = ln1 + ln2 + ln3 + ln4
        labs = [l.get_label() for l in lns]

        ax.tick_params(axis='x', rotation=90)
        ax.set_ylabel('Confirmed cases')

        ax_twin.legend(lns, labs, loc='upper left', frameon=False, fontsize=10)
        ax_twin.spines.right.set_visible(True) #This was set as False by default in the .mpstyle file
        ax_twin.tick_params(axis='x', rotation=90)
        ax_twin.set_ylabel('Deaths')