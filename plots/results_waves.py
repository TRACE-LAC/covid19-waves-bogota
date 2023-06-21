# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""     
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'waves')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'waves')
# The roots were selected by visual inspection of /figures/roots_confirmed_cases.png 
# and set by default in config.yml
b = config['WAVES']['GAUSSIAN_KERNEL'] 

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

sys.path.append(UTILS_PATH)
import utilities_waves as ut

# Read confirmed cases and waves information
df_confirmed_bogota = pd.read_csv(DATA_PATH + 'confirmed_cases_' + UPDATE + '.csv')
df_confirmed_bogota['onset'] = pd.to_datetime(df_confirmed_bogota['onset'], errors='coerce')

df_waves = pd.read_csv(OUT_PATH + 'waves.csv')
df_waves['start_date'] = pd.to_datetime(df_waves['start_date'], errors='coerce')
df_waves['end_date'] = pd.to_datetime(df_waves['end_date'], errors='coerce')

df_roots = pd.read_csv(OUT_PATH + 'roots_confirmed_cases.csv')
df_roots['date'] = pd.to_datetime(df_roots['date'])

# Calculate counts
df_counts = ut.counts('onset',df_confirmed_bogota)
df_counts['cases_gs'] = ut.gaussian_smoothing(df_counts, 'cases', 'date', b)
df_counts['cases_gs_diff'] = df_counts.cases_gs.diff()  
df_counts = df_counts[df_counts.cases_gs_diff.notnull()]
df_counts['cases_gs_diff_gs'] = ut.gaussian_smoothing(df_counts, 'cases_gs_diff', 'date', b)

# Function to plot waves figure
def plot_waves(ax):
    # Confirmed cases curve
    ln1 = ax.plot(df_counts['date'], df_counts['cases'], color = colors[0], label = 'confirmed cases')
    ln2 = ax.plot(df_counts['date'], df_counts['cases_gs'], color= colors[2], ls = '--')
    date_no_dup = []
    for d in range(len(df_waves)):
            ax.axvline(x = df_waves['start_date'].iloc[d], color = 'black', alpha = 0.8, ls = '--')
            ax.axvline(x = df_waves['end_date'].iloc[d], color = 'black', alpha = 0.8, ls = '--')
            date_no_dup.append(d)
    ax.set_ylabel('Confirmed cases')
    # First difference plot
    ax2 = ax.twinx()
    ln3 = ax2.plot(df_counts['date'],df_counts['cases_gs_diff'], colors[1], label = 'diff(confirmed cases)')
    ln4 = ax2.plot(df_counts['date'],df_counts['cases_gs_diff_gs'], colors[2], ls = '--', label = 'gaussian smoothing')

    ax2.axhline(y = 0, color = 'black', alpha = 0.8, ls = '--')
    ax2.set_ylabel('diff(Confirmed cases)')
    ax2.spines.right.set_visible(True) #This was set as False by default in the .mpstyle file

    ##Legend
    lns = ln1 + ln2 + ln3 + ln4
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc = 0)
    
def draw_waves(ax):
        date_no_dup = []
        for d in range(len(df_waves)):
                ax.axvline(x = df_waves['start_date'].iloc[d], color = 'grey', alpha = 0.5, ls = '--')
                ax.axvline(x = df_waves['end_date'].iloc[d], color = 'grey', alpha = 0.5, ls = '--')
                date_no_dup.append(d)