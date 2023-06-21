# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 2022
@author: davidsantiagoquevedo
@author: ntorresd
"""     
import yaml
import sys
import pandas as pd
import datetime
from datetime import datetime
import matplotlib.pyplot as plt
import mpl_axes_aligner as mpla
import seaborn as sns
from met_brewer import met_brew

config = yaml.load(open("config.yml", "r"))["default"]

UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'genomics')
DATE_GENOMICS = config['UPDATE_DATES']['GENOMICS']
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'severe_outcomes')

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# Import useful functions
sys.path.append(UTILS_PATH)
import utilities_severity as ut

# Load data
## Variants
df_variants = pd.read_csv(DATA_PATH + 'variants-ic-bog_'+ DATE_GENOMICS + '.csv')
df_variants['date'] = pd.to_datetime(df_variants['date'])
## Results from multinomial analysis
df_results = pd.read_csv(OUT_PATH + 'theta.csv')
df_pivot =  pd.read_csv(DATA_PATH + 'variants_pivot.csv').rename(columns = {'t' : 'week'})

df_results = df_results.merge(df_pivot[['week', 'week_name']], on = 'week')
df_results['week_date'] = df_results.week_name.apply(lambda date: datetime.strptime(date + '-1', "%Y-%W-%w"))
df_results = df_results.sort_values(by = 'week_date')
## Cases
df_confirmed_bogota = pd.read_csv(DATA_PATH + 'confirmed_cases_' + UPDATE + '.csv')
df_confirmed_bogota = df_confirmed_bogota.astype({'age':int})
df_confirmed_bogota['onset'] = pd.to_datetime(df_confirmed_bogota['onset'], errors='coerce')
df_confirmed_bogota['death'] = pd.to_datetime(df_confirmed_bogota['death'], errors='coerce')

n_ticks = 10     
        
def plot_multinomial(ax, limits):
        
        # Weekly counts - Incidence
        date1 = datetime.strptime(limits[0] + '-1', '%Y-%W-%u')
        date2 = datetime.strptime(limits[1] + '-1', '%Y-%W-%u')
        df_incidence = df_confirmed_bogota[df_confirmed_bogota['onset'].between(date1,date2)]
        df_incidence['week'] = df_confirmed_bogota['onset'].apply(lambda date: date.strftime('%Y-%W'))
        df_incidence_count = ut.counts(df_incidence, var='week',columns=['week', 'cases'])
        df_incidence_count['week_date'] = df_incidence_count['week'].apply(lambda date: datetime.strptime(date + '-1', '%Y-%W-%u'))
        df_incidence_count = df_incidence_count.sort_values(by = 'week_date')
        
        l1 = ax.bar(df_incidence_count['week_date'], df_incidence_count['cases'], 
               alpha = 0.3, color = "grey", width=5.3, label = 'Weekly new cases')
        ax.set_ylabel('Incidence')
        ax.tick_params(axis='x', rotation=90)
        ax.xaxis.set_major_locator(plt.MaxNLocator(n_ticks))
        labs = ['Weekly new cases'] 
        
        # Prevalence - Observed        
        ax1 = ax.twinx()
        stat = "mean"
        
        n = 0 
        variant = 'Alpha'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant, marker = '')
        ax1.fill_between(df_results[mask].week_date, df_results[mask1].theta, df_results[mask2].theta,
                        color  = colors[n], alpha = 0.2)
        labs.append(variant)
        
        n = 1 
        variant = 'Delta'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant, marker = '')
        ax1.fill_between(df_results[mask].week_date, df_results[mask1].theta, df_results[mask2].theta,
                        color  = colors[n], alpha = 0.2)
        labs.append(variant)
        
        n = 2
        variant = 'Gamma'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant, marker = '')
        ax1.fill_between(df_results[mask].week_date, df_results[mask1].theta, df_results[mask2].theta,
                        color  = colors[n], alpha = 0.2)
        labs.append(variant)
        
        n = 3
        variant = 'Mu'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant, marker = '')
        ax1.fill_between(df_results[mask].week_date, df_results[mask1].theta, df_results[mask2].theta,
                        color  = colors[n], alpha = 0.2)
        labs.append(variant)
        
        n = 4
        variant = 'Omicron'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant, marker = '')
        ax1.fill_between(df_results[mask].week_date, df_results[mask1].theta, df_results[mask2].theta,
                        color  = colors[n], alpha = 0.2)
        labs.append(variant)
        
        # Multinomial model
        n = 0 
        variant = 'Alpha'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
                color  = colors[n], marker = '^', linestyle = '')
        
        n = 1
        variant = 'Delta'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        ax1.plot(df_results[mask].week_date, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
                color  = colors[n], marker = '^', linestyle = '')
        n = 2
        variant = 'Gamma'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        l4 = ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant)
        ax1.plot(df_results[mask].week_date, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
                color  = colors[n], marker = '^', linestyle = '')
    
        n = 3
        variant = 'Mu'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        l5 = ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant)
        ax1.plot(df_results[mask].week_date, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
                color  = colors[n], marker = '^', linestyle = '')
        
        n = 4
        variant = 'Omicron'
        mask = (df_results.stat == stat) & (df_results.variant == variant)
        mask1 = (df_results.stat == 'q025') & (df_results.variant == variant)
        mask2 = (df_results.stat == 'q975') & (df_results.variant == variant)
        l6 = ax1.plot(df_results[mask].week_date, df_results[mask].theta, color  = colors[n], label = variant)
        ax1.plot(df_results[mask].week_date, df_results[mask][variant]/df_results[mask].weekly_count_variants, 
                color  = colors[n], marker = '^', linestyle = '')
                
        
        ax1.tick_params(axis='x', rotation=90)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(n_ticks))
        ax1.spines.right.set_visible(True)
        ax1.set_ylabel('Prevalence', rotation = 270, labelpad = 15)
        
        handles = ax.containers + ax1.get_lines()
        labels1, handles1 = ax1.get_legend_handles_labels()
        
        ax1.legend(handles, labs, loc='upper right', frameon=False, fontsize=12, ncol = 2)
        mpla.align.yaxes(ax, 0, ax1, 0, 0.03)

# Prevalence histogram
def plot_prevalence(ax):
        agg_df = df_variants.pivot_table(index = 'date', 
                                        columns = 'lineage', 
                                        values = 'PointEst', 
                                        aggfunc = 'max').reset_index()
        initial_date = min(df_variants['date'])
        final_date = max(df_variants['date'])
        bins = len(df_variants['date'].unique())
        variants_hist = sns.histplot(data = df_variants, ax=ax, 
                                     multiple = 'stack', 
                                     weights = 'PointEst', bins = bins, x = 'date', hue = 'lineage', 
                                     legend = True)
        ax.set_xticks(pd.date_range(start = min(df_variants['date']), end = final_date, freq = "M"))
        ax.set_xlim(left = min(df_variants['date']), right = final_date)
        sns.move_legend(variants_hist, 'lower right', bbox_to_anchor = (0.92, 0.95), ncol = 3, title=None)

# Advantage heatmap
def plot_heatmap(ax, n):
        df_mean = pd.read_csv(OUT_PATH + 'advantage_mean.csv')
        df_mean = df_mean.set_index('pivot_variant')
        color_list = met_brew('Cassatt1', n = n, brew_type='continuous')
        sns.heatmap(data = df_mean, annot = False, cmap = color_list, ax = ax)
        ax.set_ylabel('')
        ax.tick_params(axis='x', rotation=90)