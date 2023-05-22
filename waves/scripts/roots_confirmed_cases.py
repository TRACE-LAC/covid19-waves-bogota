# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

b=10 # Gaussian kernel width

config = yaml.load(open("config.yml", "r"))["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'waves')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'waves')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'waves')
UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']

plt.style.use(config['PATHS']['PLOT_STYLE'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

sys.path.append(UTILS_PATH)
import utilities_waves as ut

# Confirmed cases
df_confirmed_bogota = pd.read_csv(DATA_PATH + f'confirmed_cases_{UPDATE}.csv')
df_confirmed_bogota['onset'] = pd.to_datetime(df_confirmed_bogota['onset'], errors='coerce')

# Data processing
## Notice that we perform gaussian smoothing twice, before and after differentiate the epidemic curve
df_counts = ut.counts('onset',df_confirmed_bogota)
df_counts['cases_gs'] = ut.gaussian_smoothing(df_counts, 'cases', 'date', b)
df_counts['cases_gs_diff'] = df_counts.cases_gs.diff()  
df_counts = df_counts[df_counts.cases_gs_diff.notnull()]
df_counts['cases_gs_diff_gs'] = ut.gaussian_smoothing(df_counts, 'cases_gs_diff', 'date', b)

# Roots of the first diff
x = df_counts['date'].to_numpy()
y = df_counts['cases_gs_diff_gs'].to_numpy()
roots = ut.find_roots(x,y)

df_roots = pd.DataFrame(pd.Series(roots), columns = ['date'])
# df_roots['date'] = df_roots['date'].dt.date
df_roots['date'] = pd.to_datetime(df_roots['date'])
df_roots = df_roots.drop_duplicates().reset_index(drop = True)
df_roots.to_csv(OUT_PATH + 'roots_confirmed_cases.csv',index = False)