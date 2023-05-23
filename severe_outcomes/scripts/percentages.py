# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt

config = yaml.load(open("config.yml", "r"))["default"]

# Paths
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'severe_outcomes')

# Import useful functions
sys.path.append(UTILS_PATH)
import utilities_severity as ut

# Read data
df_hosp = pd.read_csv(DATA_PATH+'hosp_waves_bog.csv')
df_icu = pd.read_csv(DATA_PATH+'icu_waves_bog.csv')
df_death = pd.read_csv(DATA_PATH+'death_waves_bog.csv')

round = 3
# Calculate counts and percentages
df_hosp_perc = ut.calculate_percentage(df_hosp, perc_var_name = 'HOSP-%', nobs = 'hosp', round = round)
df_icu_perc = ut.calculate_percentage(df_icu, perc_var_name = 'ICU-%', nobs = 'icu', round = round)
df_death_perc = ut.calculate_percentage(df_death, perc_var_name = 'DEATH-%', nobs = 'deaths', round = round)  

# Calculate binomial confidence intervals
df_hosp_perc = ut.calculate_confint(df_hosp_perc, var_name = 'HOSP-%', nobs = 'hosp', round = round)
df_icu_perc = ut.calculate_confint(df_icu_perc, var_name = 'ICU-%', nobs = 'icu', round = round)
df_death_perc = ut.calculate_confint(df_death_perc, var_name = 'DEATH-%', nobs = 'deaths', round = round)

# Merge all outcomes
df_percentages = df_hosp_perc.drop(columns = ['counts', 'hosp'])\
                            .merge(df_icu_perc.drop(columns = ['counts', 'icu']), how = 'left', on = ['wave', 'age_group'])\
                            .merge(df_death_perc.drop(columns = ['counts', 'deaths']), how = 'left', on = ['wave', 'age_group'])

# Save data
#TODO: Create table with total counts per wave for all outcomes
df_percentages.to_csv(OUT_PATH+'percentages.csv', index=False)