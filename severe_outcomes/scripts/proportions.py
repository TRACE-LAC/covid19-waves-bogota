# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint

config = yaml.load(open("config.yml", "r"))["default"]

# Paths
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
OUT_PATH_WAVES = config['PATHS']['OUT_PATH'].format(dir = 'waves')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'severe_outcomes')
UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']

# Import useful functions
sys.path.append(UTILS_PATH)
import utilities_severity as ut
df_confirmed_bogota = pd.read_csv(DATA_PATH + f'confirmed_cases_waves_{UPDATE}.csv')
df_hosp = pd.read_csv(DATA_PATH+'hosp_waves_bog.csv')
df_icu = pd.read_csv(DATA_PATH+'icu_waves_bog.csv')
df_death = pd.read_csv(DATA_PATH+'death_waves_bog.csv')

df_waves = pd.read_csv(OUT_PATH_WAVES+'waves.csv')

# Changing the age groups
age_group_dic = {
    '0-9':'<60',
    '10-19':'<60',
    '20-29':'<60', 
    '30-39':'<60', 
    '40-49':'<60', 
    '50-59':'<60',
    '60-69':'60+',
    '70-79':'60+', 
    '80+':'60+'
}
df_hosp = df_hosp.replace({'age_group':age_group_dic})
df_icu = df_icu.replace({'age_group':age_group_dic})
df_death = df_death.replace({'age_group':age_group_dic})

# df_confirmed_bogota = ut.get_wave(df_confirmed_bogota, start_date='onset', waves=df_waves)
# %%
df_confirmed_bogota = ut.age_group_60(df=df_confirmed_bogota, var='age', var_unit='age_unit')
df_confirmed_bogota.dropna(subset=['wave'], inplace=True)

data_all = []
data_60p = []

strat='wave'
# Counting cases/hosp/icu/deat by wave for all the population  
cases_all = ut.size_by_strat(df_confirmed_bogota, strat=strat)
hosp_all = ut.size_by_strat(df_hosp, strat=strat)
icu_all = ut.size_by_strat(df_icu, strat=strat)
death_all = ut.size_by_strat(df_death, strat=strat)

# Counting cases/hosp/icu/deat by wave for the population older than 60 years 
cases_60p = ut.size_by_strat(df_confirmed_bogota[df_confirmed_bogota['age_group']=='60+'], strat=strat)
hosp_60p = ut.size_by_strat(df_hosp[df_hosp['age_group']=='60+'], strat=strat)
icu_60p = ut.size_by_strat(df_icu[df_icu['age_group']=='60+'], strat=strat)
death_60p = ut.size_by_strat(df_death[df_death['age_group']=='60+'], strat=strat)

# Computing the confidence interval as a binomial proportion
alpha = 0.05 #significance level

hosp_all_err = proportion_confint(count=hosp_all, nobs=cases_all, alpha=alpha)
icu_all_err = proportion_confint(count=icu_all, nobs=cases_all, alpha=alpha)
death_all_err = proportion_confint(count=death_all, nobs=cases_all, alpha=alpha)

hosp_60p_err = proportion_confint(count=hosp_60p, nobs=cases_60p, alpha=alpha)
icu_60p_err = proportion_confint(count=icu_60p, nobs=cases_60p, alpha=alpha)
death_60p_err = proportion_confint(count=death_60p, nobs=cases_60p, alpha=alpha)

# Constructing the dataframes with the results
data_all = {
    'wave':df_waves['wave'],
    'cases':cases_all, 
    'hosp':hosp_all/cases_all, 
    'hosp_lower':hosp_all_err[0].tolist(), 
    'hosp_upper':hosp_all_err[1].tolist(),
    'icu':icu_all/cases_all, 
    'icu_lower':icu_all_err[0].tolist(), 
    'icu_upper':icu_all_err[1].tolist(),
    'death':death_all/cases_all, 
    'death_lower':death_all_err[0].tolist(), 
    'death_upper':death_all_err[1].tolist()
    }
df_counts_all = pd.DataFrame(data=data_all, dtype=float)
df_counts_all.to_csv(OUT_PATH + 'proportions_all.csv',index = False)

data_60p = {
    'wave':df_waves['wave'],
    'cases':cases_60p, 
    'hosp':hosp_60p/cases_60p, 
    'hosp_lower':hosp_60p_err[0].tolist(), 
    'hosp_upper':hosp_60p_err[1].tolist(),
    'icu':icu_60p/cases_60p, 
    'icu_lower':icu_60p_err[0].tolist(), 
    'icu_upper':icu_60p_err[1].tolist(),
    'death':death_60p/cases_60p, 
    'death_lower':death_60p_err[0].tolist(), 
    'death_upper':death_60p_err[1].tolist()
    }
df_counts_60p = pd.DataFrame(data=data_60p, dtype=float)
df_counts_60p.to_csv(OUT_PATH + 'proportions_60p.csv',index = False)