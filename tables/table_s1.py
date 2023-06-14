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
import pandas as pd

config = yaml.load(open("config.yml", "r"))["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
WAVES_PATH = config['PATHS']['OUT_PATH'].format(dir = 'waves')
RT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'rt')
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'tables')
UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']

# Waves information
df_waves = pd.read_csv(WAVES_PATH + "waves.csv")
df_waves.set_index("wave", inplace = True)
# Confirmed cases waves
df_confirmed = pd.read_csv(DATA_PATH + f'confirmed_cases_waves_{UPDATE}.csv')
# Outcomes waves
df_hosp = pd.read_csv(DATA_PATH + "hosp_waves_bog.csv")
df_icu = pd.read_csv(DATA_PATH + "icu_waves_bog.csv")
df_death = pd.read_csv(DATA_PATH + "death_waves_bog.csv")
# Rt results
df_rt = pd.read_csv(RT_PATH + "rt_all_ages.csv")

# Counts by wave
def get_size(df):
    df_size = df[["wave"]].groupby("wave").size()    
    return df_size.sort_index()

# Rt by wave
def get_rt_wave(ref_date = "window_start"):
    rt_wave = []    
    for wave in df_waves.index:
        st = df_waves["start_date"].loc[wave]
        ed = df_waves["end_date"].loc[wave]
        df_temp = df_rt[df_rt[ref_date].between(st,ed)]
        max = round(df_temp["Mean(R)"].max(), 2)
        q1 = round(df_temp["Quantile.0.025(R)"].max(), 2)
        q2 = round(df_temp["Quantile.0.975(R)"].max(), 2)
        rt_wave.append(str(max) + "(" + str(q1) + ", " + str(q2) + ")")
    return rt_wave

df_waves["Confirmed cases"] = get_size(df_confirmed)
df_waves["Hospitalizations"] = get_size(df_hosp)
df_waves["ICU"] = get_size(df_icu)
df_waves["Deaths"] = get_size(df_death)
df_waves["R(t)"] = get_rt_wave()


df_waves["start_date"] = pd.to_datetime(df_waves["start_date"]).dt.strftime("%Y-%m-%d")
df_waves["end_date"] = pd.to_datetime(df_waves["end_date"]).dt.strftime("%Y-%m-%d")

df_waves = df_waves.rename(columns = {"wave" : "Wave",
                                      "start_date" : "Start",
                                      "end_date" : "End"
                                      })

df_waves.to_csv(OUT_PATH + "table_s1.csv")