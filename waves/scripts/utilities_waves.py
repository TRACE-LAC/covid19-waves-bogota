# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import yaml
import pandas as pd
import numpy as np

config = yaml.load(open("config.yml", "r"))["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'waves')

# Count of confirmed cases for epidemic curve
def counts(var,df):
    df_counts = df[var].value_counts().reset_index()
    df_counts.columns = ['date','cases']
    df_counts = df_counts.sort_values(by = 'date')
    df_counts = df_counts.reset_index(drop = True)
    df_counts = df_counts.reset_index()
    return df_counts

# Gaussian smoothing https://towardsdatascience.com/gaussian-smoothing-in-time-series-data-c6801f8a4dc3
def gaussian_smoothing(df, var_cases, var_date, b):
    smoothed_cases = []
    for date in sorted(df[var_date]):
        df['gkv'] = np.exp(
            -(((df[var_date] - date).apply(lambda x: x.days)) ** 2) / (2 * (b ** 2))
        )
        df['gkv'] /= df['gkv'].sum()
        smoothed_cases.append(round(df[var_cases] * df['gkv']).sum())
    del(df['gkv'])
    return smoothed_cases

# Roots - Simple linear interpolation
def find_roots(x,y):
    s = np.abs(np.diff(np.sign(y))).astype(bool)
    return x[:-1][s] + np.diff(x)[s]/(np.abs(y[1:][s]/y[:-1][s])+1)