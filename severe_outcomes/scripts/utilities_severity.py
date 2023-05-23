# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import numpy as np

# calculate counts and percentages by wave and age_group
def calculate_percentage(df, perc_var_name, nobs, strat='wave', group='age_group', print_sum=False, round=4):
    df_percentage = df.groupby(by=[strat, group]).size().reset_index(name='counts')
    df_percentage[nobs] = 0
    df_percentage[perc_var_name] = 0
    strat_list = df_percentage[strat].unique()
    for strat_ in strat_list:
        mask = (df_percentage[strat]==strat_)
        df_percentage.loc[mask, nobs] = df_percentage.loc[mask, 'counts'].sum()
        df_percentage.loc[mask, perc_var_name] = 100 * df_percentage.loc[mask, 'counts'] / df_percentage.loc[mask, nobs]
        if print_sum:
            print('The sum over the percentages gives: ', df_percentage.loc[mask ,perc_var_name].sum())
    df_percentage[perc_var_name] = df_percentage[perc_var_name].round(round)
    return df_percentage

def get_pp_wave(ref_dt, waves):
    df_wave=waves[(waves['start_date']<=ref_dt)&(waves['end_date']>ref_dt)]
    if len(df_wave)>0:
        wave=df_wave.wave.iloc[0]
    else:
        wave=np.nan
    return wave

def get_wave(df, start_date, waves):
    df['wave']=df[start_date].apply(lambda x: get_pp_wave(x, waves=waves))
    return df

def age_group_60(df,var,var_unit): 
    conditions=[((df[var]<60) & (df[var_unit]==1))|(df[var_unit].isin([2,3]))] + [(df[var]>=60) & (df[var_unit]==1)]
    choices=['<60', '60+']
    df['age_group']=np.select(conditions,choices,default=np.nan)
    return df  

def age_group_dec(df,var,var_unit,dic): 
    conditions=[(df[var]<10)|(df[var_unit].isin([2,3]))] + [(df[var]>= dic[ge][0]) & (df[var] <= dic[ge][1]) for ge in list(dic.keys())[1:]]
    choices=[ge for ge in list(dic.keys())]
    df['age_group']=np.select(conditions,choices,default=np.nan)
    return df  

def size_by_strat(df, strat='wave'):
    strat_list = sorted(df[strat].unique())
    data = []
    for strat_ in strat_list:
        size = len(df[df[strat]==strat_].index)
        data.append(size)
    return np.array(data)

# Count of confirmed cases for epidemic curve
def counts(df, var='onset', columns=['date', 'counts']):
    df_counts = df[var].value_counts().reset_index()
    df_counts.columns = columns
    df_counts = df_counts.sort_values(by=columns[0])
    return df_counts

def cumulative(df, var='onset', columns=['date', 'counts']):
    df_counts = df[var].value_counts().reset_index()
    df_counts.columns = columns
    df_counts = df_counts.sort_values(by=columns[0])

    df_cum = df_counts.cumsum()

# Calculate binomial confidence interval for the specified stratification
from statsmodels.stats.proportion import proportion_confint
def calculate_confint(df, var_name, nobs, counts='counts', strat='wave', alpha=0.95, method='normal', round=4):
    prop = proportion_confint(count = df[counts], 
                              nobs = df[nobs],
                              alpha = alpha,
                              method = method)
    
    var_name_lower = var_name + '_lower'
    var_name_upper = var_name + '_upper'
    df[var_name_lower] = 100 * prop[0]
    df[var_name_upper] = 100 * prop[1]

    df[var_name_lower] = df[var_name_lower].round(round)
    df[var_name_upper] = df[var_name_upper].round(round)
    return df
