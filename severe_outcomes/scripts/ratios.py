# -*- coding: utf-8 -*-
"""
@author: davidsantiagoquevedo
@author: ntorresd
"""
import sys
import yaml
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

config = yaml.load(open("config.yml", "r"))["default"]
DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'severe_outcomes')
SCRIPTS_PATH = config['PATHS']['SCRIPTS_PATH'].format(dir = 'severe_outcomes')
UTILS_PATH = config['PATHS']['UTILS_PATH'].format(dir = 'severe_outcomes')
UPDATE = config['UPDATE_DATES']['CONFIRMED_CASES']

# Import useful functions
sys.path.append(UTILS_PATH)
import utilities_severity as ut
df_confirmed_bogota = pd.read_csv(DATA_PATH + f'confirmed_cases_waves_{UPDATE}.csv')

# Change age groups
age_group_dic={'0-9':[0,9],
              '10-19':[10,19],
              '20-29':[20,29], 
              '30-39':[30,39], 
              '40-49':[40,49], 
              '50-59':[50,59],
              '60-69':[60,69],
              '70-79':[70,79], 
              '80+':[80,9999]}

df_confirmed_bogota = ut.age_group_dec(df_confirmed_bogota, 'age', 'age_unit', age_group_dic)

# Counts of cases and outcomes by wave

#### Cases
cases_wave = df_confirmed_bogota[['wave', 'age_group']].value_counts()\
            .reset_index().rename(columns = {0 : 'cases'})

cases_all = df_confirmed_bogota[['wave']].value_counts()\
            .reset_index().rename(columns = {0 : 'cases'})
cases_all['age_group'] = 'all'

cases_wave = pd.concat([cases_wave, cases_all])


#### Hospitalisations
hosp = \
df_confirmed_bogota[(df_confirmed_bogota.hospitalization.notnull())]\
                    [['wave', 'age_group']].value_counts()\
                    .reset_index().rename(columns = {0 : 'hosp'})
                    
hosp_all = df_confirmed_bogota[(df_confirmed_bogota.hospitalization.notnull())]\
                    [['wave']].value_counts()\
                    .reset_index().rename(columns = {0 : 'hosp'})
hosp_all['age_group'] = 'all'

hosp = pd.concat([hosp, hosp_all])

#### ICUs                    
icu = \
df_confirmed_bogota[(df_confirmed_bogota.icu.notnull())]\
                    [['wave', 'age_group']].value_counts()\
                    .reset_index().rename(columns = {0 : 'icu'})

icu_all = df_confirmed_bogota[(df_confirmed_bogota.icu.notnull())]\
                    [['wave']].value_counts()\
                    .reset_index().rename(columns = {0 : 'icu'})
icu_all['age_group'] = 'all'


icu = pd.concat([icu, icu_all])

#### HOSP - DEATH
hosp_death = \
df_confirmed_bogota[(df_confirmed_bogota.hospitalization.notnull()) &
                    (df_confirmed_bogota.condition == 'Fallecido') & 
                    (df_confirmed_bogota.death.notnull())]\
                    [['wave', 'age_group']].value_counts()\
                    .reset_index().rename(columns = {0 : 'hosp_death'})

h_d_all = df_confirmed_bogota[(df_confirmed_bogota.hospitalization.notnull())&
                                (df_confirmed_bogota.condition == 'Fallecido') & 
                                (df_confirmed_bogota.death.notnull())]\
                    [['wave']].value_counts()\
                    .reset_index().rename(columns = {0 : 'hosp_death'})
h_d_all['age_group'] = 'all'

hosp_death = pd.concat([hosp_death, h_d_all])

#### ICU - DEATH                    
icu_death = \
df_confirmed_bogota[(df_confirmed_bogota.icu.notnull()) &
                    (df_confirmed_bogota.condition == 'Fallecido') & 
                    (df_confirmed_bogota.death.notnull())]\
                    [['wave', 'age_group']].value_counts()\
                    .reset_index().rename(columns = {0 : 'icu_death'})

i_d_all = df_confirmed_bogota[(df_confirmed_bogota.icu.notnull())&
                                (df_confirmed_bogota.condition == 'Fallecido') & 
                                (df_confirmed_bogota.death.notnull())]\
                    [['wave']].value_counts()\
                    .reset_index().rename(columns = {0 : 'icu_death'})
i_d_all['age_group'] = 'all'


icu_death = pd.concat([icu_death, i_d_all])

#### DEATHS                    
deaths = \
df_confirmed_bogota[(df_confirmed_bogota.condition == 'Fallecido') & 
                    (df_confirmed_bogota.death.notnull())]\
                    [['wave', 'age_group']].value_counts()\
                    .reset_index().rename(columns = {0 : 'deaths'})

death_all = df_confirmed_bogota[(df_confirmed_bogota.condition == 'Fallecido') & 
                                (df_confirmed_bogota.death.notnull())]\
                    [['wave']].value_counts()\
                    .reset_index().rename(columns = {0 : 'deaths'})
death_all['age_group'] = 'all'

deaths = pd.concat([deaths, death_all])

#### Merge all the outcomes and cases
cases_wave = \
cases_wave.merge(hosp, how = 'left', on = ['wave', 'age_group'])\
                .merge(icu, how = 'left', on = ['wave', 'age_group'])\
                .merge(hosp_death, how = 'left', on = ['wave', 'age_group'])\
                .merge(icu_death, how = 'left', on = ['wave', 'age_group'])\
                .merge(deaths, how = 'left', on = ['wave', 'age_group'])

cat_type = pd.api.types.CategoricalDtype(categories = ['all', '60+']+list(age_group_dic.keys()),
                                         ordered=True)
cases_wave['age_group'] = cases_wave['age_group'].astype(cat_type)
cases_wave = cases_wave.sort_values(by = ['wave', 'age_group'])

# Ratios and confidence intervals using binomial proportions
alpha = 0.05 #significance level
method = 'normal'

#CFR
cases_wave['CFR'] =  (cases_wave.deaths / cases_wave.cases).round(4)
prop = proportion_confint(count = cases_wave.deaths,
                         nobs = cases_wave.cases,
                         alpha = alpha,
                         method = method)
cases_wave['CFR_lower'] = prop[0].round(4)
cases_wave['CFR_upper'] = prop[1].round(4)

#HCR
cases_wave['HCR'] =  (cases_wave.hosp / cases_wave.cases).round(4)
prop = proportion_confint(count = cases_wave.hosp,
                         nobs = cases_wave.cases,
                         alpha = alpha,
                         method = method)
cases_wave['HCR_lower'] = prop[0].round(4)
cases_wave['HCR_upper'] = prop[1].round(4)

#ICU-CR
cases_wave['ICU-CR'] =  (cases_wave.icu / cases_wave.cases).round(4)
prop = proportion_confint(count = cases_wave.icu,
                         nobs = cases_wave.cases,
                         alpha = alpha,
                         method = method)
cases_wave['ICU-CR_lower'] = prop[0].round(4)
cases_wave['ICU-CR_upper'] = prop[1].round(4)

#HFR
cases_wave['HFR'] =  (cases_wave.hosp_death / cases_wave.hosp).round(4)
prop = proportion_confint(count = cases_wave.hosp_death,
                         nobs = cases_wave.hosp,
                         alpha = alpha,
                         method = method)
cases_wave['HFR_lower'] = prop[0].round(4)
cases_wave['HFR_upper'] = prop[1].round(4)

#ICU-FR
cases_wave['ICU-FR'] =  (cases_wave.icu_death / cases_wave.icu).round(4)
prop = proportion_confint(count = cases_wave.icu_death,
                         nobs = cases_wave.icu,
                         alpha = alpha,
                         method = method)
cases_wave['ICU-FR_lower'] = prop[0].round(4)
cases_wave['ICU-FR_upper'] = prop[1].round(4)


rates = ['CFR', 'HCR', 'ICU-CR', 'HFR', 'ICU-FR']
for col in rates:
    cases_wave[col] = cases_wave[col].fillna(0)
    cases_wave[col+'_lower'] = cases_wave[col+'_lower'].fillna(0)
    cases_wave[col+'_upper'] = cases_wave[col+'_upper'].fillna(0)

cases_wave.to_csv(OUT_PATH + 'ratios.csv', index = False)