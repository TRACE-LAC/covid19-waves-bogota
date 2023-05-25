# -*- coding: utf-8 -*-
"""
Adapted from: https://github.com/mrc-ide/Brazil_COVID19_distributions
@author: davidsantiagoquevedo
@author: ntorresd
"""
import yaml
import pandas as pd
import numpy as np
from scipy import stats, special, integrate
ymlfile = open("config.yml", "r")
cfg = yaml.load(ymlfile)
config = cfg["default"]

DATA_PATH = config['PATHS']['DATA_PATH']
OUT_PATH = config['PATHS']['OUT_PATH'].format(dir = 'epidemiological_distributions')
FIG_PATH = config['PATHS']['FIG_PATH'].format(dir = 'epidemiological_distributions')
SCRIPTS_PATH = config['PATHS']['SCRIPTS_PATH'].format(dir = 'epidemiological_distributions')

SEED = config['MODELS']['SEED']
ITER = config['MODELS']['ITER']
CHAINS = config['MODELS']['CHAINS']
MIN_VAL = config['MODELS']['MIN_VAL']
MAX_VAL = config['MODELS']['MAX_VAL']

######################################################################################
######################################################################################
# Data processing functions 
def prepare_confirmed_cases_data(strat = 'wave'):
    drop_columns = ['Start_date', 'End_date']

    df_icu_stay = pd.read_csv(DATA_PATH + 'icu_stay_bog.csv')
    df_icu_stay = df_icu_stay[(df_icu_stay['icu_stay'] > MIN_VAL) & (df_icu_stay['icu_stay'] <= MAX_VAL)]

    df_hosp_stay = pd.read_csv(DATA_PATH + 'hosp_stay_bog.csv')
    df_hosp_stay = df_hosp_stay[(df_hosp_stay['hosp_stay'] > MIN_VAL) & (df_hosp_stay['hosp_stay'].abs() <= MAX_VAL)]

    df_onset_icu = pd.read_csv(DATA_PATH + 'onset_icu_bog.csv')
    df_onset_icu = df_onset_icu[(df_onset_icu['onset_icu'] > MIN_VAL) & (df_onset_icu['onset_icu'] <= MAX_VAL)]

    df_onset_hosp = pd.read_csv(DATA_PATH + 'onset_hosp_bog.csv')
    df_onset_hosp = df_onset_hosp[(df_onset_hosp['onset_hosp'] > MIN_VAL) & (df_onset_hosp['onset_hosp'] <= MAX_VAL)]

    df_onset_death = pd.read_csv(DATA_PATH + 'onset_death_bog.csv')
    df_onset_death = df_onset_death[(df_onset_death['onset_death'] > MIN_VAL) & (df_onset_death['onset_death'] <= MAX_VAL)]

    all_dfs = [df_icu_stay, df_hosp_stay, df_onset_icu, df_onset_hosp, df_onset_death]

    # clean the data and prepare some the variables list 'columns'
    for df in all_dfs:
        df.dropna(inplace=True)
        
    strat_ages = df_onset_icu['age_group'].unique()
    strat_sex = df_onset_icu['sex'].unique()
    strat_wave= df_onset_icu['wave'].unique()

    strat_sex.sort()
    strat_ages.sort()
    strat_wave.sort()

    strat_sex_map = dict(zip(strat_sex, list(range(1, len(strat_sex)+1))))
    strat_sex = list(range(1, len(strat_sex)+1))

    strat_ages_map = dict(zip(strat_ages, list(range(1, len(strat_ages)+1))))
    strat_ages = list(range(1, len(strat_ages)+1))

    strat_wave_map = dict(zip(strat_wave, list(range(1, len(strat_wave)+1))))
    strat_wave = list(range(1, len(strat_wave)+1))

    if strat == 'wave':
        strat_= strat_wave
    elif strat == 'age':
        strat_= strat_ages
    elif strat == 'sex':
        strat_ = strat_sex

    columns = []
    for df in all_dfs:
        df.dropna(inplace=True) # remove the rows with nan values
        try:
            df.drop(columns = drop_columns, inplace = True)
        except:
            print('No columns to drop')
        col = str(df.columns[4])
        columns.append(col)
        df['age_group_id'] = df['age_group'].map(strat_ages_map)
        df['sex_id'] = df['sex'].map(strat_sex_map)
        df['wave_id'] = df['wave'].astype(int)
    
    return all_dfs, columns

def load_samples(stat = 'mean'):
    samp_posteriors  = {'icu_stay':{},
                        'hosp_stay':{},
                        'onset_icu':{},
                        'onset_hosp':{},
                        'onset_death':{}
                    }

    for col in list(samp_posteriors.keys()):
        samp_posteriors[col].update({'Gamma': pd.read_csv(OUT_PATH + col +'-samples-gamma.csv').agg(stat)})
        samp_posteriors[col].update({'Lognormal': pd.read_csv(OUT_PATH + col +'-samples-log_normal.csv').agg(stat)})
        samp_posteriors[col].update({'Weibull': pd.read_csv(OUT_PATH + col +'-samples-weibull.csv').agg(stat)})
        samp_posteriors[col].update({'Exponential': pd.read_csv(OUT_PATH + col +'-samples-exponential.csv').agg(stat)})
        samp_posteriors[col].update({'Gen Lognormal': pd.read_csv(OUT_PATH + col +'-samples-gln.csv').agg(stat)})
    return samp_posteriors

######################################################################################
######################################################################################
# Fit functions - Bayes inference
def fit_district(values, list_of_params, model):
    stdata = values
    stan_data = {'N': len(stdata), 'y': stdata}
    fit = model.sampling(data = stan_data, iter = ITER, seed = SEED,
                        chains = CHAINS, n_jobs = -1)
    print(fit)                            
    df = fit.to_dataframe()
    df = df[list_of_params]
    return df

def get_posteriors_district(param_list, columns, all_dfs, model):
    district_posteriors = {}
    for i in range(len(columns)):
        df = all_dfs[i]
        col = columns[i]
        print(col)
        vals = df[col].values
        # watch out here!!! we're shifting the data!!!!
        vals = vals + 0.5
        posterior = fit_district(vals, param_list, model)
        district_posteriors.update({col: posterior})  
        
    return district_posteriors

def fit_partial_pooling(df, col, model, params, priors, n_strats, strat_name,):
    stan_pp_data = {'K' : n_strats, 
                    'N' : df.shape[0], 
                    'X' : df[col].values + 0.5,
                    strat_name : df[strat_name+'_id'].values}
    for param in params:
        param_mean = priors[col][param].mean()
        stan_pp_data.update({param + '_prior': param_mean})
    fit = model.sampling(data = stan_pp_data, 
                         iter = ITER,
                         seed = SEED, 
                         chains=CHAINS, n_jobs=-1,
                         control={'adapt_delta': 0.8})
    print(fit)
    posterior_df = fit.to_dataframe()
    params_columns = posterior_df.columns.str.startswith(tuple(params))\
                   + posterior_df.columns.str.startswith('sigma_')
    posterior_df = posterior_df.loc[:,params_columns]
    return posterior_df

def get_posteriors_pooling(all_dfs, columns, model, model_name, param_list, priors, n_strats, strat_name):
    for i in range(len(columns)):
        df = all_dfs[i]
        col = columns[i]
        print(col)
        posteriors_pooling = {}
        posterior = fit_partial_pooling(df, col, model, param_list, priors, n_strats, strat_name)
        # add national estimates
        posterior = pd.concat([posterior, priors[col]], axis = 1, sort = False)
        posteriors_pooling.update({col: posterior})
        # save the output
        posterior.to_csv(OUT_PATH + col + f'-samples-{model_name}.csv', index = False)
        posteriors_pooling.update({col: posterior})

def best_model():
    ep_distributions  = {'icu_stay':{},
                        'hosp_stay':{},
                        'onset_icu':{},
                        'onset_hosp':{},
                        'onset_death':{}
                        }
    best_models = pd.DataFrame({})

    for dist in ep_distributions.keys():
        df = pd.read_csv(OUT_PATH + 'bf_'+dist+'.csv')
        df = df.set_index(df.columns[0])
        ep_distributions[dist].update({'bf' : df,
                                       'best model' : df[ df >= 0].dropna()})
        best_models = pd.concat([best_models, ep_distributions[dist]['best model']])
    cols = ['Epidemilogical distribution'] + best_models.columns.tolist()
    best_models['Epidemilogical distribution'] = list(ep_distributions.keys()) 
    best_models = best_models[cols]
    best_models = best_models.set_index('Epidemilogical distribution')
    return best_models

######################################################################################
######################################################################################
# Statistical functions    
def gln_pdf(x, mu, sigma, g):
    """PDF of the generalised log-normal distribution"""
    k = g / (2**((g+1)/g) * sigma * special.gamma(1/g))
    return k/x * np.exp(-0.5 * np.abs((np.log(x)-mu)/sigma)**g)

def gln_lpdf(x, mu, sigma, g):
    """Log-PDF of the generalised log-normal distribution"""
    logk = np.log(g) - ( ((g+1)/g)*np.log(2) + np.log(sigma) + np.log(special.gamma(1/g)))
    return logk - np.log(x) - 0.5 * np.abs((np.log(x)-mu)/sigma)**g

def gln_cdf_help(x, mu, sigma, g):
    """CDF of the generalised log-normal distribution"""
    m = x
    result = integrate.quad(lambda x: gln_pdf(x,2,0.5,2.5), 0, m)
    tmp = result[0]
    return tmp
    
def gln_cdf(x, mu, sigma, g):
    c = np.vectorize(gln_cdf_help)
    return c(x, mu, sigma, g)

def LogLaplaceCovariance(posterior, col):
    result = 1/2 * len(posterior[col]['mu']) * np.log(2*np.pi) 
    result += 1/2 * np.log(np.linalg.det(posterior[col]['cov']))
    result += posterior[col]['Logf']
    return result

def mean_exponential(beta):
    return beta

def var_exponential(beta):
    return beta**2

def mean_gamma(alpha, beta):
    return alpha/beta

def var_gamma(alpha, beta):
    return alpha/(beta**2)

def mean_weibull(alpha, sigma):
    return sigma * special.gamma(1+(1/alpha))

def var_weibull(alpha, sigma):
    return special.gamma(1+(2/alpha)) - (special.gamma(1+(1/alpha)))**2

def mean_log_normal(mu, sigma):
    return np.exp(mu+0.5*sigma**2)

def var_log_normal(mu, sigma):
    return (np.exp(sigma**2)-1)*np.exp(2*mu + sigma**2)


def sum_term_gln(r, mu, sigma, g, inf):
    sum = 0
    s_k = 0
    for j in range(1, inf): 
        s_k = (r * sigma)**j
        s_k = s_k * (1 + (-1)**j) * 2**(j/g)
        s_k = s_k * (special.gamma((j+1)/g) / special.gamma(j+1))
        s_k = s_k * 1/2 * 1/special.gamma(1/g)
        sum += s_k
    if s_k >=  0.0000001:
        return np.nan
    return sum

def mean_gln(mu, sigma, g, inf = 150):
    return np.exp(mu)*(1 + sum_term_gln(r = 1, mu = mu, sigma = sigma, g = g, inf = inf))

def var_gln(mu, sigma, g, inf = 150):
    return np.exp(2*mu)*(1 + sum_term_gln(r = 2, mu = mu, sigma = sigma, g = g, inf = inf)) - mean_gln(mu, sigma, g, inf = 150)

def q025(x):
    return x.quantile(0.025)

def q975(x):
    return x.quantile(0.975)