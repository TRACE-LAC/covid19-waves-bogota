rm(list = ls())

library(rstudioapi)
library(EpiEstim)
library(dplyr)
library(tidyverse)

library(config)
library(here)

# set the working directory
setwd(here())

# paths
config <- config::get('PATHS')
DATA_PATH = config$DATA_PATH
OUT_PATH = config$OUT_PATH %>% str_replace('\\{dir\\}', 'rt')

config <- config::get('UPDATE_DATES')
UPDATE = config$CONFIRMED_CASES

# parameters 
incubation_period <- 5
rt_window <- 14

# read data
df_confirmed_bogota <- read_csv(paste0(DATA_PATH, 'confirmed_cases_', UPDATE,'.csv'))
df_confirmed_bogota <- df_confirmed_bogota[, c('onset', 'age', 'age_unit', 'contagion_source')]
df_confirmed_bogota$onset <- as.Date(df_confirmed_bogota$onset)
df_confirmed_bogota$infection <- df_confirmed_bogota$onset - incubation_period

# define age groups
df_confirmed_bogota <- df_confirmed_bogota %>% 
  mutate(
    age_group = case_when(
      age_unit == 1 & age < 60 ~ "<60",
      age_unit == 2 ~ "<60",
      age_unit == 3 ~ "<60",
      age_unit == 1 & age >= 60 ~ "60+"
    ),
    age_group = factor(age_group, level=c("<60", "60+"))
  )

# compute incidence
get_incidence <- function(df) {
  df_local <- df %>%
    filter(contagion_source == "local") %>%
    group_by(infection) %>%
    summarise(local = n())


  df_imported <- df %>%
    filter(contagion_source == "imported") %>%
    group_by(infection) %>%
    summarise(imported = n())

  df_incidence <- df_local %>% full_join(df_imported, by = "infection")
  return(df_incidence)
}

# complete missing dates
complete_dates <- function(df) {
  all_dates <- data.frame(infection = seq.Date(min(df$infection),
                                               max(df$infection),
                                               by = "day")
  )

  df <- left_join(all_dates, df)
  df <- df %>% replace(is.na(df), 0)
  return(df)
}

df_incidence <- get_incidence(df_confirmed_bogota) %>% complete_dates()
mask = (df_confirmed_bogota$age_group == "60+")
df_incidence_60p <- get_incidence(df_confirmed_bogota[mask, ]) %>% complete_dates()

# compute Rt
compute_rt <- function(df_incidence,
                       method = "parametric_si",
                       mean_si = 6.48,
                       std_si = 3.83,
                       imported_ = TRUE) {
  t_start <- seq(incubation_period, nrow(df_incidence) - rt_window)
  t_end <- t_start + rt_window
  if(imported_){
    local <- df_incidence
  }
  rt_data <- estimate_R(df_incidence,
                        method = method,
                        config = make_config(list(
                          mean_si = mean_si,
                          std_si = std_si,
                          t_start = t_start,
                          t_end = t_end)))
  
  df_rt <- rt_data$R
  df_rt$window_start <- min(df_incidence$infection) + df_rt$t_start
  df_rt$window_end <- min(df_incidence$infection) + df_rt$t_end

  return(df_rt)
}
df_rt <- df_incidence %>% compute_rt()
df_rt_60p <- df_incidence_60p %>% compute_rt()

# save data
write_csv(df_rt, paste0(OUT_PATH, "rt_all_ages.csv"))
write_csv(df_rt_60p, paste0(OUT_PATH, "rt_60_plus.csv"))