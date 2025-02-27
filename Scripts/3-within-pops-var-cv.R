#######################
# 1 - Within-pops variation on CV

## association with the environmental factors 
# Variables:
# Latitude
# Distance to range edge
# Spatial variation
# Temporal variation
#######################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
cv_data_date <- '2025-01-08'
dir <- '~/Documents/MorphCVLatitude/' #change as needed
# dir <- '~/Work/Research/Projects/MorphCVLatitude/'



# load packages -----------------------------------------------------------

library(cmdstanr)
library(MCMCvis)
library(tidyverse)
library(tidyr)
library(shinystan)



# read in data -----------------

cv_data_covs <- readRDS(paste0(dir, 'Data/L2/cv_data_covs-', cv_data_date, '.rds'))



# ------------------- MASS --------------------------------

#filter mass data
cv_data_mass <- subset(cv_data_covs, trait == 'mass')

#negative values get 0
cv_data_mass$distance_edge_km_buffer10 <- ifelse(cv_data_mass$distance_edge_km_buffer10 < 0, 0, cv_data_mass$distance_edge_km_buffer10)
#log distance to range edge
cv_data_mass$ldre <- ifelse(is.na(cv_data_mass$distance_edge_km_buffer10), log(1), log(cv_data_mass$distance_edge_km_buffer10 + 1))

data_mass_str <- cv_data_mass %>% 
  dplyr::group_by(sp_id) %>%
  dplyr::mutate(lat_sc = scale(lat, scale = FALSE), #not scale, just center
                ldre_sc = scale(ldre, scale = FALSE), #not scale, just center
                DHI_spat_sc = scale(log(DHI_spatial_var_10km), scale = FALSE),
                DHI_temp_sc = scale(log(DHI_temporal_var_10km), scale = FALSE)) %>% 
  dplyr::ungroup() %>%
  #fill NA with 0 for one species above
  tidyr::replace_na(list(ldre = 0))


## Create data objects for Stan using PCA
str(DATA_mass <- list(N_obs = NROW(data_mass_str),
                  N_species = length(unique(data_mass_str$sp_id)),
                  CVobs = data_mass_str$cv_within_post_mean*1000,
                  sigmaCV = data_mass_str$cv_within_post_sd*1000,
                  lat = data_mass_str$lat_sc[,1], 
                  distrange = data_mass_str$ldre_sc[,1],
                  spatVar = data_mass_str$DHI_spat_sc[,1], #using only DHI
                  tempVar = data_mass_str$DHI_temp_sc[,1], #using only DHI
                  sp_id = data_mass_str$sp_id))



# Call Stan model -----------------

options(mc.cores = parallel::detectCores())
options("cmdstanr_verbose"=FALSE)

DELTA <- 0.99
#TREE_DEPTH <- 12
#STEP_SIZE <- 0.01
CHAINS <- 4
ITER <- 5000

mass_mod <- cmdstanr::cmdstan_model(paste0(dir, 'Scripts/model_files/3-within-pops-var-cv.stan'),
                                    force_recompile = TRUE)

#sample
fitMass <- mass_mod$sample(
  data = DATA_mass,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 1000,
  adapt_delta = DELTA)#,
  #step_size = STEP_SIZE)



# check results ----------------

shinystan::launch_shinystan(fitMass)

#take a look on results: 
fitMassSummary <- MCMCsummary(fitMass, c('mu_gamma', 'mu_theta1', 'mu_theta2', 'mu_theta3',
                                         'mu_theta4'), pg0=T)



# save summary -----------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fitMass,
                  round = 4,
                  file_name = paste0('fitMass-within-pops-CV-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('fitMass-within-pops-CV-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('fitMass-fit-', run_date),
                  add_obj = list(DATA_mass, data_mass_str),
                  add_obj_names = c(paste0('fitMass-data-', run_date), paste0('raw-mass-data-', run_date)),
                  cp_file = c('Scripts/model_files/3-within-pops-var-cv.stan',
                              'Scripts/3-within-pops-var-cv.R'),
                  cp_file_names = c(paste0('3-within-pops-var-cv', run_date, '.stan'),
                                    paste0('3-within-pops-var-cv', run_date, '.R')))



# Plot all effects --------------------------------------------

param_labels <- c('Latitude', 
                  'Distance range edge', 
                  'Spatial variation of productivity', 
                  'Temporal variation of productivity')
# Plot adding labels
Mass <- MCMCplot(fitMass, 
                 params = c('mu_theta1', 'mu_theta2', 'mu_theta3', 'mu_theta4'), 
                 main = 'CV body mass',
                 ci = c(50, 89),
                 labels = param_labels, sz_labels=1.05, sz_med =1.1)




# ------------------ WING ---------------------------------

#filter wing data
cv_data_wing <- subset(cv_data_covs, trait == "wing")

#negative values get 0
cv_data_wing$distance_edge_km_buffer10 <- ifelse(cv_data_wing$distance_edge_km_buffer10 < 0, 0, cv_data_wing$distance_edge_km_buffer10)
#log distance to range edge
cv_data_wing$ldre <- ifelse(is.na(cv_data_wing$distance_edge_km_buffer10), log(1), log(cv_data_wing$distance_edge_km_buffer10 + 1))

data_wing_str <- cv_data_wing %>% 
  dplyr::group_by(sp_id) %>%
  dplyr::mutate(lat_sc = scale(lat, scale = FALSE),
                ldre_sc = scale(ldre, scale = FALSE),
                DHI_spat_sc = scale(log(DHI_spatial_var_10km), scale = FALSE),
                DHI_temp_sc = scale(log(DHI_temporal_var_10km), scale = FALSE)) %>%
  dplyr::ungroup() %>%
  #fill NA with 0 for one species above
  tidyr::replace_na(list(ldre = 0))


## Create data objects for Stan
str(DATA_wing <- list(N_obs = NROW(data_wing_str),
                  N_species = length(unique(data_wing_str$sp_id)),
                  CVobs = data_wing_str$cv_within_post_mean*1000,
                  sigmaCV = data_wing_str$cv_within_post_sd*1000,
                  lat = data_wing_str$lat_sc[,1], 
                  distrange = data_wing_str$ldre_sc[,1],
                  spatVar = data_wing_str$DHI_spat_sc[,1], #using only DHI
                  tempVar = data_wing_str$DHI_temp_sc[,1], #using only DHI
                  sp_id = data_wing_str$sp_id))



# Call Stan model -----------------

options(mc.cores = parallel::detectCores())
options("cmdstanr_verbose"=FALSE)

DELTA <- 0.99
#TREE_DEPTH <- 12
#STEP_SIZE <- 0.0005
CHAINS <- 4
ITER <- 5000

wing_mod <- cmdstanr::cmdstan_model(paste0(dir, 'Scripts/model_files/3-within-pops-var-cv.stan'),
                                    force_recompile = TRUE)

#sample
fitWing <- wing_mod$sample(
  data = DATA_wing,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 1000,
  adapt_delta = DELTA)#, 
  #step_size = STEP_SIZE)#,
  #max_treedepth = TREE_DEPTH)



# check results -----------------

shinystan::launch_shinystan(fitWing)

#take a look on results: 
fitWingSummary <- MCMCsummary(fitWing, c('mu_gamma', 'mu_theta1', 'mu_theta2', 'mu_theta3',
                                         'mu_theta4'), pg0=T)



# save summary -----------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fitWing,
                  round = 4,
                  file_name = paste0('fitWing-within-pops-CV-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('fitWing-within-pops-CV-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('fitWing-fit-', run_date),
                  add_obj = list(DATA_wing, data_wing_str),
                  add_obj_names = c(paste0('fitWing-data-', run_date), paste0('raw-wing-data-', run_date)),
                  cp_file = c('Scripts/model_files/3-within-pops-var-cv.stan',
                              'Scripts/3-within-pops-var-cv.R'),
                  cp_file_names = c(paste0('3-within-pops-var-cv', run_date, '.stan'),
                                    paste0('3-within-pops-var-cv', run_date, '.R')))



# Plot all effects --------------------------------------------

param_labels <- c('Latitude', 
                  'Distance range edge', 
                  'Spatial variation of productivity', 
                  'Temporal variation of productivity')
MCMCvis::MCMCplot(fitWing,
         params = c('mu_theta1', 'mu_theta2', 'mu_theta3', 'mu_theta4'), 
                 main = 'CV wing length',
                 ci = c(50, 89),
                 labels = param_labels, sz_labels=1.05, sz_med =1.1)


