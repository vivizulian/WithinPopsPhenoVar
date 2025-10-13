#######################
# 2 - Within-pops variation on CV among species

#association with the environmental factors 
#code to estimate the effect of range size, generation length, HWI, 
#and Migratory status on within-population CV across species. 
#code also account for residuals with phylogenetic signal. 
#######################


# load packages -----------------------------------------------------------

library(cmdstanr)
library(MCMCvis)
library(tidyverse)
library(shinystan)
library(ape)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
cv_data_date <- '2025-01-08'
tree_date <- '2025-01-10'
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------

cv_data_covs <- readRDS(paste0(dir, 'Data/L2/cv_data_covs-', cv_data_date, '.rds'))

tree <- readRDS(paste0(dir, 'Data/L1/consensus_tree-', tree_date, '.rds'))


#get corr matrix of tree
Rho <- ape::vcv.phylo(tree, corr = TRUE)

#matrix: 
#Diagonal elements represent the total branch length from the root to each tip (variance).
#Off-diagonal elements are covariances for tip taxa.


# ------------------ MASS ---------------------------------

#filter mass data
cv_data_mass <- subset(cv_data_covs, trait == "mass")

cv_data_mass <- cv_data_mass  %>% 
  dplyr::distinct(sp_id, .keep_all = TRUE) %>%
  dplyr::mutate(sci_name = gsub(' ', '_', sci_name))

## Need to make sure the order of the species is the same as the mass dataset
# order of species from mass data set
cv_sp_order <- cv_data_mass$sci_name
cv_sp_order <- cv_sp_order[cv_sp_order %in% rownames(Rho)]

# Reorder the 'V' matrix based on the species order in 'cv_sp_order'
Rho <- Rho[cv_sp_order, cv_sp_order]

cv_data_mass$migstatus <- ifelse(as.numeric(cv_data_mass$migration_dist_km)==0, 0, 1)
# 0 = species that do not migrate (n=11)
# 1 = species that migrate (n=88)


#log transform skewed data
cv_data_mass$range_size_log <- log(as.numeric(cv_data_mass$range_size_km2))
cv_data_mass$GenLength_log <- log(as.numeric(cv_data_mass$GenLength))
cv_data_mass$HWI_log <- log(as.numeric(cv_data_mass$HWI))

data_mass_str <- cv_data_mass %>% 
  dplyr::mutate(range_size_str = scale(range_size_log, scale = TRUE),
                GenLength_str = scale(GenLength_log, scale = TRUE),
                HWI_str = scale(HWI_log, scale = TRUE))


## Create data objects for Stan
str(DATA_mass <- list(N_sp = length(unique(data_mass_str$sp_id)),
                      CVobs = data_mass_str$cv_within_among_sp_post_mean*1000,
                      sigmaCV = data_mass_str$cv_within_among_sp_post_sd*1000,
                      rangesize = data_mass_str$range_size_str[,1],
                      genlength = data_mass_str$GenLength_str[,1],
                      HWI = data_mass_str$HWI_str[,1],
                      migstatus = data_mass_str$migstatus,
                      Rho = Rho))


# Call Stan model -----------------

options(mc.cores = parallel::detectCores())
options("cmdstanr_verbose"=FALSE)

DELTA <- 0.99
#TREE_DEPTH <- 12
#STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 5000

mass_mod <- cmdstanr::cmdstan_model(paste0('Scripts/model_files/4-within-pops-among-sp.stan'),
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
  #max_treedepth = TREE_DEPTH)



# check results ----------------

shinystan::launch_shinystan(fitMass)

fitMassSummary <- MCMCsummary(fitMass, c('mu_gamma', 'theta1', 'theta2', 'theta3',
                                         'theta4'), pg0=T)



# save summary -----------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fitMass,
                  round = 4,
                  file_name = paste0('fitMassCV-within-among-sp-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('fitMassCV-within-among-sp-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('fitMass-fit-', run_date),
                  add_obj = list(DATA_mass, data_mass_str),
                  add_obj_names = c(paste0('fitMass-data-', run_date), paste0('raw-mass-data-', run_date)),
                  cp_file = c('Scripts/model_files/4-within-pops-among-sp.stan',
                              'Scripts/4-within-pops-among-sp-var-cv.R'),
                  cp_file_names = c(paste0('4-within-pops-among-sp-', run_date, '.stan'),
                                    paste0('4-within-pops-among-sp-', run_date, '.R')))



# plot results -----------------

# Plot all effects:
param_labels <- c("Range Size", 
                  "Generation Length", 
                  "Hand Wing Index", 
                  "Migratory status = 1")

#pdf(paste0(dir, "Results/CV_sp.pdf"), width = 12, height = 8)

# Plot adding labels
Mass <- MCMCplot(fitMass, 
         params = c("theta1", "theta2", "theta3", "theta4"), 
         main = "CV body mass",
         ci = c(50, 89),
         labels = param_labels, sz_labels=1.05, sz_med =1.1)#, sz_ax =1, 
         #sz_ax_txt=1, sz_tick_txt=1, sz_main_txt=1)

#dev.off()





# ------------------ WING ---------------------------------

#filter wing data
cv_data_wing <- subset(cv_data_covs, trait == "wing")

cv_data_wing <- cv_data_wing  %>% 
  dplyr::distinct(sp_id, .keep_all = TRUE) %>%
  dplyr::mutate(sci_name = gsub(' ', '_', sci_name))

## Need to make sure the order of the species is the same as the mass dataset
# order of species from mass data set
cv_sp_order <- cv_data_wing$sci_name
cv_sp_order <- cv_sp_order[cv_sp_order %in% rownames(Rho)]

# Reorder the 'V' matrix based on the species order in 'cv_sp_order'
Rho <- Rho[cv_sp_order, cv_sp_order]


cv_data_wing$migstatus <- ifelse(as.numeric(cv_data_wing$migration_dist_km)==0, 0, 1)
# 0 = species that do not migrate (n=11)
# 1 = species that migrate (n=88)


#log transform skewed data
cv_data_wing$range_size_log <- log(as.numeric(cv_data_wing$range_size_km2))
cv_data_wing$GenLength_log <- log(as.numeric(cv_data_wing$GenLength))
cv_data_wing$HWI_log <- log(as.numeric(cv_data_wing$HWI))

data_wing_str <- cv_data_wing %>% 
  dplyr::mutate(range_size_str = scale(range_size_log, scale = TRUE),
                GenLength_str = scale(GenLength_log, scale = TRUE),
                HWI_str = scale(HWI_log, scale = TRUE))


## Create data objects for Stan
str(DATA_wing <- list(N_sp = length(unique(data_wing_str$sp_id)),
                      CVobs = data_wing_str$cv_within_among_sp_post_mean*1000,
                      sigmaCV = data_wing_str$cv_within_among_sp_post_sd*1000,
                      rangesize = data_wing_str$range_size_str[,1],
                      genlength = data_wing_str$GenLength_str[,1],
                      migstatus = data_wing_str$migstatus,
                      HWI = data_wing_str$HWI_str[,1],
                      Rho = Rho))


# Call Stan model -----------------

options(mc.cores = parallel::detectCores())
options("cmdstanr_verbose"=FALSE)

DELTA <- 0.99
#TREE_DEPTH <- 12
#STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 5000

wing_mod <- cmdstanr::cmdstan_model(paste0(dir, 'Scripts/model_files/4-within-pops-among-sp.stan'),
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



# check results ----------------

shinystan::launch_shinystan(fitWing)

fitWingSummary <- MCMCsummary(fitWing, c("mu_gamma", "theta1", "theta2", "theta3",
                                         "theta4"), pg0=T)


 
# save summary -----------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fitWing,
                  round = 4,
                  file_name = paste0('fitWingCV-within-among-sp-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('fitWingCV-within-among-sp-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('fitWing-fit-', run_date),
                  add_obj = list(DATA_wing, data_wing_str),
                  add_obj_names = c(paste0('fitWing-data-', run_date), paste0('raw-wing-data-', run_date)),
                  cp_file = c('Scripts/model_files/4-within-pops-among-sp.stan',
                              'Scripts/4-within-pops-among-sp-var-cv.R'),
                  cp_file_names = c(paste0('4-within-pops-among-sp-', run_date, '.stan'),
                                    paste0('4-within-pops-among-sp-', run_date, '.R')))



# plot results -----------------

# Plot all effects:
param_labels <- c("Range Size", 
                  "Generation Length", 
                  "Hand Wing Index", 
                  "Migratory status = 1")

#pdf(paste0(dir, "Results/CVper_spWing.pdf"), width = 8, height = 10)

# Plot adding labels
Wing <- MCMCplot(fitWing, 
         params = c("theta1", "theta2", "theta3", "theta4"), 
         main = "CV wing length", 
         labels = param_labels, sz_labels=1.05, sz_med =1.1)#, sz_ax =1, 
         #sz_ax_txt=1, sz_tick_txt=1, sz_main_txt=1)

#dev.off()



