################
## Code to derive metrics related to trait measurement and cv in different levels of biological organization from
## NIMBLE model output

# Mean and uncertainty for trait measurement and cv for within-pops for station/species and for within-pops among species
# Mean and uncertainty for trait measurement and cv for among species (within each location) 
# Mean and uncertainty for trait measurement and cv for 

################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
model_run_date <- "2024-12-12"
dir <- '~/Documents/MorphCVLatitude/' #change as needed


# load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(tidyr)


# read in data -------------------------------------------------------------

# MASS results--------------------

fitMass <- readRDS(paste0(dir, 'Results/MassCV-output-', model_run_date, '/MassCV-fit-', model_run_date, '.rds')) 
massData <- readRDS(paste0(dir, 'Results/MassCV-output-', model_run_date, '/MassCV-data-', model_run_date, '.rds'))


# WING results--------------------

fitWing <- readRDS(paste0(dir, 'Results/WingCV-output-', model_run_date, '/WingCV-fit-', model_run_date, '.rds')) 
wingData <- readRDS(paste0(dir, 'Results/WingCV-output-', model_run_date, '/WingCV-data-', model_run_date, '.rds'))



# derive metrics -------------------------------------------------------------

# function to extract all values of interest from model files

cv_fun <- function(fit_object, data_object, trait_name){
  
  # WITHIN-POPS
  
  # 1) Measurements (trait values)
  
  # 1.a) Extract the posterior estimates for mean_measurement within-pops:
  
  mean_meas <- MCMCvis::MCMCchains(fit_object, params = 'mean_measurement') #extract the chains for each species/station - 'mean_measurement' (mean_within_pops)
  mean_meas_within_pops_post_mean <- apply(mean_meas, 2, FUN = mean)
  mean_meas_within_pops_post_sd <- apply(mean_meas, 2, FUN = sd)
  
  
  # 1.b) extract estimated sd measurement: 'sd_measurement' (sd_within_pops) 
  
  sd_meas_within_pops <- MCMCvis::MCMCchains(fit_object, params = 'sd_measurement')
  sd_meas_within_pops_post_mean <- apply(sd_meas_within_pops, 2, FUN = mean)
  sd_meas_within_pops_post_sd <- apply(sd_meas_within_pops, 2, FUN = sd)
  
  
  
  # 2) CV
  
  # 2.a) Extract the posterior estimates for CV within-pops:
  
  cv_within <- MCMCvis::MCMCchains(fit_object, params = 'cv') #extract the chains for each species/station, mean, and sd:
  cv_within_pops_post_mean <- apply(cv_within, 2, FUN = mean)
  cv_within_pops_post_sd <- apply(cv_within, 2, FUN = sd)
  
  
  # 2.b) Join results in a data.frame
  
  within_pops <- cbind.data.frame(cn_id = as.numeric(unique(data_object$cn_id)), 
                                  sp_id = as.numeric(data_object$sp_id),
                                  station_id = as.numeric(data_object$station_id),
                                  trait = rep(trait_name, length(unique(data_object$cn_id))),
                                  mean_meas_within_post_mean = mean_meas_within_pops_post_mean, 
                                  mean_meas_within_post_sd = mean_meas_within_pops_post_sd,
                                  sd_meas_within_post_mean = sd_meas_within_pops_post_mean, 
                                  sd_meas_within_post_sd = sd_meas_within_pops_post_sd,
                                  cv_within_post_mean = cv_within_pops_post_mean, 
                                  cv_within_post_sd = cv_within_pops_post_sd)
  
  
  # 2.c) CV within-pops for species (among sp)
  
  cv_within_among_sp <- MCMCvis::MCMCchains(fit_object, params = 'cv')
  row.names(cv_within_among_sp) <- seq(1:nrow(cv_within_among_sp))
  
  #transpose so each line is a species/station and columns are the iterations:
  cv_within_among_sp_t <- cbind.data.frame(t(cv_within_among_sp), sp_id = data_object$sp_id)
  
  
  # 3) Join all within pops results and export as .csv
  
  cv_within <- cv_within_among_sp_t %>%
    tidyr::pivot_longer(-sp_id, names_to = "variable", values_to = "value") %>%
    dplyr::group_by(sp_id) %>%
    dplyr::summarise(cv_within_among_sp_post_mean = mean(value), 
                     cv_within_among_sp_post_sd = sd(value)) %>%
    dplyr::left_join(within_pops, by = "sp_id") %>%
    dplyr::relocate(c(cv_within_among_sp_post_mean, cv_within_among_sp_post_sd), .after = cv_within_post_sd) 
  
  
  # Export results in a .cvs
  
  write.csv(cv_within, file=paste0(dir, 'Data/L2/cv_within_pops_', trait_name, run_date,'.csv'))


  
  # ACROSS-POPS
  
  # 1) Measurements (trait values)
  
  # 1.a) Extract the posterior estimates for the mean (mu_k) and uncertainty (sigma_k) of mean measurement across-pops:
  
  #mu_k
  
  mu_mean_measure_sp <- MCMCvis::MCMCchains(fit_object, params = 'mu') #estimated species-specific mean measurement for wing length or mass (mu_k).
  mu_meas_across_pops_post_mean <- apply(mu_mean_measure_sp, 2, FUN = mean) #mu
  mu_meas_across_pops_post_sd <- apply(mu_mean_measure_sp, 2, FUN = sd) #sd mu
  
  
  #sigma_k
  
  sd_mean_meas_across_pops <- MCMCvis::MCMCchains(fit_object, params = 'sigma')
  sigma_meas_across_pops_post_mean <- apply(sd_mean_meas_across_pops, 2, FUN = mean) #mean sigma
  sigma_meas_across_pops_post_sd <- apply(sd_mean_meas_across_pops, 2, FUN = sd) #sd sigma
  
  
  
  # 1.b) Extract the posterior estimates for the uncertainty of mean measurement across-pops:
  
  #mu_sd
  
  mean_uncert_sd_measure <- MCMCvis::MCMCchains(fit_object, params = 'mu_sd')
  mu_sd_meas_across_pops_post_mean <- apply(mean_uncert_sd_measure, 2, FUN = mean) #mean mu_sd
  mu_sd_meas_across_pops_post_sd <- apply(mean_uncert_sd_measure, 2, FUN = sd) #sd mu_sd
  
  #sigma_sd
  
  sd_uncert_sd_measure <- MCMCvis::MCMCchains(fit_object, params = 'sigma_sd')
  sigma_sd_meas_across_pops_post_mean <- apply(sd_uncert_sd_measure, 2, FUN = mean) #mean sigma_sd
  sigma_sd_meas_across_pops_post_sd <- apply(sd_uncert_sd_measure, 2, FUN = sd) #sd sigma_sd
  
  
  
  # 2) CV
  
  # 2.a) Extract the posterior estimates, mean, and sd for cv across-pops:
  
  cv_across_pops <- MCMCvis::MCMCchains(fit_object, params = 'cv_acrosspops')
  cv_across_pops_post_mean <- apply(cv_across_pops, 2, FUN = mean)
  cv_across_pops_post_sd <- apply(cv_across_pops, 2, FUN = sd)
  
  
  
  # 3) Join all within pops results and export as .csv
  
  cv_across_pops <- cbind.data.frame(sp_id = unique(data_object$sp_id),
                                          mu_meas_across_post_mean = mu_meas_across_pops_post_mean, 
                                          mu_meas_across_post_sd = mu_meas_across_pops_post_sd,
                                          sigma_meas_across_post_mean = sigma_meas_across_pops_post_mean, 
                                          sigma_meas_across_post_sd = sigma_meas_across_pops_post_sd,
                                          mu_sd_meas_across_post_mean = mu_sd_meas_across_pops_post_mean, 
                                          mu_sd_meas_across_post_sd = mu_sd_meas_across_pops_post_sd,
                                          sigma_sd_meas_across_post_mean = sigma_sd_meas_across_pops_post_mean, 
                                          sigma_sd_meas_across_post_sd = sigma_sd_meas_across_pops_post_sd,
                                          cv_across_post_mean = cv_across_pops_post_mean, 
                                          cv_across_post_sd = cv_across_pops_post_sd)
  
  rownames(cv_across_pops) <- NULL #remove row names
  
  
  # Export results in a .cvs
  
  write.csv(cv_across_pops, file=paste0(dir, 'Data/L2/cv_across_pops_', trait_name, run_date,'.csv'))
  
  
  
  # AMONG SPECIES
  
  # 1) Measurements (trait values) and CV together
  
  # 1.a) Extract the posterior estimates for the mean and uncertainty of mean measurement among-sp:
  
  mean_measurement_ch <- MCMCvis::MCMCchains(fit_object, params = 'mean_measurement')
  row.names(mean_measurement_ch) <- seq(1:nrow(mean_measurement_ch)) #substitute the row names by just a sequence from 1:number of iterations
  mean_measurement_among <- cbind.data.frame(t(mean_measurement_ch), station_id = data_object$station_id) #transpose so each line is a different sp_id and columns are the iterations:
  
  
  # 1.b) group by station_id and iteration and calculate the mean and sd measurement and cv:
  
  data_among_sp_long <- mean_measurement_among %>%
    pivot_longer(-station_id, names_to = "variable", values_to = "value") %>%
    group_by(variable, station_id) %>%
    summarise(mean_meas_among_sp = mean(value, na.rm = TRUE), #here summarizing mean and sd measur per iteration and station
              sd_meas_among_sp = sd(value, na.rm = TRUE), #stations with only one sp have NA on sd - should be 160 stations.
              cv_among_sp = sd_meas_among_sp/mean_meas_among_sp) 
  
  
  # 1.c) extract the mean and uncertainty of measurement and cv over all the iterations for each station:
  
  cv_among_sp_summ <- data_among_sp_long %>%
    group_by(station_id) %>%
    summarise(mean_meas_among_sp_post_mean = mean(mean_meas_among_sp, na.rm = TRUE),
              mean_meas_among_sp_post_sd = sd(mean_meas_among_sp, na.rm = TRUE),
              sd_meas_among_sp_post_mean = mean(sd_meas_among_sp, na.rm = TRUE),
              sd_meas_among_sp_post_sd = sd(sd_meas_among_sp, na.rm = TRUE),
              cv_among_sp_post_mean = mean(cv_among_sp, na.rm = TRUE), 
              cv_among_sp_post_sd = sd(cv_among_sp, na.rm = TRUE))
  
  
  # Export results in a .cvs
  
  write.csv(cv_among_sp_summ, file=paste0(dir, 'Data/L2/cv_among_sp_', trait_name, run_date,'.csv'))
  
  
  
  # COMMUNITY
  
  # 1) Measurements (trait values) and CV together
  
  # 1.a) Extract the posterior estimates for the mean and uncertainty of mean measurement:
  
  mean_measurement_ch <- MCMCvis::MCMCchains(fit_object, params = 'mean_measurement')
  row.names(mean_measurement_ch) <- seq(1:nrow(mean_measurement_ch)) #substitute the row names by just a sequence from 1:number of iterations
  mean_measurement_commun <- cbind.data.frame(t(mean_measurement_ch), station_id = data_object$station_id) #transpose so each line is a different sp_id and columns are the iterations:
  
  
  # 1.b) group by station_id and iteration and calculate the mean and sd measurement, and cv:
  
  mean_commun_ch <- mean_measurement_commun %>%
    pivot_longer(-station_id, names_to = "variable", values_to = "value") %>%
    group_by(variable, station_id) %>%
    summarise(mean_community = mean(value, na.rm = TRUE), #here summarizing mean and sd measur per iteration and stations
              sd_community = sd(value, na.rm = TRUE),
              cv_community = sd_community/mean_community)  %>% 
    ungroup() 
  
  
  # 1.c) extract the mean and uncertainty of measurement and cv over all the stations
  
  cv_community_summ <- mean_commun_ch %>%
    dplyr::summarise(mean_meas_comm_post_mean = mean(mean_community, na.rm = TRUE),
                     mean_meas_comm_post_sd = sd(mean_community, na.rm = TRUE),
                     sd_meas_comm_post_mean = mean(sd_community, na.rm = TRUE),
                     sd_meas_comm_post_sd = sd(sd_community, na.rm = TRUE),
                     cv_comm_post_mean = mean(cv_community, na.rm = TRUE), 
                     cv_comm_post_sd = sd(cv_community, na.rm = TRUE))
  
  
  # Export results in a .cvs
  
  write.csv(cv_community_summ, file=paste0(dir, 'Data/L2/cv_community_', trait_name, run_date,'.csv'))
  
  
  # Compile all results in a single data frame:
  
  cv_data <- cv_within %>% 
    dplyr::left_join(cv_across_pops, by = "sp_id") %>%
    dplyr::left_join(cv_among_sp_summ, by = "station_id") %>%
    dplyr::mutate(across(everything(), ~ .)) %>%
    cbind(cv_community_summ)

  
  return(cv_data)
  
}


# Run function for mass --------

cv_data_mass <- cv_fun(fitMass, massData, "mass")

saveRDS(cv_data_mass, file=paste0(dir, 'Data/L2/cv_data_mass-', run_date,'.rds'))



# Run function for wing --------

cv_data_wing <- cv_fun(fitWing, wingData, "wing")

saveRDS(cv_data_wing, file=paste0(dir, 'Data/L2/cv_data_wing-', run_date,'.rds'))



  
  
