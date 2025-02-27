#######################
# Joining covariates data with CV data set

#for within-population across space for wing and mass
#for within-population among species for wing and mass

#######################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
cv_data_date <- '2025-01-08'
trait_date <- '2024-12-13'
dhi_date <- '2024-12-13'
rang_date <- '2024-12-15'
dir <- '~/Documents/MorphCVLatitude/' #change as needed


# load packages -----------------------------------------------------------

library(dplyr)


# read in data -----------------------------------------------------------

#cv data
cv_data_mass <- readRDS(paste0(dir, 'Data/L2/cv_data_mass-', cv_data_date,'.rds'))
cv_data_wing <- readRDS(paste0(dir, 'Data/L2/cv_data_wing-', cv_data_date,'.rds'))

#dhi variation
dhi_variation <- readRDS(paste0(dir, 'Data/L2/dhi_variation-', dhi_date,'.rds'))

#distance to range edge, range size, and migration distance
range_metrics <- readRDS(paste0(dir, 'Data/L2/dist_edge_size_migdist-', rang_date,'.rds'))

#trait data
traits_data <- readRDS(paste0(dir, 'Data/L2/trait_data-', trait_date,'.rds'))



# join wing data frames -----------------------------------------------------------

cv_data_wing <- cv_data_wing %>%
  dplyr::left_join(dhi_variation, by = 'cn_id') %>%
  dplyr::left_join(range_metrics, by = 'cn_id') %>%
  dplyr::left_join(traits_data %>% dplyr::select(-sci_name), by = 'sp_id') %>%
  dplyr::relocate(c(sci_name, lng, lat), .after = station_id)




# join mass data frames -----------------------------------------------------------

cv_data_mass <- cv_data_mass %>%
  dplyr::left_join(dhi_variation, by = 'cn_id') %>%
  dplyr::left_join(range_metrics, by = 'cn_id') %>%
  dplyr::left_join(traits_data %>% dplyr::select(-sci_name), by = 'sp_id') %>%
  dplyr::relocate(c(sci_name, lng, lat), .after = station_id)



# write out mass file
write.csv(cv_data_mass, file = paste0(dir, 'Data/L2/cv_data_covs_mass-', run_date,'.csv'),
          row.names = FALSE)

saveRDS(cv_data_mass, file = paste0(dir, 'Data/L2/cv_data_covs_mass-', run_date,'.rds'))


# write out wing file
write.csv(cv_data_wing, file = paste0(dir, 'Data/L2/cv_data_covs_wing-', run_date,'.csv'),
          row.names = FALSE)

saveRDS(cv_data_wing, file = paste0(dir, 'Data/L2/cv_data_covs_wing-', run_date,'.rds'))



# join both mass and wing data frames
cv_data_covs <- rbind(cv_data_wing, cv_data_mass)

# write out file

saveRDS(cv_data_covs, file = paste0(dir, 'Data/L2/cv_data_covs-', run_date,'.rds'))


