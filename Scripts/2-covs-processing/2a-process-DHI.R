################
## DHI

# Source: https://silvis.forest.wisc.edu/data/dhis/
# Dowload of NDVI16 (Normalized Difference Vegetation Index) for the period 2003-2015
# Resolution is 1 x 1km. 

# The DHIs are designed for biodiversity assessments and to describe habitats of different species. 
# Three individual indices comprise the DHIs:
# DHI cum – cumulative DHI, i.e., the area under the phenological curve of a year
# DHI min – minimum DHI, i.e., the minimum value of the phenological curve of a year
# DHI var – seasonality DHI, i.e., the coefficient of variation of the phenological curve of a year

# After downloading all the NDVI .tif files from https://silvis.forest.wisc.edu/data/dhis/:
################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
maps_data_date <- "2024-12-12" 
dhi_date <- "2024-12-13"       #date of processing DHI files
dir <- 'XXXX' #location of working diretory - change as needed
data_dir <- "XXXX" #location of data - change as needed


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)


# read in data -------------------------------------------------

#list files
#source: https://silvis.forest.wisc.edu/data/dhis/
#Radeloff et al. 2019 Remote Sensing of Environment
#2003-2015
#1-km resolution
tif_files <- list.files(data_dir, pattern = "\\.tif$", full.names = TRUE)

#read in all tifs
main_rast <- terra::rast(tif_files)

#read in morph data
st_data <- readRDS(paste0(dir, 'Data/L1/MAPS-master-', maps_data_date, '.rds'))


# process DHI rasters ------------------------------------------------------------
# need to run this only once

#indices for each metric
idx <- 1:(length(tif_files) * 3)
# DHI cum – cumulative DHI, i.e., the area under the phenological curve of a year
cum_idx <- idx[seq(1, max(idx), 3)]

#separate rasters
cum_rast <- main_rast[[cum_idx]]

#CV over time of cumulative NDVI
sd_cum <- terra::app(cum_rast, fun = sd)
mn_cum <- terra::app(cum_rast, fun = mean)
cv_cum <- sd_cum / mn_cum

#combine
dhi_met <- c(mn_cum, cv_cum)

#change names (dhi_cum_mean, dhi_cv_year, dhi_cv_season)
names(dhi_met) <- c('dhi_cum_mean', 'dhi_cv_year')

#plot
#plot(dhi_met) #very small cv_cum

#write to file
#dhi_cum_mean = mean DHI for each pixel
#dhi_cv_year = temporal CV for DHI for each pixel
terra::writeRaster(dhi_met, paste0(dir, 'Data/L1/DHI_rast-', run_date,'.tif'), overwrite = TRUE)


# reproject raster and stations data -----------------------------------------------------

# mean DHI for each station/species buffer (mean_dhi_cum_mean)
# sd DHI for each station/species (sd_dhi_cum_mean)
# cv_dhi_cum_mean = sd_dhi_cum_mean / mean_dhi_cum_mean

dhi_met <- terra::rast(paste0(dir, 'Data/L1/DHI_rast-', dhi_date,'.tif'))

#albers_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
eqArea <- "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"  #EPGS code for Equal Area projection

# Reproject DHI raster - Takes a long time on laptop
dhi_met_ea <- terra::project(dhi_met, eqArea)

# Save reprojected raster
terra::writeRaster(dhi_met_ea, paste0(dir, 'Data/L1/DHI_rast_ea-', run_date, '.tif'), 
                   overwrite = TRUE)

# Reproject stations 
stations_data <- dplyr::select(st_data, lng, lat, cn_id) %>%
  dplyr::distinct() %>%
  # Convert the station locations to an sf object with the initial CRS
  sf::st_as_sf(coords = c("lng", "lat"), crs = 4326) %>%
  #reproject station lat/lon
  sf::st_transform(eqArea)


# spatial variability -----------------------------------------------------

# load reprojected DHI data
dhi_met_ea <- terra::rast(paste0(dir, 'Data/L1/DHI_rast_ea-', dhi_date, '.tif'))


# Function to extract DHI mean, sd, and calculate cv

dhi_fun <- function(raster_object, stations_object, buffer_size, names_cols){
  
  # Create buffers around each station
  buffers <- sf::st_buffer(stations_object, dist = buffer_size)
  
  #mean dhi vals for each buffer
  mean_values <- terra::extract(raster_object, buffers, mean, na.rm=TRUE)
  colnames(mean_values) <- c("cn_id", "mean_dhi_cum_mean", "mean_dhi_cv_year")
  
  #sd dhi vals for each buffer
  sd_values <- terra::extract(raster_object, buffers, sd, na.rm=TRUE)
  colnames(sd_values) <- c("cn_id", "sd_dhi_cum_mean", "sd_dhi_cv_year")
  
  #calc cv
  dhi_spat_cv <- dplyr::left_join(mean_values, sd_values, by = "cn_id") %>%
    dplyr::mutate(cv_dhi_cum_mean = sd_dhi_cum_mean / mean_dhi_cum_mean)
  
  colnames(dhi_spat_cv) <- names_cols
  
  return(dhi_spat_cv)
  
}


#10km buffer
dhi_spat_cv_10km <- dhi_fun(raster_object = dhi_met_ea, 
                            stations_object = stations_data, 
                            buffer_size = 10000,
                            names_cols = c('cn_id', 
                                           'mean_dhi_cum_mean_10km', 
                                           'mean_dhi_cv_year_10km',
                                           'sd_dhi_cum_mean_10km', 
                                           'sd_dhi_cv_year_10km',
                                           'cv_dhi_cum_mean_10km'))
  


#merge spatial and temporal data
#cv_dhi_cum_mean_XXXX = spatial
#mean_dhi_cv_year_XXXX = temporal
dhi_variation <- dhi_spat_cv_10km %>%
  dplyr::select(cn_id,
                dplyr::starts_with('cv_dhi_cum_mean'),
                dplyr::starts_with('mean_dhi_cv_year')) %>%
  dplyr::rename(DHI_spatial_var_10km = cv_dhi_cum_mean_10km,
                DHI_temporal_var_10km = mean_dhi_cv_year_10km)


# write out -----------------------------------

write.csv(dhi_variation, file = paste0(dir, 'Data/L2/dhi_variation-', run_date,'.csv'),
          row.names = FALSE)

saveRDS(dhi_variation, file = paste0(dir, 'Data/L2/dhi_variation-', run_date,'.rds'))


