#######################
# Read data files and plot components of Fig. 1

#Ice sheets shapefiles downloaded from:
#Dalton, A.S., Dulfer, H.E., Margold, M., Heyman, J., Clague, J.J., Froese, D.G., et al. (2023). 
#Deglaciation of the north American ice sheet complex in calendar years based on a 
#comprehensive database of chronological data: NADI-1. 
#Quaternary Science Reviews, 321, 108345.
#######################


# load packages -----------------------------------------------------------

library(sf)
library(rnaturalearth)
library(tmap)
library(smoothr)
library(terra)
library(ggplot2)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
cv_data_date <- '2025-01-08'
bl_run_date <- '2024-12-13'  #species ranges
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------

#cv data
cv_data_covs <- readRDS(paste0(dir, 'Data/L2/cv_data_covs-', cv_data_date, '.rds'))

#range of species
BL_data_filt <- readRDS(paste0(dir, 'Data/L1/BOTW_ranges-', bl_run_date, '.rds'))

#ice sheet shapefiles
maxGla22ka <- sf::read_sf(paste0(dir, 'Data/L0/IceSheets_Dalton_etal_2023/22ka_cal_OPTIMAL_NADI-1_Dalton_etal_QSR.shp'))
maxGla14ka <- sf::read_sf(paste0(dir, 'Data/L0/IceSheets_Dalton_etal_2023/14ka_cal_OPTIMAL_NADI-1_Dalton_etal_QSR.shp'))
maxGla10ka <- sf::read_sf(paste0(dir, 'Data/L0/IceSheets_Dalton_etal_2023/10ka_cal_OPTIMAL_NADI-1_Dalton_etal_QSR.shp'))
maxGla6ka <- sf::read_sf(paste0(dir, 'Data/L0/IceSheets_Dalton_etal_2023/6ka_cal_OPTIMAL_NADI-1_Dalton_etal_QSR.shp'))


#creating objects for plot with tmap
world_data <- rnaturalearth::ne_countries(scale = 'medium', returnclass = 'sf')
NorthAme <- world_data %>%
  dplyr::filter(continent == 'North America')

# Lambert Conformal Conic projection
lcc <- "+proj=lcc +lon_0=-90 +lat_0=60 +lat_1=33 +lat_2=45"

north_america_lcc <- NorthAme %>%
  sf::st_as_sf() %>%
  sf::st_transform(crs = lcc)

# Define bounding box
bb_lcc <- sf::st_bbox(c(xmin = -137, xmax = -65, ymin = 20, ymax = 75), crs = st_crs(4326)) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(crs = sf::st_crs(lcc)) %>%
  sf::st_bbox()


# Panel A) Age of populations (latitude)

# Smooth shapefiles
maxGla22ka_smooth <- smoothr::smooth(maxGla22ka, method = "ksmooth", smoothness = 200)
maxGla14ka_smooth <- smoothr::smooth(maxGla14ka, method = "ksmooth", smoothness = 200)
maxGla10ka_smooth <- smoothr::smooth(maxGla10ka, method = "ksmooth", smoothness = 200)
maxGla6ka_smooth <- smoothr::smooth(maxGla6ka, method = "ksmooth", smoothness = 400)

# Transform spatial data
maxGla22ka_lcc <- sf::st_transform(maxGla22ka_smooth, crs = lcc)
maxGla14ka_lcc <- sf::st_transform(maxGla14ka_smooth, crs = lcc)
maxGla10ka_lcc <- sf::st_transform(maxGla10ka_smooth, crs = lcc)
maxGla6ka_lcc <- sf::st_transform(maxGla6ka_smooth, crs = lcc)

# check that everything is valid
maxGla22ka_lcc <- sf::st_make_valid(maxGla22ka_lcc)
maxGla14ka_lcc <- sf::st_make_valid(maxGla14ka_lcc)
maxGla10ka_lcc <- sf::st_make_valid(maxGla10ka_lcc)
maxGla6ka_lcc <- sf::st_make_valid(maxGla6ka_lcc)

# Plot with tmap
ice_sheet_map <- tmap::tm_shape(north_america_lcc, bbox = bb_lcc) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(maxGla22ka_lcc) +
  tmap::tm_fill(col = "white", alpha = 0.5) +  
  tmap::tm_borders(col = "black", lwd = 0.8) +
  tmap::tm_shape(maxGla14ka_lcc) +
  tmap::tm_fill(col = "lightblue", alpha = 0.5) +
  tmap::tm_borders(col = "black", lwd = 0.8) +
  tmap::tm_shape(maxGla10ka_lcc) +
  tmap::tm_fill(col = "lightpink", alpha = 0.5) +
  tmap::tm_borders(col = "black", lwd = 0.8) +
  tmap::tm_shape(maxGla6ka_lcc) +
  tmap::tm_fill(col = "lightgreen", alpha = 0.5) +
  tmap::tm_borders(col = "black", lwd = 0.8) +
  tmap::tm_graticules(alpha = 0.6,
                      x = seq(-180, -50, by = 20),  
                      y = seq(10, 80, by = 10),  
                      lwd = 0.4)


# Save as PDF for Illustrator
tmap::tmap_save(ice_sheet_map, filename = paste0(dir, 'Results/ice_sheet_map-', run_date, '.pdf'), 
                width = 8, height = 7)



# Panel B) Distance to range edge

#choose one species to represent in the map
sp_range <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Thryothorus ludovicianus')

#get the coordinates of the species
stations_data <- cv_data_covs %>%
  dplyr::filter(sci_name == 'Thryothorus ludovicianus') %>%
  dplyr::distinct(station_id, .keep_all = TRUE) %>%
  dplyr::filter(distance_edge_km_buffer0 > 20 & distance_edge_km_buffer0 < 25 | 
                  distance_edge_km_buffer0 == max(distance_edge_km_buffer0)) %>%
  sf::st_as_sf(coords = c('lng', 'lat'), crs = 4326) %>%
  sf::st_transform(crs = lcc)
  
# Transform spatial data
sp_range_lcc <- sf::st_transform(sp_range, crs = lcc)

# Plot with tmap
dist_range_map <- tmap::tm_shape(north_america_lcc, bbox = bb_lcc) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(sp_range_lcc) +
  tmap::tm_fill(col = 'black', alpha = 0.2) +
  tmap::tm_borders(col = "black", lwd = 0.3) +
  tmap::tm_shape(stations_data) +
  tmap::tm_dots(size = 0.6, 
                col = 'black') +  # stations
  tmap::tm_graticules(alpha = 0.6,
                      x = seq(-180, -50, by = 20),  
                      y = seq(10, 80, by = 10),  
                      lwd = 0.4)


# Save as PDF for Illustrator
tmap::tmap_save(dist_range_map, filename = paste0(dir, 'Results/dist_range_map-', run_date, '.pdf'), 
                width = 8, height = 7)



# Panel C) Spatial variation

#high variation
highVar <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
# Assign random values
values(highVar) <- runif(terra::ncell(highVar), min = 0, max = 100)

#low variation
lowVar <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
# Assign random values
values(lowVar) <- runif(terra::ncell(lowVar), min = 0, max = 40)

# Convert rasters to data frames for ggplot2
highVar_df <- as.data.frame(highVar, xy = TRUE, na.rm = TRUE) %>%
  dplyr::mutate(raster_type = "High Variation")

lowVar_df <- as.data.frame(lowVar, xy = TRUE, na.rm = TRUE) %>%
  dplyr::mutate(raster_type = "Low Variation")

# Combine the data frames into one
raster_df <- dplyr::bind_rows(highVar_df, lowVar_df)

# Define the global min and max values for the color scale
min_val <- min(min(values(highVar), na.rm = TRUE), min(values(lowVar), na.rm = TRUE))
max_val <- max(max(values(highVar), na.rm = TRUE), max(values(lowVar), na.rm = TRUE))

# Create the plot
spatVar <- ggplot(raster_df) +
  geom_tile(aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_gradient(low = "darksalmon", high = "brown4", limits = c(min_val, max_val)) +  # Apply the same color scale for both
  facet_wrap(~raster_type, ncol = 2,
             labeller = labeller(raster_type = c("High Variation" = "High spatial variation",
                                                 "Low Variation" = "Low spatial variation"))) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'none',
        strip.text = element_text(size = 16, margin = margin(b = 2)),
        axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        panel.grid = element_blank())


# Save as PDF for Illustrator
ggplot2::ggsave(paste0(dir, 'Results/spat_var-', run_date, '.pdf'), spatVar,
                device = 'pdf', width = 8, height = 5)



# Panel D) Temporal variation

# Generate time series data
time_points <- seq(1, 100, by = 1)

# Generate a random time series with temporal variation (e.g., trend + seasonality + noise)
trend <- 0.001 * (as.numeric(time_points))
HighNoise <- rnorm(length(time_points), mean = 0, sd = 5)  
LowNoise <- rnorm(length(time_points), mean = 0, sd = 1)
  
# Combine the components to create the final time series
valueHigh <- trend + HighNoise
valueLow <- trend + LowNoise

# Create data frames and binding them for ggplot
High <- data.frame(ID = seq(1:length(time_points)),
                   tempvar = rep('highvar', length(time_points)),
                   value = valueHigh)

Low <- data.frame(ID = seq(1:length(time_points)), 
                  tempvar = rep('lowvar', length(time_points)),
                  value = valueLow)

time_series_data <- dplyr::bind_rows(High, Low)

# Plot the time series
tempVar <- ggplot(time_series_data, aes(x = ID, y = value)) +
  geom_line() +  
  facet_wrap(~tempvar, ncol = 2,
             labeller = labeller(tempvar = c('highvar' = "High temporal variation",
                                             'lowvar' = "Low temporal variation"))) +
  labs(x = NULL, 
       y = NULL) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(size = 16, margin = margin(b = 5)),
        axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


# Save as PDF for Illustrator
ggplot2::ggsave(paste0(dir, 'Results/temp_var-', run_date, '.pdf'), tempVar,
                device = 'pdf', width = 8, height = 5)



# Panel G) Range size

# Define bounding box
bb_lcc <- sf::st_bbox(c(xmin = -137, xmax = -60, ymin = 35, ymax = 55), crs = st_crs(4326)) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(crs = sf::st_crs(lcc)) %>%
  sf::st_bbox()


#choose two species to represent in the map for range size
big_range <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Cardinalis cardinalis') %>%
  sf::st_transform(crs = sf::st_crs(lcc))

small_range <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Poecile rufescens') %>%
  sf::st_transform(crs = sf::st_crs(lcc))

# Plot large range sp with tmap
big_map <- tmap::tm_shape(north_america_lcc, bbox = bb_lcc) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(big_range) +
  tmap::tm_fill(col = 'black', alpha = 0.2) +
  tmap::tm_borders(col = "black", lwd = 0) 

# Plot small range sp with tmap
small_map <- tmap::tm_shape(north_america_lcc, bbox = bb_lcc) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(small_range) +
  tmap::tm_fill(col = 'black', alpha = 0.2) +
  tmap::tm_borders(col = 'black', lwd = 0)
  

#Arrange both maps side by side
range_size <- tmap::tmap_arrange(big_map, small_map, ncol = 2)

# Save as PDF for Illustrator
tmap::tmap_save(range_size, filename = paste0(dir, 'Results/range_size_map-', run_date, '.pdf'), 
                width = 8, height = 5)



# Panel H) Migratory status

# Lambert Conformal Conic projection
lcc <- "+proj=lcc +lon_0=-90 +lat_0=60 +lat_1=33 +lat_2=45"

#choose one species to represent in the map for migratory
# 2	Breeding Season
# 3	Non-breeding Season
sp_range_breed <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Icterus galbula') %>%
  dplyr::filter(seasonal == 2) %>%
  sf::st_transform(crs = sf::st_crs(lcc))

sp_range_nonbreed <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Icterus galbula') %>%
  dplyr::filter(seasonal == 3) %>%
  sf::st_transform(crs = sf::st_crs(lcc))

#creating objects for plot with tmap
world_data <- rnaturalearth::ne_countries(scale = 'medium', returnclass = 'sf')
americas <- world_data %>%
  dplyr::filter(continent == 'North America' | continent == 'South America') %>%
  sf::st_as_sf() %>%
  sf::st_transform(crs = lcc)

# Define bounding box
bb_lcc_mig <- sf::st_bbox(c(xmin = -130, xmax = -65, ymin = -5, ymax = 75), crs = st_crs(4326)) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(crs = sf::st_crs(lcc)) %>%
  sf::st_bbox()

#choose one species to represent in the map for non-migratory
sp_range_nonmigr <- BL_data_filt %>%
  dplyr::filter(sci_name == 'Poecile gambeli') %>%
  sf::st_transform(crs = sf::st_crs(lcc))


# Plot migratory sp with tmap
mig_map <- tmap::tm_shape(americas, bbox = bb_lcc_mig) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(sp_range_breed) +
  tmap::tm_fill(col = 'red', alpha = 0.2) +
  tmap::tm_borders(col = "black", lwd = 0) +
  tmap::tm_shape(sp_range_nonbreed) +
  tmap::tm_fill(col = 'blue', alpha = 0.2) +
  tmap::tm_borders(col = "black", lwd = 0) +
  tmap::tm_graticules(alpha = 0.6,
                      x = seq(-180, -50, by = 20),  
                      y = seq(-5, 80, by = 15),  
                      lwd = 0.4)

# Plot non-migratory sp with tmap
nonmig_map <- tmap::tm_shape(americas, bbox = bb_lcc_mig) +
  tmap::tm_fill(col = "grey85") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(sp_range_nonmigr) +
  tmap::tm_fill(col = 'blue', alpha = 0.2) +
  tmap::tm_borders(col = 'black', lwd = 0.1) +
  tmap::tm_graticules(alpha = 0.6,
                      x = seq(-180, -50, by = 20),  
                      y = seq(-5, 80, by = 15),  
                      lwd = 0.4)


#Arrange both maps side by side
migr_status <- tmap::tmap_arrange(mig_map, nonmig_map, ncol = 2)

# Save as PDF for Illustrator
tmap::tmap_save(migr_status, filename = paste0(dir, 'Results/migr_status_map-', run_date, '.pdf'), 
                width = 8, height = 5)



