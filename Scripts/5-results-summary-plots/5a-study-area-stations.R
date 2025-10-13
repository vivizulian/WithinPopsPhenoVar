
#######################
# Read data files and plot coordinates for each banding station

#######################


# load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(raster)
library(tmap)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
cv_data_date <- '2025-01-08'
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------

#cv data
cv_data_covs <- readRDS(paste0(dir, 'Data/L2/cv_data_covs-', cv_data_date, '.rds'))

#filter unique stations data
stations_data <- cv_data_covs %>%
  dplyr::distinct(station_id, .keep_all = TRUE)


#creating objects for plot with tmap

world_data <- rnaturalearth::ne_countries(scale = 'medium', returnclass = 'sf')
NorthAme <- world_data %>%
  dplyr::filter(continent == 'North America')

north_america_sf <- sf::st_as_sf(NorthAme)

# Convert points to sf
cll_sf <- sf::st_as_sf(stations_data, coords = c('lng', 'lat'), crs = 4326)

# Lambert Conformal Conic projection
lcc <- "+proj=lcc +lon_0=-90 +lat_0=60 +lat_1=33 +lat_2=45"

# Transform spatial data
north_america_lcc <- sf::st_transform(north_america_sf, crs = lcc)
cll_sf_lcc <- sf::st_transform(cll_sf, crs = lcc)

# Define bounding box
bb_lcc <- sf::st_bbox(c(xmin = -140, xmax = -60, ymin = 20, ymax = 75), crs = st_crs(4326)) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(crs = sf::st_crs(lcc)) %>%
  sf::st_bbox()


# Plot with tmap
stations_map <- tmap::tm_shape(north_america_lcc, bbox = bb_lcc) +
  tmap::tm_fill(col = "grey90") +
  tmap::tm_borders(col = "black", lwd = 0.3) +  # Country borders
  tmap::tm_shape(cll_sf_lcc) +
  tmap::tm_dots(size = 0.3, fill_alpha = 0.7, 
                col = 'black') +  # stations
  tmap::tm_graticules(alpha = 0.6,
                      x = seq(-180, -50, by = 20),  
                      y = seq(10, 80, by = 10),  
                      lwd = 0.4) +
  tmap::tm_layout(legend.bg.color = "white",
                  legend.position = c("left", "bottom"),
                  legend.text.size = 1,
                  legend.title.size = 1.2,
                  legend.frame = FALSE)


# Save as PDF for Illustrator
tmap::tmap_save(stations_map, filename = paste0(dir, 'Results/stations_map-', run_date, '.pdf'), 
                width = 8, height = 7)


