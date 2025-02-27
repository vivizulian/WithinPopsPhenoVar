#######################
# Code to process variables related to range of species
# Data source: BirdLife International and Handbook of the Birds of the World (2024) 
# Bird species distribution maps of the world

# Citation: BirdLife International and Handbook of the Birds of 
# the World (2024) Bird species distribution maps of the world. Version 2024.2. 
# Available at http://datazone.birdlife.org/species/requestdis.
# Download date: 2024-11-22

# Distance to range edge
# Range size
# Migration distance
#######################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
maps_data_date <- "2024-12-12"
bl_run_date <- "2024-12-13"

dir <- '~/Documents/MorphCVLatitude/' #change as needed


# load packages -----------------------------------------------------------

library(dplyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(geosphere)


# read in data -------------------------------------------------------

#maps data
st_data <- readRDS(paste0(dir, 'Data/L1/MAPS-master-', maps_data_date, '.rds'))

#filter unique species/station combinations
station_data <- st_data %>%
  dplyr::distinct(cn_id, .keep_all = TRUE)

#range data - takes some time to read in
BL_data <- sf::st_read(dsn = paste0(dir, 'Data/L0/BOTW_2024/BOTW.gdb'))


# prepare and save data -------------------------------------------------------------

#unique species on BOTW data set
bl_usp <- unique(BL_data$sci_name)

#unique species on MAPS data set
MAPS_sp <- unique(gsub('_', ' ', station_data$sci_name)) #remove the dash between the genus and species

#where matches are
nin <- which(MAPS_sp %in% bl_usp)

#names df
names_df <- data.frame(MAPS_sci_name = MAPS_sp, BL_sci_name = NA)
names_df$BL_sci_name[nin] <- MAPS_sp[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(BL_sci_name))


# Update the non matching name on the range map:

#Icterus bullockiorum -> Icterus bullockii
BL_data$sci_name <- gsub("Icterus bullockiorum", "Icterus bullockii", BL_data$sci_name)

#unique species on BOTW data set
bl_usp <- unique(BL_data$sci_name)

#where matches are
nin <- which(MAPS_sp %in% bl_usp)

names_df <- data.frame(MAPS_sci_name = MAPS_sp, BL_sci_name = NA)
names_df$BL_sci_name[nin] <- MAPS_sp[nin]

#no matches - empty after correcting names
(nm <- dplyr::filter(names_df, is.na(BL_sci_name)))

#filter ranges
BL_data_filt <- dplyr::filter(BL_data, sci_name %in% names_df$BL_sci_name)


# Resolving issues with Acanthis flammea - has a distribution that extends to Europe/Asia 
# want to crop only the North America

curr.range <- BL_data_filt %>% dplyr::filter(sci_name == "Acanthis flammea")

world_data <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
NorthAme <- world_data %>%
  filter(continent == "North America")

# Check if the geometry is valid
valid_curr.range <- sf::st_is_valid(curr.range)
valid_NorthAme <- sf::st_is_valid(NorthAme)

# If either is not valid, make it valid
if (!all(valid_curr.range)) {
  curr.range <- sf::st_make_valid(curr.range)
}

if (!all(valid_NorthAme)) {
  NorthAme <- sf::st_make_valid(NorthAme)
}

# Perform the clipping using st_intersection 
cropped_shape <- sf::st_intersection(curr.range, NorthAme)

# Plot to check that it is working
world_data <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
GlobMap <- world_data %>%
  ggplot() +
  geom_sf()  #GlobMap is a gg ggplot file

GlobMap +
  geom_sf(fill = "grey95") +
  geom_sf(data = cropped_shape, fill = "lightgreen", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add the cropped shape to the object BL_data_filt

#remove the old data corresponding to "Acanthis flammea"
BL_data_filt <- BL_data_filt[BL_data_filt$sci_name != "Acanthis flammea", ]

#add the new data to BL_data_filt
BL_data_filt <- rbind(BL_data_filt, cropped_shape[,c(1:18,187)])


# Save the cleaned files

#range maps
saveRDS(BL_data_filt, paste0(dir, 'Data/L1/BOTW_ranges-', run_date, '.rds'))

#species names key
write.csv(names_df, file = paste0(dir, 'Data/L1/range_names_key-', 
                                  run_date, '.csv'), row.names = FALSE)



# calculate metrics -------------------------------------------------------------
# range size, distance to range edge and migration distance

# Key for 'seasonal' column (from BOTW):
# 1	Resident	The species is/was known or thought very likely to be resident throughout the year.
# 2	Breeding Season	The species is/was known or thought very likely to occur regularly during the breeding season and to breed or to be capable of breeding.
# 3	Non-breeding Season 	The species is/was known or thought very likely to occur regularly during the non-breeding season. In the Eurasian and North American contexts, this encompasses ‘winter’.
# 4	Passage	The species is/was known or thought very likely to occur regularly during a relatively short period(s) of the year on migration between breeding and non-breeding ranges.
# 5	Seasonal Occurrence Uncertain	The species is/was present, but it is not known if it is present during part or all of the year.

# To calculate the range size, the range polygon should include Resident, Breeding, and Non-breeding polygons.


# Load filtered data set
BL_data_filt <- readRDS(paste0(dir, 'Data/L1/BOTW_ranges-', bl_run_date, '.rds'))


#read in modified csv (names key)
names_df <- read.csv(file = paste0(dir, 
                                       'Data/L1/range_names_key-', bl_run_date, '.csv'))

# get unique species
spID <- names_df$BL_sci_name
station_data$sci_name <- gsub('_', ' ', station_data$sci_name) #remove the dash between the genus and species


range_fun <- function(sp_ID, range_data, buffer_size_holes, station_data, name_col_buffer){
  
  #create objects to hold results
  #range size
  range_size <- data.frame(matrix(NA, nrow = length(spID), ncol = 2))
  colnames(range_size) <- c('sci_name', 'range_size_km2')
  
  #distance to range edge
  dist_range_edge <- vector("list", length(spID))
  names(dist_range_edge) <- spID 
  
  #migration distance
  migrat_dist <- data.frame(matrix(NA, nrow = length(spID), ncol = 2))
  colnames(migrat_dist) <- c('sci_name', 'migration_dist_km')
  
  # Enter loop for range size, distance to range edge, and migration distance
  counter <- 1
  for (i in 1:length(spID)){
    print(paste0("Currently on species ", i, " out of ", length(spID)))
    curr.sp <- spID[i]
    curr.range <- range_data %>% dplyr::filter(sci_name == curr.sp)
    
    eqArea <- "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"  #EPGS code for Equal Area projection
    curr.range <- sf::st_transform(curr.range, crs = eqArea)
    
    #if more than one polygon (resident + breeding ranges), merge them
    if (NROW(curr.range) > 1){
      curr.range <- dplyr::filter(curr.range, seasonal <= 3)
      
      # Check and fix invalid geometries
      if (any(!sf::st_is_valid(curr.range))) {
        curr.range <- sf::st_make_valid(curr.range)
      }
      
      # Save CRS for later
      tcrs <- sf::st_crs(curr.range)
      sf::st_crs(curr.range) <- NA  # Set CRS to missing
      
      # Attempt to union the geometries
      curr.range2 <- try(sf::st_union(curr.range), silent = TRUE)
      
      # Reassign CRS
      if (!inherits(curr.range2, "try-error")) {
        sf::st_crs(curr.range2) <- tcrs
      } else {
        # If union failed, attempt to fix invalid geometries and retry
        curr.range <- sf::st_make_valid(curr.range)
        curr.range2 <- try(sf::st_union(curr.range), silent = TRUE)
        
        if (!inherits(curr.range2, "try-error")) {
          sf::st_crs(curr.range2) <- tcrs
        } else {
          # If it still fails, disable S2 processing and try again
          sf::sf_use_s2(FALSE)
          curr.range2 <- try(sf::st_union(curr.range), silent = TRUE)
          sf::sf_use_s2(TRUE)
          
          if (!inherits(curr.range2, "try-error")) {
            sf::st_crs(curr.range2) <- tcrs
          } else {
            stop("Failed to union geometries: ", attr(curr.range2, "condition")$message)
          }
        }
      }
    } else {
      curr.range2 <- curr.range
    }

    # range size - km^2
    rsize_km2 <- as.numeric(round(sf::st_area(curr.range2) / 1000^2, 0))
    
    #fill df - current id, range size
    range_size[counter,] <- c(curr.sp, rsize_km2)
    
    
    # distance to range edge - km
    stations <- dplyr::filter(station_data, sci_name == spID[i])
    
    # Specify the initial CRS for the station locations
    initial_crs <- 4326  # EPSG code for WGS 84
    
    # Select only the necessary columns (longitude, latitude, and any identifier)
    stations <- stations[, c("sci_name","lng", "lat", "cn_id")]  
    
    # Convert the station locations to an sf object with the initial CRS
    stations_sf <- sf::st_as_sf(stations, coords = c("lng", "lat"), crs = initial_crs)
    
    # Transform the bird locations to the same CRS
    station_sf_transformed <- sf::st_transform(stations_sf, crs = eqArea)
    
    curr.range.valid <- sf::st_make_valid(curr.range2)
    
    # Remove the holes by applying a positive buffer - close gaps and holes
    curr.range2.no.holes <- sf::st_buffer(curr.range.valid, dist = buffer_size_holes)  
    
    if (!all(st_is_valid(station_sf_transformed))) {
      station_sf_transformed <- sf::st_make_valid(station_sf_transformed)
    }
    
    if (!all(st_is_valid(curr.range2.no.holes))) {
      curr.range2.no.holes <- sf::st_make_valid(curr.range2.no.holes)
    }
    
    exterior_boundary <- sf::st_boundary(curr.range2.no.holes)
    
    # Calculate the shortest distance to the boundary
    distances_to_boundary <- sf::st_distance(station_sf_transformed, exterior_boundary)
    min_distances_to_boundary <- apply(distances_to_boundary, 1, min)
    
    station_sf_transformed$min_distance_to_boundary <- min_distances_to_boundary
    
    # Check which points are inside the polygon
    inside <- sf::st_intersects(station_sf_transformed, curr.range2.no.holes, sparse = FALSE)
    inside <- apply(inside, 1, any)
    
    # subtract the buffer size and get distance in km
    stations$distance_to_range_edge_km <- (as.numeric(min_distances_to_boundary) - buffer_size_holes) / 1000
    
    # set distance to 0 for the points that are outside the range
    stations[which(!inside), "distance_to_range_edge_km"] <- 0
    dist_range_edge[[i]] <- stations
    
    station_sf_transformed$min_distance_to_boundary <- dist_range_edge[[i]]$distance_to_range_edge_km
    
    
    # migration distance - km
    
    #if more than one polygon (resident + breeding ranges), merge them
    if (NROW(curr.range) > 1){
      curr.rangeNB <- dplyr::filter(curr.range, seasonal == 1 | seasonal == 3) #non-breeding
      curr.rangeB <- dplyr::filter(curr.range, seasonal == 1 | seasonal == 2) #breeding
      
      # non-breeding range
      if (NROW(curr.rangeNB) > 1) {
        curr.rangeNB <- try(sf::st_union(sf::st_make_valid(curr.rangeNB)))
        if (!inherits(curr.rangeNB, "try-error")) {
          sf::st_crs(curr.rangeNB) <- tcrs
        }
      } 
      else {
        sf::st_crs(curr.rangeNB) <- tcrs
      }
      
      # breeding range
      if (NROW(curr.rangeB) > 1) {
        curr.rangeB <- try(sf::st_union(sf::st_make_valid(curr.rangeB)))
        if (!inherits(curr.rangeB, "try-error")) {
          sf::st_crs(curr.rangeB) <- tcrs
        }
      } 
      else {
        sf::st_crs(curr.rangeB) <- tcrs
      }
    } 
    else {
      curr.rangeNB <- curr.range
      curr.rangeB <- curr.range
      sf::st_crs(curr.range) <- tcrs
    }
    
    # get centroid - ignore warning
    cen_NB <- sf::st_make_valid(sf::st_centroid(curr.rangeNB)) %>%
      sf::st_transform(4326) %>%
      sf::st_coordinates() 
    cen_NB <- ggplot2::fortify(as.data.frame(cen_NB))
    
    # get centroid - ignore warning
    cen_B <- sf::st_make_valid(sf::st_centroid(curr.rangeB)) %>%
      sf::st_transform(4326) %>%
      sf::st_coordinates()
    cen_B <- ggplot2::fortify(as.data.frame(cen_B))
    
    # Calculate migration distance in km
    migrat_dist[counter,] <- c(curr.sp, geosphere::distm(cen_NB, cen_B) / 1000)
    
    # mapSP2 <- GlobMap +
    #   geom_sf(fill = "grey95") +
    #   geom_sf(data = curr.rangeB, fill = "lightgreen", alpha = 0.5) +
    #   geom_sf(data = curr.rangeNB, fill = "lightpink", alpha = 0.5) +
    #   geom_point(data = cen_B, aes(x = X, y = Y), col = "darkgreen") +
    #   geom_point(data = cen_NB, aes(x = X, y = Y), col = "red") +
    #   theme_bw() +
    #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    #ggsave(filename = paste0("Results/BredNonBred_range_", spID[i], "_map.png"), 
           #mapSP2, width = 10, height = 6)
    
    counter <- counter + 1
    rm(curr.range, curr.range2, curr.rangeB, curr.rangeNB)
    gc()
    
  }
  
  # compile results
  
  resultsRangeEdge <- do.call(rbind, dist_range_edge)
  
  resultsRangeEdgeSize <- resultsRangeEdge %>% 
    dplyr::left_join(range_size, by = 'sci_name') %>%
    dplyr::left_join(migrat_dist, by = 'sci_name')
  
  
  # rename column according to buffer size:
  colnames(resultsRangeEdgeSize)[5] <- name_col_buffer
  
  return(resultsRangeEdgeSize)

}


# Run function 

# buffer - size 10km
dist_edge_size_buffer10km <- range_fun(sp_ID = spID, 
                                      range_data = BL_data_filt, 
                                      buffer_size_holes = 10000, 
                                      station_data = station_data,
                                      name_col_buffer = 'distance_edge_km_buffer10')



# Create an object with values

dist_edge_size <- dist_edge_size_buffer10km %>% 
  dplyr::relocate(c(range_size_km2, migration_dist_km), .after = cn_id)



# write out files

write.csv(dist_edge_size, file = paste0(dir, 'Data/L2/dist_edge_size_migdist-', run_date, '.csv'))

saveRDS(dist_edge_size, file = paste0(dir, 'Data/L2/dist_edge_size_migdist-', run_date, '.rds'))


