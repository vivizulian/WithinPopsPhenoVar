#######################
## Code to compile species traits ##
# Data sets downloaded from:
# 1) AVONET: Tobias et al, 2022, Eco Letters (https://doi.org/10.1111/ele.13898)
# 2) Bird et al, 2020, Cons Bio (https://doi.org/10.1111/cobi.13486)
#######################


# Load packages ------------------------------------------

library(dplyr)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
maps_run_date <- "2024-12-12" 
dir <- '~/Documents/MorphCVLatitude/'


# read in data -----------------------------------------------------------

# species data
maps_data <- readRDS(paste0(dir, 'Data/L1/MAPS-master-', maps_run_date, '.rds'))
sci_name <- unique(maps_data$sci_name)


# Data set 1: AVONET - Tobias et al 2022 (HWI) ----------------

tobias_eco_letters <- read.csv(paste0(dir, "Data/L0/Traits/Tobias_et_al_2022_Eco_Letters.csv"))

tobias_sci_name <- unique(tobias_eco_letters$Species1) #check how names are written

# Capitalize the first letter of the first word of species name:
tobias_eco_letters$Species1 <- gsub("^(\\w)(.*)", "\\U\\1\\E\\2", tobias_eco_letters$Species1, perl = TRUE)
tobias_sci_name <- unique(tobias_eco_letters$Species1)

#where matches are
nin <- which(sci_name %in% tobias_sci_name)

#names df
names_df <- data.frame(MAPS_sci_name = sci_name, 
                       TOBIAS_sci_name = NA)

names_df$TOBIAS_sci_name[nin] <- sci_name[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(TOBIAS_sci_name))

# nm[1] = Regulus calendula -> Corthylio calendula
tobias_sci_name[grep('Regulus calendula', tobias_sci_name)]
idx2 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[1])
names_df$TOBIAS_sci_name[idx2] <- 'Corthylio calendula'

tobias_eco_letters[tobias_eco_letters$Species1=='Regulus calendula', 'Species1'] <- 'Corthylio calendula'


# nm[1] = Icterus bullockiorum -> Icterus bullockii
tobias_sci_name[grep('Icterus bullockiorum', tobias_sci_name)]
idx2 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[2])
names_df$TOBIAS_sci_name[idx2] <- 'Icterus bullockii'

tobias_eco_letters[tobias_eco_letters$Species1=='Icterus bullockiorum', 'Species1'] <- 'Icterus bullockii'

# Trait data for subset of species
trait_data_tobias <- dplyr::filter(tobias_eco_letters, Species1 %in% names_df$TOBIAS_sci_name)


# Extract Hand.Wing.Index - rename to HWI
# there are other traits that could be extracted if needed

trait_data_tobias <- trait_data_tobias %>% dplyr::select(Species1,
                                                         Hand.Wing.Index)
colnames(trait_data_tobias) <- c('sci_name', "HWI")




# Data set 2: Bird et al. 2020 (generation length) ----------------

bird_con_bio <- read.csv(paste0(dir, "Data/L0/Traits/Bird_et_al_2020_Con_Bio.csv"))

bird_sci_name <- unique(bird_con_bio$Sci_name) #check how names are written

#where matches are
nin <- which(sci_name %in% bird_sci_name)

#names df
names_df <- data.frame(MAPS_sci_name = sci_name, BIRD_sci_name = NA)
names_df$BIRD_sci_name[nin] <- sci_name[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(BIRD_sci_name))

# Regulus calendula -> Corthylio calendula
bird_sci_name[grep('Regulus calendula', bird_sci_name)]
idx2 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[1])
names_df$BIRD_sci_name[idx2] <- 'Corthylio calendula'

bird_con_bio[bird_con_bio$Sci_name=='Regulus calendula', 'Sci_name'] <- 'Corthylio calendula'


# Icterus bullockiorum -> Icterus bullockii
bird_sci_name[grep('Icterus bullockiorum', bird_sci_name)]
idx2 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[2])
names_df$BIRD_sci_name[idx2] <- 'Icterus bullockii'

bird_con_bio[bird_con_bio$Sci_name=='Icterus bullockiorum', 'Sci_name'] <- 'Icterus bullockii'


# Trait data for subset of species
trait_data_bird <- dplyr::filter(bird_con_bio, Sci_name %in% names_df$BIRD_sci_name)


# Extract GenLength and rename columns

trait_data_bird <- trait_data_bird %>% dplyr::select(Sci_name, GenLength)

colnames(trait_data_bird) <- c('sci_name', 'GenLength')

colnames(names_df)[1] <- 'sci_name'

traitData <- names_df %>%
  dplyr::select(-BIRD_sci_name) %>%
  dplyr::left_join(trait_data_tobias, by = 'sci_name') %>%
  dplyr::left_join(trait_data_bird, by = 'sci_name')


#match sci_name from traitData with sp_id from maps_data

maps_sci_name_id <- unique(maps_data[c('sci_name', 'sp_id')])

traitDataID <- maps_sci_name_id %>%
  dplyr::left_join(traitData, by = 'sci_name') 



# write out ---------------------------------------------------------------

# Export trait data into a .cvs file
write.csv(traitDataID, file=paste0(dir, 'Data/L2/trait_data-', run_date,'.csv'))

# Export trait data into a .rds file
saveRDS(traitDataID, file=paste0(dir, 'Data/L2/trait_data-', run_date,'.rds'))



