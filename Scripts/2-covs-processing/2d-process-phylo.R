#######################
# Phylogenetic data

# Download consensus tree from BirdTree
#https://data.vertlife.org/birdtree/Stage2/EricsonStage2_0001_1000.zip
#Download date: 2024-12-13
#######################


# load packages -----------------------------------------------------------

library(ape)
library(phytools)
library(dplyr)
library(tidyr)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
#tree_date <- '2024-12-13'
maps_data_date <- '2024-12-12'
#phylo_date <- '2024-12-15'
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------------------------------------------------

#read in morph data
st_data <- readRDS(paste0(dir, 'Data/L1/MAPS-master-', maps_data_date, '.rds'))

#load trees
BT_tree <- ape::read.tree(paste0(dir, 'Data/L0/mnt/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/AllBirdsEricson1.tre'))

#load taxonomy info - FROM BIRDNET
bt_phylo <- read.csv(paste0(dir, 'Data/L0/BLIOCPhyloMasterTax.csv'))

#get unique species
MAPS_usp <- unique(gsub(' ', '_', st_data$sci_name))


# change names to match birdtree.org --------------------------------------

#where matches are
nin <- which(MAPS_usp %in% bt_phylo$TipLabel)

#names df
names_df <- data.frame(MAPS_sci_name = MAPS_usp, BT_sci_name = NA)
names_df$BT_sci_name[nin] <- MAPS_usp[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(BT_sci_name))

#data.frame of matches
match <- data.frame(MAPS_sci_name = nm$MAPS_sci_name,
                    BT_sci_name = 
                      c('Carduelis_flammea',
                        'Wilsonia_canadensis',
                        'Wilsonia_pusilla',
                        'Regulus_calendula',
                        'Picoides_nuttallii',
                        'Picoides_pubescens',
                        'Oporornis_formosus',
                        'Oporornis_philadelphia',
                        'Oporornis_tolmiei',
                        'Carpodacus_cassinii',
                        'Carpodacus_mexicanus',
                        'Carpodacus_purpureus',
                        'Vermivora_celata',
                        'Vermivora_luciae',
                        'Vermivora_peregrina',
                        'Vermivora_ruficapilla',
                        'Vermivora_virginiae',
                        'Pipilo_crissalis',
                        'Seiurus_motacilla',
                        'Seiurus_noveboracensis',
                        'Parus_atricapillus',
                        'Parus_gambeli',
                        'Parus_rufescens',
                        'Parula_americana',
                        'Wilsonia_citrina',
                        'Dendroica_coronata',
                        'Dendroica_discolor',
                        'Dendroica_magnolia',
                        'Dendroica_occidentalis',
                        'Dendroica_pensylvanica',
                        'Dendroica_petechia',
                        'Dendroica_townsendi',
                        'Carduelis_pinus',
                        'Carduelis_psaltria',
                        'Troglodytes_troglodytes',
                        'Vermivora_pinus'))

#fill in gaps
for (i in 1:NROW(match)){
  #i <- 1
  names_df$BT_sci_name[which(names_df$MAPS_sci_name == match$MAPS_sci_name[i])] <- match$BT_sci_name[i]
}


sci_names_phylo <- dplyr::left_join(names_df, 
                                    bt_phylo, 
                                    by = c('BT_sci_name' = 'TipLabel')) %>%
  dplyr::select(MAPS_sci_name, BT_sci_name, 
                common_name = English, 
                family_latin = BLFamilyLatin, 
                family_english = BLFamilyEnglish, 
                order = IOCOrder)




# prune trees ----------------------------------------------------------------

prune_trees <- list()

for(i in 1:length(BT_tree)){
  
  #species to species to drop from tree 
  nm <- setdiff(BT_tree[[i]]$tip.label, sci_names_phylo$BT_sci_name)
  
  #prune specified tips from tree
  prune_trees[[i]] <- ape::drop.tip(BT_tree[[i]], nm)
  
  #print(length(prune_trees[[i]]$tip.label)) #make sure all trees have the number of species (99)
  
}


#get consensus tree
consen_tree <- phytools::consensus.edges(prune_trees)

#check that names match
table(sci_names_phylo$BT_sci_name %in% consen_tree$tip.label)


#change names of tree to match the maps data

# Match names to keep the same order
match_names <- match(consen_tree$tip.label, sci_names_phylo$BT_sci_name)

#double check that the order is correct before replacing names
tiptest <- sci_names_phylo$MAPS_sci_name[match_names]
cbind(consen_tree$tip.label, tiptest) 

# Replace tree tip label to maps name 
consen_tree$tip.label <- sci_names_phylo$MAPS_sci_name[match_names]
cbind(consen_tree$tip.label, tiptest) #another check for updated names


# write out ---------------------------------------------------------------

write.csv(sci_names_phylo, file = paste0(dir, 'Data/L1/phylo_names-', 
                                         run_date, '.csv'), row.names = FALSE)

saveRDS(consen_tree, file = paste0(dir, 'Data/L1/consensus_tree-', run_date, '.rds'))



