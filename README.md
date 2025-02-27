# WithinPopsPhenoVar

Repository containing the analysis to investigate how different mechanisms shape the within-population intraspecific morphological variation in North American birds.

## Associated publications:

Zulian, V and C Youngflesh. Drivers of phenotypic variability: How multiple mechanisms shape variation within populations of North American birds. In Review

## Repo structure:

* `Scripts/`
  * `1-CV-calculation/`
    * `1a-cv-modeling-NIMBLE.R` - read in `MAPS-master.rds` file and estimate hierarchically the CV for each species/station combinations using a NIMBLE model.
    * `1b-cv-processing.R` - process the results from `1a-cv-modeling-NIMBLE.R`, calculate the CV for multiple levels of biological organization (within sp, across sp, across stations, etc).
  * `2-covs-processing/`
    * `2a-process-DHI.R` - read in DHI files from Radeloff et al. 2019 *Remote Sensing of Environment*. Calculate spatial and temporal variation as `dhi-YYYY-MM-DD-.csv`
    * `2b-process-range-size-distances.R` - read in BOTW species range, update names and calculate range size, migration distance, and distance to range edge (using different buffer sizes to remove the lakes and other holes in the distribution. Results saved as `dist_edge_size_migdist-YYYY-MM-DD.rds`.
    * `2c-process-species-traits.R` - read in AVONET data (Tobias et al. 2022 *Ecology Letters*) and Bird et al. 2020 *Conservation Biology* and compile Hand Wing Index (HWI) and generation length of each species as `trait_data-YYYY-MM-DD.csv`.
    * `2d-process-phylo.R` - read in consensus trees from BirdTree (https://data.vertlife.org/birdtree/Stage2/EricsonStage2_0001_1000.zip), check names, prune trees and save a consensus tree for the species of interest as `L1/consensus_tree-YYYY-MM-DD.rds`.
    * `2e-join-cov-cv.R` - read in cv data and covs data, run PCA for spatial and temporal variation, join all in a single file as `cv_data_covs-YYYY-MM-DD.rds`.
  * `3-within-pops-var-cv.R` - run the within species variation on CV for mass and wing.
  * `4-within-pops-among-sp-var-cv.R` - run the within species among species variation on CV for mass and wing.
  * `5-results-summary-plots/` - scripts to summarize results and plot
    * `5a-study-area-stations.R`
    * `5b-results-within-pops.R`
    * `5c-results_within-pops-among-sp.R`
    * `5d-components-fig1.R`
  * `model_files/`
    * `3-within-pops-var-cv.stan` - Stan model file `3-within-pops-var-cv.R`
    * `4-within-pops-among-sp.stan` - Stan model file for `4-within-pops-among-sp-var-cv.R`
  * `5-results-processing-figs.R` - plot results, save figures
* `Data/` (ignored)
  * `L0/`
    * `BLIOCPhyloMasterTax.csv` - taxonomy info from BIRDNET. Used to process the phylogenetic trees.
    * `BOTW_2024/` - range maps for all species from http://datazone.birdlife.org/species/requestdis. Download date: 2024-11-22.
    * `mnt/` - phylogenetic trees from https://data.vertlife.org/birdtree/Stage2/EricsonStage2_0001_1000.zip. Download date: 2024-12-13.
    * `DHI/` - Dynamic Habitat Indices from Radeloff et al. 2019 *Remote Sensing of Environment* (https://doi.org/10.1016/j.rse.2018.12.009)
    * `Traits/` - different .csv files containing trait information.
     * `Bird_et_al_2020_Con_Bio.csv` - generation length data from Bird et al. 2020 *Conservation Biology* (https://doi.org/10.1111/cobi.13486)
     * `Tobias_et_al_2022_Eco_Letters.csv` - AVONET trait database from Tobias et al. 2022 *Ecology Letters* (https://doi.org/10.1111/ele.13898)
  * `L1/`
    * `MAPS-master-YYYY-MM-DD.rds` - MAPS data
    * `DHI_rast-2024-12-12.tif` - DHI raster, mean productivity and temporal CV productivity. Intermediate product `2a-process-DHI.R`.
    * `DHI_rast-ea-2024-12-12.tif` - Same as `DHI_rast-2024-12-12.tif` but equal area projection. Intermediate product `2a-process-DHI.R`.
    * `consensus_tree-YYYY-MM-DD.rds` - consensus tree for the species of interest, from `2d-process-phylo.R`.
    * `phylo_names-YYYY-MM-DD.csv` - name key for filtering species in the consensus tree, from `2d-process-phylo.R`.
  * `L2/`
    * `dhi_variation-YYYY-MM-DD-.csv` - Spatial and temporal variability DHI from `2a-process-DHI.R`.
    * `dist_edge_size_migdist-YYYY-MM-DD.csv` - range metrics from `2b-process-range-size-distances.R`.
    * `trait_data.csv` - HWI and generation length of each species from `2c-process-species-traits.R`.
    * `cv_data_covs.rds` - complete data set, with morph and covs data. Derived from `2e-join-cov-cv.R`.
* `Results/` (ignored)
