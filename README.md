# Geography, environmental conditions, and life history shape patterns of within-population phenotypic variation in North American birds


## Overview

Repository containing the analysis to investigate how different mechanisms shape within-population intraspecific phenotypic variation in North American birds.

## Associated publications

Zulian, V. and C. Youngflesh. [Geography, environmental conditions, and life history shape patterns of within-population phenotypic variation in North American birds](https://onlinelibrary.wiley.com/doi/10.1111/ele.70244). **_Ecology Letters_** 28:e70244.

## Data

* Bird, J. P., R. Martin, H. R. Akçakaya, J. Gilroy, I. J. Burfield, S. Garnett, A. Symes, J. Taylor, Ç. H. Şekercioğlu, and S. H. M. Butchart. 2020. Generation lengths of the world’s birds and their implications for extinction risk. *Conservation Biology*. 34(5), 1252-1261. DOI: [doi.org/10.1111/cobi.13486](doi.org/10.1111/cobi.13486)

* Tobias, J. A., C. Sheard, A. L. Pigot, A. J. M. Devenish, J. Yang, F. Sayol, M. H. C. Neate‐Clegg, N. Alioravainen, T. L. Weeks, R. A. Barber, P. A. Walkden, H. E. A. MacGregor, S. E. I. Jones, C. Vincent, A. G. Phillips, N. M. Marples, F. A. Montaño‐Centellas, V. Leandro‐Silva, S. Claramunt, B. Darski, B. G. Freeman, T. P. Bregman, C. R. Cooney, E. C. Hughes, E. J. R. Capp, Z. K. Varley, N. R. Friedman, H. Korntheuer, A. Corrales‐Vargas, C. H. Trisos, B. C. Weeks, D. M. Hanz, T. Töpfer, G. A. Bravo, V. Remeš, L. Nowak, L. S. Carneiro, A. J. Moncada R., B. Matysioková, D. T. Baldassarre, A. Martínez‐Salinas, J. D. Wolfe, P. M. Chapman, B. G. Daly, M. C. Sorensen, A. Neu, M. A. Ford, R. J. Mayhew, L. Fabio Silveira, D. J. Kelly, N. N. D. Annorbah, H. S. Pollock, A. M. Grabowska‐Zhang, J. P. McEntee, J. Carlos T. Gonzalez, C. G. Meneses, M. C. Muñoz, L. L. Powell, G. A. Jamie, T. J. Matthews, O. Johnson, G. R. R. Brito, K. Zyskowski, R. Crates, M. G. Harvey, M. Jurado Zevallos, P. A. Hosner, T. Bradfer‐Lawrence, J. M. Maley, F. G. Stiles, H. S. Lima, K. L. Provost, M. Chibesa, M. Mashao, J. T. Howard, E. Mlamba, M. A. H. Chua, B. Li, M. I. Gómez, N. C. García, M. Päckert, J. Fuchs, J. R. Ali, E. P. Derryberry, M. L. Carlson, R. C. Urriza, K. E. Brzeski, D. M. Prawiradilaga, M. J. Rayner, E. T. Miller, R. C. K. Bowie, R. Lafontaine, R. P. Scofield, Y. Lou, L. Somarathna, D. Lepage, M. Illif, E. L. Neuschulz, M. Templin, D. M. Dehling, J. C. Cooper, O. S. G. Pauwels, K. Analuddin, J. Fjeldså, N. Seddon, P. R. Sweet, F. A. J. DeClerck, L. N. Naka, J. D. Brawn, A. Aleixo, K. Böhning‐Gaese, C. Rahbek, S. A. Fritz, G. H. Thomas, and M. Schleuning. 2022. AVONET: morphological, ecological and geographical data for all birds. *Ecology Letters* 25:581–597. DOI: [https://doi.org/10.1111/ele.13898](https://doi.org/10.1111/ele.13898)

* Radeloff, V. C., M. Dubinin, N. C. Coops, A. M. Allen, T. M. Brooks, M. K. Clayton, G. C. Costa, C. H. Graham, D. P. Helmers, A. R. Ives, D. Kolesov, A. M. Pidgeon, G. Rapacciuolo, E. Razenkova, N. Suttidate, B. E. Young, L. Zhu, and M. L. Hobi. 2019. The Dynamic Habitat Indices (DHIs) from MODIS and global biodiversity. *Remote Sensing of Environment* 222:204–214. DOI: [https://doi.org/10.1016/j.rse.2018.12.009]

* BirdLife International and Handbook of the Birds of the World. (2022). Bird species distribution maps of the world. Version 2022.2. Available at [http://datazone.birdlife.org/species/requestdis](http://datazone.birdlife.org/species/requestdis).
  
* Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. *Nature* 491:444–448. DOI: [https://doi.org/10.1038/nature11631](https://doi.org/10.1038/nature11631)

## Workflow

The workflow for this repository involves data read in as raw data at Level 0 (L0). Level 1 (L1) data represent cleaned L0 data. Level 2 (L2) data are data merged or otherwise derived from two or more L1 data, etc. Data citations with associated DOIs are given above. Data must be organized according to the below Repository structure (i.e., in `Data/`). Scripts are designed to be run in sequence from `1-XXXX.R` -> `2-XXXX.R` -> ... The directory where the data are located (i.e., the parent directory of `Data/`) as well as the run date for each data product need to be specified at the top of associated scripts (in the `set dirs` section) before running.


## Repository structure

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
* `Data/` (ignored)
  * `L0/`
    * `BLIOCPhyloMasterTax.csv` - taxonomy info from BirdTree.org. Used to process the phylogenetic trees.
    * `BOTW_2024/` - range maps for all species from http://datazone.birdlife.org/species/requestdis. Download date: 2024-11-22.
    * `mnt/` - phylogenetic trees from BirdTree.org. Download date: 2024-12-13.
    * `DHI/` - Dynamic Habitat Indices from Radeloff et al. 2019 *Remote Sensing of Environment* (https://doi.org/10.1016/j.rse.2018.12.009)
    * `Traits/` - different .csv files containing trait information.
      * `Bird_et_al_2020_Con_Bio.csv` - generation length data from Bird et al. 2020 *Conservation Biology* (https://doi.org/10.1111/cobi.13486)
      * `Tobias_et_al_2022_Eco_Letters.csv` - AVONET trait database from Tobias et al. 2022 *Ecology Letters* (https://doi.org/10.1111/ele.13898)
  * `L1/`
    * `MAPS-master-YYYY-MM-DD.rds` - MAPS data
    * `DHI_rast-YYYY-MM-DD.tif` - DHI raster, mean productivity and temporal CV productivity. Intermediate product `2a-process-DHI.R`.
    * `DHI_rast-ea-YYYY-MM-DD.tif` - Same as `DHI_rast-2024-12-12.tif` but equal area projection. Intermediate product `2a-process-DHI.R`.
    * `consensus_tree-YYYY-MM-DD.rds` - consensus tree for the species of interest, from `2d-process-phylo.R`.
    * `phylo_names-YYYY-MM-DD.csv` - name key for filtering species in the consensus tree, from `2d-process-phylo.R`.
  * `L2/`
    * `dhi_variation-YYYY-MM-DD.rds` - Spatial and temporal variability DHI from `2a-process-DHI.R`.
    * `dist_edge_size_migdist-YYYY-MM-DD.rds` - range metrics from `2b-process-range-size-distances.R`.
    * `trait_data-YYYY-MM-DD.rds` - HWI and generation length of each species from `2c-process-species-traits.R`.
    * `cv_data_covs-YYYY-MM-DD.rds` - complete data set, with morph and covs data. Derived from `2e-join-cov-cv.R`.
* `Results/` (ignored)

## Contact Information

Viviane Zulian - zulian.vi@gmail.com

Casey Youngflesh - cyoungf@clemson.edu
