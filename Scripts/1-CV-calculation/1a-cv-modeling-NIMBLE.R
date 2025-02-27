
### Code to calculate the CV for mass and wing translated from Stan:
### Code to estimate:
### 1. the mean and standard deviation of individual measurements of wing length and mass within populations and species:
### i = individuals
### j = stations (pops)
### k = species
### measure_ijk ~ normal(mean_measurement_jk, sd_measurement_jk);

### hierarchical estimation of mean_measurement_jk and sd_measurement_jk with a species mean
### mean_measurement_jk ~ normal(mu_k, sigma_k); 
### sd_measurement_jk ~ normal(mu_sd_k, sigma_sd_k); 
  
### Vectorized format, using cn_id = Species/station id for each observation index
  
### Individual measurements are log scaled
  
### 2. morphological variation (cv) estimation as a derived quantity of sd/mean estimated above for each each species/station combinations:
### cv_jk = sd_measurement_jk / mean_measurement_jk
    
### 3. across populations cv
### cv_across_k = sigma_k / mu_k;


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
maps_date <- "2024-12-12"
dir <- '~/Documents/MorphCVLatitude/' #change as needed


# load packages -----------------------------------------------------------

library(nimble)
library(MCMCvis)


# read in data -------------------------------------------------------------

#read in data - already processed and filtered
meas_data <- readRDS(paste0(dir, 'Data/L1/MAPS-master-', maps_date, '.rds'))

# Define the model code
modelCV <- nimbleCode({
  for (i in 1:N_obs) {
    measure[i] ~ dnorm(mean_measurement[cn_id[i]], sd = sd_measurement[cn_id[i]])
  }

  for (j in 1:N_cn_id) {  #N_cn_id: Species/station id for each observation index
    mean_measurement[j] ~ dnorm(mu[sp_id[j]], sd = sigma[sp_id[j]])    #species specific mean and sigma measurement
    sd_measurement[j] ~ dnorm(mu_sd[sp_id[j]], sd = sigma_sd[sp_id[j]])
  }

  for (i in 1:N_cn_id) {
    cv[i] <- sd_measurement[i] / mean_measurement[i]
  }

  for (k in 1:N_species) {
    mu[k] ~ dnorm(mm, sd = ms)
    sigma[k] ~ dnorm(msi, sd = msis)
    mu_sd[k] ~ dnorm(msd, sd = msdd)
    sigma_sd[k] ~ dnorm(smu, sd = smd)
    cv_acrosspops[k] <- sigma[k] / mu[k]
  }

  # Priors
  mm ~ dgamma(4, 1)
  ms ~ dgamma(1, 2)
  msi ~ dgamma(1, 2)
  msis ~ dgamma(1, 5)
  msd ~ dgamma(1, 2)
  msdd ~ dgamma(1, 5)
  smu ~ dgamma(1, 2)
  smd ~ dgamma(1, 5)
})



### MASS ---------------------------------

# Define the data
dataM <- list(measure = meas_data$log_mass)

# Create sp index for each cn_id
cn_sp <- unique(meas_data[,c('sp_id', 'cn_id', 'station_id')])

# Define the constants
constantsM <- list(N_obs = nrow(meas_data),
  N_cn_id = length(unique(meas_data$cn_id)), #sp/station combinations
  cn_id = meas_data$cn_id, #sp/station index
  N_species = length(unique(meas_data$sp_id)), 
  sp_id = cn_sp$sp_id, #sp index for each cn_id
  station_id = cn_sp$station_id) #keep station_id even though not used (important for data processing later)

# Define the initial values
initsM <- list(sd_measurement = rep(0.5, length(unique(meas_data$cn_id))),
              mean_measurement = rep(3, length(unique(meas_data$cn_id))),
              mu = rep(2, length(unique(meas_data$sp_id))),
              sigma = rep(0.5, length(unique(meas_data$sp_id))), 
              mu_sd = rep(0.5, length(unique(meas_data$sp_id))),
              sigma_sd = rep(0.5, length(unique(meas_data$sp_id))), 
              mm = 4,
              ms = 1,
              msi = 0.5,
              msis = 1,
              msd = 0.5,
              msdd = 1,
              smu = 1,
              smd = 1)

# Create the NIMBLE model
CVmodM <- nimbleModel(
  code = modelCV,
  data = dataM,
  constants = constantsM,
  inits = initsM,
  name = "modelCV"
)

# Compile the NIMBLE model
CVmodcM <- compileNimble(CVmodM)

# Configure MCMC
confM <- configureMCMC(CVmodM, monitors = c("mean_measurement", "sd_measurement", "mu", "sigma", "mu_sd", 
                                          "sigma_sd", "cv",  "cv_acrosspops", "mm", "ms", "msi", "msis", 
                                          "msd", "msdd", "smu", "smd"))

# Build the MCMC
MCMC_M <- buildMCMC(confM)

# Compile the MCMC
cMCMC_M <- compileNimble(MCMC_M, project = CVmodcM)

# Run the MCMC for mass
outputMass <- runMCMC(cMCMC_M, nchains = 4, niter = 50000, nburnin = 20000, thin = 20, 
                  summary = TRUE)


# Check the output summary -------
print(outputMass$summary)
summary <- MCMCvis::MCMCsummary(outputMass$samples, params = "mean_measurement")
summarySDMu <- MCMCvis::MCMCsummary(outputMass$samples, params = "mu_sd")
summaryMu <- MCMCvis::MCMCsummary(outputMass$samples, params = "mu")
summaryCVacross <- MCMCvis::MCMCsummary(outputMass$samples, params = "cv_acrosspops")

#Good Rhat and large n.eff.


# save summary ------

#save out summary, model fit, data
MCMCvis::MCMCdiag(outputMass$samples,
                  round = 4,
                  file_name = paste0('MassCV-output-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('MassCV-output-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('MassCV-fit-', run_date),
                  add_obj = list(c(dataM, constantsM, initsM)),
                  add_obj_names = paste0('MassCV-data-', run_date),
                  cp_file = paste0(dir,'Scripts/2-CV-calculation/2a-cv-modeling-NIMBLE.R'),
                  cp_file_names = paste0('2a-cv-modeling-NIMBLE', run_date, '.R'))


# PPO ---------------------------------------------------------------------

#Simulate values for the priors.
#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mm = rgamma(10000, 4, 1),
                    ms = rgamma(10000, 1, 2),
                    msi = rgamma(10000, 1, 2),
                    msis = rgamma(10000, 1, 5),
                    msd = rgamma(10000, 1, 2),
                    msdd = rgamma(10000, 1, 5),
                    smu = rgamma(10000, 1, 2),
                    smd = rgamma(10000, 1, 5))

MCMCtrace(outputMass$sample, 
          params = c("mm", "ms", "msi", "msis", 
                     "msd", "msdd", "smu", "smd"),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE)




  
### WING -----------------------------------
# Run the MCMC for Wing

# Define data
dataW <- list(measure = meas_data$log_wing)

# Create sp index for each cn_id
cn_sp <- unique(meas_data[,c('sp_id', 'cn_id', 'station_id')])

# Define the constants
constantsW <- list(N_obs = nrow(meas_data),
                  N_cn_id = length(unique(meas_data$cn_id)), #sp/station combinations
                  cn_id = meas_data$cn_id, #sp/station index
                  N_species = length(unique(meas_data$sp_id)), 
                  sp_id = cn_sp$sp_id,  #sp index for each obs
                  station_id = cn_sp$station_id) 

# Define the initial values
initsW <- list(sd_measurement = rep(0.5, length(unique(meas_data$cn_id))),
              mean_measurement = rep(4, length(unique(meas_data$cn_id))),
              mu = rep(2, length(unique(meas_data$sp_id))),
              sigma = rep(0.5, length(unique(meas_data$sp_id))), 
              mu_sd = rep(0.5, length(unique(meas_data$sp_id))),
              sigma_sd = rep(0.5, length(unique(meas_data$sp_id))), 
              mm = 4,
              ms = 1,
              msi = 0.5,
              msis = 1,
              msd = 0.5,
              msdd = 1,
              smu = 1,
              smd = 1)

# Create the NIMBLE model
CVmodW <- nimbleModel(
  code = modelCV,
  data = dataW,
  constants = constantsW,
  inits = initsW,
  name = "modelCV"
)

# Compile the NIMBLE model
CVmodcW <- compileNimble(CVmodW)

# Configure MCMC
confW <- configureMCMC(CVmodW, monitors = c("mean_measurement", "sd_measurement", "mu", "sigma", "mu_sd", 
                                          "sigma_sd", "cv", "cv_acrosspops", "mm", "ms", "msi", "msis", 
                                          "msd", "msdd", "smu", "smd"))

# Build the MCMC
MCMC_W <- buildMCMC(confW)

# Compile the MCMC
cMCMC_W <- compileNimble(MCMC_W, project = CVmodcW)

# Run the MCMC for mass
outputWing <- runMCMC(cMCMC_W, nchains = 4, niter = 50000, nburnin = 20000, thin = 20, 
                      summary = TRUE)

# Check the output summary
print(outputWing$summary)
MCMCvis::MCMCsummary(outputWing$samples, params = "mean_measurement")
summary <- MCMCvis::MCMCsummary(outputWing$samples, params = "cv")
summaryMu <- MCMCvis::MCMCsummary(outputWing$samples, params = "mu")

#Good Rhat and large n.eff.


# save summary ------

#save out summary, model fit, data
MCMCvis::MCMCdiag(outputWing$samples,
                  round = 4,
                  file_name = paste0('WingCV-output-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('WingCV-output-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('WingCV-fit-', run_date),
                  add_obj = list(c(dataW, constantsW, initsW)),
                  add_obj_names = paste0('WingCV-data-', run_date),
                  cp_file = paste0(dir,'Scripts/2-CV-calculation/2a-cv-modeling-NIMBLE.R'),
                  cp_file_names = paste0('2a-cv-modeling-NIMBLE', run_date, '.R'))



# PPO ---------------------------------------------------------------------

#Simulate values for the priors.
#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mm = rgamma(10000, 4, 1),
                    ms = rgamma(10000, 1, 2),
                    msi = rgamma(10000, 1, 2),
                    msis = rgamma(10000, 1, 5),
                    msd = rgamma(10000, 1, 2),
                    msdd = rgamma(10000, 1, 5),
                    smu = rgamma(10000, 1, 2),
                    smd = rgamma(10000, 1, 5))

MCMCtrace(outputWing$sample, 
          params = c("mm", "ms", "msi", "msis", 
                     "msd", "msdd", "smu", "smd"),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE)

