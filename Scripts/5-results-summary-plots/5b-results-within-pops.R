#######################
# Read files and plot results from model #4 for mass and wing

# Model #4: Within pops variation on CV

# mu_cv = gamma[sp_id] + 
#   theta1[sp_id] * lat +           #latitude
#   theta2[sp_id] * distrange +     #distance to range edge (buffer 10km)
#   theta3[sp_id] * spatVar +       #spatial variation on DHI (productivity, buffer of 10km)
#   theta4[sp_id] * tempVar         #temporal variation on DHI (productivity, buffer of 10km)

#######################


# load packages -----------------------------------------------------------

library(ggplot2)
library(MCMCvis)
library(tidyverse)
library(tidyr)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
model_run_date <- '2025-01-21'
cv_data_date <- '2025-01-08'
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------

#cv data
cv_data_covs <- readRDS(paste0(dir, 'Data/L2/cv_data_covs-', cv_data_date, '.rds'))


#mass
fitMass <- readRDS(paste0(dir, 'Results/fitMass-within-pops-CV-', model_run_date, '/fitMass-fit-', model_run_date, '.rds'))
DATA_mass <- readRDS(paste0(dir, 'Results/fitMass-within-pops-CV-', model_run_date, '/fitMass-data-', model_run_date, '.rds'))
raw_data_mass <- readRDS(paste0(dir, 'Results/fitMass-within-pops-CV-', model_run_date, '/raw-mass-data-', model_run_date, '.rds'))


#wing
fitWing <- readRDS(paste0(dir, 'Results/fitWing-within-pops-CV-', model_run_date, '/fitWing-fit-', model_run_date, '.rds'))
DATA_wing <- readRDS(paste0(dir, 'Results/fitWing-within-pops-CV-', model_run_date, '/fitWing-data-', model_run_date, '.rds'))
raw_data_wing <- readRDS(paste0(dir, 'Results/fitWing-within-pops-CV-', model_run_date, '/raw-wing-data-', model_run_date, '.rds'))



#check result summary

#mass
MCMCsummary(fitMass, c("mu_gamma", "mu_theta1", "mu_theta2", "mu_theta3", "mu_theta4"), pg0=T)

#                 mean         sd        2.5%        50%       97.5% Rhat n.eff_bulk n.eff_tail  p>0
# mu_gamma  23.8956408 0.58933525 22.70099250 23.9016000 25.0430375 1.02        449        556 1.00
# mu_theta1 -0.1469941 0.02473177 -0.19682007 -0.1466925 -0.1000129 1.00       2587       4248 0.00
# mu_theta2  0.1497571 0.06404916  0.02426606  0.1491390  0.2806598 1.00       3235       4865 0.99
# mu_theta3  0.2841467 0.16859222 -0.06322086  0.2889830  0.5964668 1.00       2512       4008 0.95
# mu_theta4  0.1826092 0.26352617 -0.32843312  0.1793040  0.7066724 1.00       4158       5430 0.76

# pg0=T -> proportion of the posterior that is greater than 0


#wing
MCMCsummary(fitWing, c("mu_gamma", "mu_theta1", "mu_theta2", "mu_theta3", "mu_theta4"), pg0=T)

#                 mean          sd        2.5%          50%        97.5% Rhat n.eff_bulk n.eff_tail  p>0
# mu_gamma   7.077562209 0.075930550  6.92635075  7.076860000  7.224942750 1.01        540       1017 1.00
# mu_theta1 -0.011928943 0.002778788 -0.01720780 -0.011962350 -0.006342587 1.00       6364       7034 0.00
# mu_theta2  0.002197474 0.007769763 -0.01376391  0.002425705  0.016988972 1.00       6131       5149 0.63
# mu_theta3 -0.028829529 0.019803350 -0.06771532 -0.028863600  0.009908435 1.00       9671       7759 0.07
# mu_theta4  0.068297093 0.028022307  0.01397861  0.068324800  0.124240000 1.00       7617       7222 0.99


# Plot caterpillar plots:

# MASS

#Extract means for each species and parameter
theta1meanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "theta1")) %>% 
  dplyr::summarise(across(everything(), mean))

theta2meanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "theta2")) %>% 
  dplyr::summarise(across(everything(), mean))

theta3meanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "theta3")) %>% 
  dplyr::summarise(across(everything(), mean))

theta4meanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "theta4")) %>% 
  dplyr::summarise(across(everything(), mean))

species_meansMass <- as.data.frame(cbind(sci_name = unique(raw_data_mass$sci_name), 
                                         theta1meanM = t(theta1meanM),
                                         theta2meanM = t(theta2meanM), 
                                         theta3meanM = t(theta3meanM),
                                         theta4meanM = t(theta4meanM)), 
                                   row.names = FALSE) %>%
  dplyr::rename_with(~ c("sci_name", "theta1meanM", "theta2meanM", "theta3meanM", "theta4meanM"))


mu_theta1M <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "mu_theta1")) %>% 
  dplyr::summarise(mu = mean(mu_theta1),
                   CI055 = quantile(mu_theta1, probs=0.055),
                   CI945 = quantile(mu_theta1, probs=0.945))

mu_theta2M <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "mu_theta2")) %>% 
  dplyr::summarise(mu = mean(mu_theta2),
                   CI055 = quantile(mu_theta2, probs=0.055),
                   CI945 = quantile(mu_theta2, probs=0.945))

mu_theta3M <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "mu_theta3")) %>% 
  dplyr::summarise(mu = mean(mu_theta3),
                   CI055 = quantile(mu_theta3, probs=0.055),
                   CI945 = quantile(mu_theta3, probs=0.945))

mu_theta4M <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = "mu_theta4")) %>% 
  dplyr::summarise(mu = mean(mu_theta4),
                   CI055 = quantile(mu_theta4, probs=0.055),
                   CI945 = quantile(mu_theta4, probs=0.945))

overall_meanMass <- as.data.frame(rbind(mu_theta1M, mu_theta2M, mu_theta3M, mu_theta4M))
colnames(overall_meanMass) <- c("mean", "q055", "q945")
overall_meanMass$parameter <- factor(c("mu_theta1", "mu_theta2", "mu_theta3", "mu_theta4"))

# Reshape species means to long format for easier plotting
species_means_Masslong <- species_meansMass %>%
  tidyr::pivot_longer(cols = starts_with("theta"), names_to = "parameter", values_to = "value")

# Map parameters in species_means to the same structure as mu_thetas
species_means_Masslong$parameter <- factor(dplyr::recode(species_means_Masslong$parameter, 
                                                         "theta1meanM" = "mu_theta1",
                                                         "theta2meanM" = "mu_theta2",
                                                         "theta3meanM" = "mu_theta3",
                                                         "theta4meanM" = "mu_theta4"))


#pdf(paste0(dir, 'Results/plotAllresultsMass-', run_date, '.pdf'), width = 8, height = 8)

# Plot
mass <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 1) +
  
  geom_jitter(data = species_means_Masslong, aes(x = factor(parameter, 
                                                            levels = c("mu_theta4", "mu_theta3", "mu_theta2", "mu_theta1")), 
                                                            y = as.numeric(value)), 
             color = "#482677FF", alpha = 0.15, size = 2,  stroke = 0, width = 0.18) +
  
  geom_pointrange(data = overall_meanMass, aes(x = factor(parameter, 
                                                          levels = c("mu_theta4", "mu_theta3", "mu_theta2", "mu_theta1")),
                                               y = mean, ymin = q055, ymax = q945), 
                  color = "#482677FF", linewidth = 1.2, fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = "Estimated effect", title = "Body mass") +
  scale_x_discrete(labels = c("Temporal variation\n of productivity",
                              "Spatial variation\n of productivity",
                              "Distance to\n range edge",
                              "Latitude"),
                   expand = c(0.2,0.2)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, margin = margin(t = 15)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +

  coord_flip() 

#dev.off()



# WING

#Extract means for each species and parameter
theta1meanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "theta1")) %>% 
  dplyr::summarise(across(everything(), mean))

theta2meanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "theta2")) %>% 
  dplyr::summarise(across(everything(), mean))

theta3meanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "theta3")) %>% 
  dplyr::summarise(across(everything(), mean))

theta4meanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "theta4")) %>% 
  dplyr::summarise(across(everything(), mean))

species_meansWing <- as.data.frame(cbind(sci_name = unique(raw_data_wing$sci_name), 
                                         theta1meanW = t(theta1meanW),
                                         theta2meanW = t(theta2meanW), 
                                         theta3meanW = t(theta3meanW),
                                         theta4meanW = t(theta4meanW)), 
                                   row.names = FALSE) %>%
  dplyr::rename_with(~ c("sci_name", "theta1meanW", "theta2meanW", "theta3meanW", "theta4meanW"))


mu_theta1W <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "mu_theta1")) %>% 
  dplyr::summarise(mu = mean(mu_theta1),
                   CI055 = quantile(mu_theta1, probs=0.055),
                   CI945 = quantile(mu_theta1, probs=0.945))

mu_theta2W <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "mu_theta2")) %>% 
  dplyr::summarise(mu = mean(mu_theta2),
                   CI055 = quantile(mu_theta2, probs=0.055),
                   CI945 = quantile(mu_theta2, probs=0.945))

mu_theta3W <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "mu_theta3")) %>% 
  dplyr::summarise(mu = mean(mu_theta3),
                   CI055 = quantile(mu_theta3, probs=0.055),
                   CI945 = quantile(mu_theta3, probs=0.945))

mu_theta4W <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = "mu_theta4")) %>% 
  dplyr::summarise(mu = mean(mu_theta4),
                   CI055 = quantile(mu_theta4, probs=0.055),
                   CI945 = quantile(mu_theta4, probs=0.945))

overall_meanWing <- as.data.frame(rbind(mu_theta1W, mu_theta2W, mu_theta3W, mu_theta4W))
colnames(overall_meanWing) <- c("mean", "q055", "q945")
overall_meanWing$parameter <- factor(c("mu_theta1", "mu_theta2", "mu_theta3", "mu_theta4"))

# Reshape species means to long format for easier plotting
species_means_Winglong <- species_meansWing %>%
  tidyr::pivot_longer(cols = starts_with("theta"), names_to = "parameter", values_to = "value")

# Map parameters in species_means to the same structure as mu_thetas
species_means_Winglong$parameter <- factor(dplyr::recode(species_means_Winglong$parameter, 
                                                         "theta1meanW" = "mu_theta1",
                                                         "theta2meanW" = "mu_theta2",
                                                         "theta3meanW" = "mu_theta3",
                                                         "theta4meanW" = "mu_theta4"))


#pdf(paste0(dir, 'Results/plotAllresultsWing-', run_date, '.pdf'), width = 8, height = 8)

# Plot
wing <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 1) +
  geom_jitter(data = species_means_Winglong, aes(x = factor(parameter, 
                                                           levels = c("mu_theta4", "mu_theta3", "mu_theta2", "mu_theta1")), 
                                                y = as.numeric(value)), 
             color = "#29AF7FFF", alpha = 0.15, size = 2,  stroke = 0, width = 0.18) +
  geom_pointrange(data = overall_meanWing, aes(x = factor(parameter, 
                                                          levels = c("mu_theta4", "mu_theta3", "mu_theta2", "mu_theta1")),
                                               y = mean, ymin = q055, ymax = q945), 
                  color = "#29AF7FFF", linewidth = 1.2, fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = "Estimated effect", title = "Wing length") +
  scale_x_discrete(labels = c("Temporal variability\n of productivity",
                              "Spatial variability\n of productivity",
                              "Distance to\n range edge",
                              "Latitude"),
                   expand = c(0.2,0.2)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, margin = margin(t = 15)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

#dev.off()


pdf(paste0('Results/Fig-results-within-pops-', run_date, '.pdf'), width = 12, height = 6)
gridExtra::grid.arrange(mass, wing, nrow=1)
dev.off()




# combination line and cat plots ------------------------------------------

#extract posteriors - mass
mass_gamma_mn <- MCMCvis::MCMCpstr(fitMass, params = 'gamma')[[1]]
mass_theta1_mn <- MCMCvis::MCMCpstr(fitMass, params = 'theta1')[[1]]
mass_theta2_mn <- MCMCvis::MCMCpstr(fitMass, params = 'theta2')[[1]]
mass_theta3_mn <- MCMCvis::MCMCpstr(fitMass, params = 'theta3')[[1]]
mass_theta4_mn <- MCMCvis::MCMCpstr(fitMass, params = 'theta4')[[1]]

mass_mu_gamma_ch <- MCMCvis::MCMCchains(fitMass, params = 'mu_gamma')
mass_mu_theta1_ch <- MCMCvis::MCMCchains(fitMass, params = 'mu_theta1')
mass_mu_theta2_ch <- MCMCvis::MCMCchains(fitMass, params = 'mu_theta2')
mass_mu_theta3_ch <- MCMCvis::MCMCchains(fitMass, params = 'mu_theta3')
mass_mu_theta4_ch <- MCMCvis::MCMCchains(fitMass, params = 'mu_theta4')

#extract posteriors - wing
wing_gamma_mn <- MCMCvis::MCMCpstr(fitWing, params = 'gamma')[[1]]
wing_theta1_mn <- MCMCvis::MCMCpstr(fitWing, params = 'theta1')[[1]]
wing_theta2_mn <- MCMCvis::MCMCpstr(fitWing, params = 'theta2')[[1]]
wing_theta3_mn <- MCMCvis::MCMCpstr(fitWing, params = 'theta3')[[1]]
wing_theta4_mn <- MCMCvis::MCMCpstr(fitWing, params = 'theta4')[[1]]

wing_mu_gamma_ch <- MCMCvis::MCMCchains(fitWing, params = 'mu_gamma')
wing_mu_theta1_ch <- MCMCvis::MCMCchains(fitWing, params = 'mu_theta1')
wing_mu_theta2_ch <- MCMCvis::MCMCchains(fitWing, params = 'mu_theta2')
wing_mu_theta3_ch <- MCMCvis::MCMCchains(fitWing, params = 'mu_theta3')
wing_mu_theta4_ch <- MCMCvis::MCMCchains(fitWing, params = 'mu_theta4')


#for each species
str(raw_data_mass)
usp <- unique(raw_data_mass$sp_id)

#sim values mass
pred_df <- data.frame(sp_id = rep(NA, 50 * length(usp)),
                      lat = NA,
                      lat_sc = NA,
                      mass_pred_lat = NA,
                      wing_pred_lat = NA,
                      ldre = NA,
                      ldre_sc = NA,
                      mass_pred_ldre = NA,
                      wing_pred_ldre = NA,
                      DHI_spat = NA,
                      DHI_spat_sc = NA,
                      mass_pred_DHI_spat = NA,
                      wing_pred_DHI_spat = NA,
                      DHI_temp = NA,
                      DHI_temp_sc = NA,
                      mass_pred_DHI_temp = NA,
                      wing_pred_DHI_temp = NA)
counter <- 1
for (i in 1:length(usp))
{
  #i <- 1
  tt <- dplyr::filter(raw_data_mass, sp_id == usp[i]) %>%
    dplyr::mutate(DHI_spat = log(DHI_spatial_var_10km), 
                  DHI_temp = log(DHI_temporal_var_10km)) %>%
    dplyr::select(lat, lat_sc, 
                  ldre, ldre_sc, 
                  DHI_spat,
                  DHI_spat_sc, 
                  DHI_temp,
                  DHI_temp_sc)
  
  #LAT
  #get sc values and sim line
  t_lat_sc_sim <- seq(min(tt$lat_sc), max(tt$lat_sc), length.out = 50)
  #to rescale because data were * 1000 to fit model
  t_lat_mass_pred <- (mass_gamma_mn[i] + mass_theta1_mn[i] * t_lat_sc_sim) / 1000
  t_lat_wing_pred <- (wing_gamma_mn[i] + wing_theta1_mn[i] * t_lat_sc_sim) / 1000
  
  #LDRE
  t_ldre_sc_sim <- seq(min(tt$ldre_sc), max(tt$ldre_sc), length.out = 50)
  #to rescale because data were * 1000 to fit model
  t_ldre_mass_pred <- (mass_gamma_mn[i] + mass_theta2_mn[i] * t_ldre_sc_sim) / 1000
  t_ldre_wing_pred <- (wing_gamma_mn[i] + wing_theta2_mn[i] * t_ldre_sc_sim) / 1000
  
  #SPAT
  t_DHI_spat_sc_sim <- seq(min(tt$DHI_spat_sc), max(tt$DHI_spat_sc), length.out = 50)
  #to rescale because data were * 1000 to fit model
  t_DHI_spat_mass_pred <- (mass_gamma_mn[i] + mass_theta3_mn[i] * t_DHI_spat_sc_sim) / 1000
  t_DHI_spat_wing_pred <- (wing_gamma_mn[i] + wing_theta3_mn[i] * t_DHI_spat_sc_sim) / 1000
  
  #TEMP
  t_DHI_temp_sc_sim <- seq(min(tt$DHI_temp_sc), max(tt$DHI_temp_sc), length.out = 50)
  #to rescale because data were * 1000 to fit model
  t_DHI_temp_mass_pred <- (mass_gamma_mn[i] + mass_theta4_mn[i] * t_DHI_temp_sc_sim) / 1000
  t_DHI_temp_wing_pred <- (wing_gamma_mn[i] + wing_theta4_mn[i] * t_DHI_temp_sc_sim) / 1000
  
  #fill df cov
  pred_df$sp_id[counter:(counter + 49)] <- i
  #transform to raw lat
  pred_df$lat[counter:(counter + 49)] <- t_lat_sc_sim + mean(tt$lat)
  pred_df$lat_sc[counter:(counter + 49)] <- t_lat_sc_sim
  
  pred_df$ldre[counter:(counter + 49)] <- t_ldre_sc_sim + mean(tt$ldre)
  pred_df$ldre_sc[counter:(counter + 49)] <- t_ldre_sc_sim
  
  pred_df$DHI_spat[counter:(counter + 49)] <- t_DHI_spat_sc_sim + mean(tt$DHI_spat)
  pred_df$DHI_spat_sc[counter:(counter + 49)] <- t_DHI_spat_sc_sim
  
  pred_df$DHI_temp[counter:(counter + 49)] <- t_DHI_temp_sc_sim + mean(tt$DHI_temp)
  pred_df$DHI_temp_sc[counter:(counter + 49)] <- t_DHI_temp_sc_sim
  
  #fill df ped
  pred_df$mass_pred_lat[counter:(counter + 49)] <- t_lat_mass_pred
  pred_df$wing_pred_lat[counter:(counter + 49)] <- t_lat_wing_pred
  
  pred_df$mass_pred_ldre[counter:(counter + 49)] <- t_ldre_mass_pred
  pred_df$wing_pred_ldre[counter:(counter + 49)] <- t_ldre_wing_pred
  
  pred_df$mass_pred_DHI_spat[counter:(counter + 49)] <- t_DHI_spat_mass_pred
  pred_df$wing_pred_DHI_spat[counter:(counter + 49)] <- t_DHI_spat_wing_pred
  
  pred_df$mass_pred_DHI_temp[counter:(counter + 49)] <- t_DHI_temp_mass_pred
  pred_df$wing_pred_DHI_temp[counter:(counter + 49)] <- t_DHI_temp_wing_pred
  counter <- counter + 50
}

#get lat values and sim line
o_lat_sim <- seq(min(pred_df$lat), max(pred_df$lat), length.out = 50)
#get lat_sc values
o_lat_sc_sim <- o_lat_sim - mean(pred_df$lat)
#to rescale because data were * 1000 to fit model
o_lat_mass_pred_mat <- matrix(NA, 
                              nrow = NROW(mass_mu_gamma_ch), 
                              ncol = 50)
o_lat_wing_pred_mat <- matrix(NA, 
                              nrow = NROW(mass_mu_gamma_ch), 
                              ncol = 50)
for (i in 1:NROW(mass_mu_gamma_ch))
{
  o_lat_mass_pred_mat[i,] <- (mass_mu_gamma_ch[i,] + 
                                mass_mu_theta1_ch[i,] * o_lat_sc_sim) / 1000  
  o_lat_wing_pred_mat[i,] <- (wing_mu_gamma_ch[i,] + 
                                wing_mu_theta1_ch[i,] * o_lat_sc_sim) / 1000  
}

#ldre
o_ldre_sim <- seq(min(pred_df$ldre), max(pred_df$ldre), length.out = 50)
o_ldre_sc_sim <- o_ldre_sim - mean(pred_df$ldre)
o_ldre_mass_pred_mat <- matrix(NA, 
                               nrow = NROW(mass_mu_gamma_ch), 
                               ncol = 50)
o_ldre_wing_pred_mat <- matrix(NA, 
                               nrow = NROW(mass_mu_gamma_ch), 
                               ncol = 50)
for (i in 1:NROW(mass_mu_gamma_ch))
{
  o_ldre_mass_pred_mat[i,] <- (mass_mu_gamma_ch[i,] + 
                                 mass_mu_theta2_ch[i,] * o_ldre_sc_sim) / 1000  
  o_ldre_wing_pred_mat[i,] <- (wing_mu_gamma_ch[i,] + 
                                 wing_mu_theta2_ch[i,] * o_ldre_sc_sim) / 1000  
}


#SPAT
o_spat_sim <- seq(min(pred_df$DHI_spat), max(pred_df$DHI_spat), length.out = 50)
o_spat_sc_sim <- o_spat_sim - mean(pred_df$DHI_spat)
o_DHI_spat_mass_pred_mat <- matrix(NA, 
                                   nrow = NROW(mass_mu_gamma_ch), 
                                   ncol = 50)
o_DHI_spat_wing_pred_mat <- matrix(NA, 
                                   nrow = NROW(mass_mu_gamma_ch), 
                                   ncol = 50)
for (i in 1:NROW(mass_mu_gamma_ch))
{
  o_DHI_spat_mass_pred_mat[i,] <- (mass_mu_gamma_ch[i,] + 
                                     mass_mu_theta3_ch[i,] * o_spat_sc_sim) / 1000  
  o_DHI_spat_wing_pred_mat[i,] <- (wing_mu_gamma_ch[i,] + 
                                     wing_mu_theta3_ch[i,] * o_spat_sc_sim) / 1000  
}



#TEMP
o_temp_sim <- seq(min(pred_df$DHI_temp), max(pred_df$DHI_temp), length.out = 50)
o_temp_sc_sim <- o_temp_sim - mean(pred_df$DHI_temp)
o_DHI_temp_mass_pred_mat <- matrix(NA, 
                                   nrow = NROW(mass_mu_gamma_ch), 
                                   ncol = 50)
o_DHI_temp_wing_pred_mat <- matrix(NA, 
                                   nrow = NROW(mass_mu_gamma_ch), 
                                   ncol = 50)
for (i in 1:NROW(mass_mu_gamma_ch))
{
  o_DHI_temp_mass_pred_mat[i,] <- (mass_mu_gamma_ch[i,] + 
                                     mass_mu_theta4_ch[i,] * o_temp_sc_sim) / 1000  
  o_DHI_temp_wing_pred_mat[i,] <- (wing_mu_gamma_ch[i,] + 
                                     wing_mu_theta4_ch[i,] * o_temp_sc_sim) / 1000  
}


#plot df for cross species effects
pred_df2 <- data.frame(lat = o_lat_sim,
                       lat_sc = o_lat_sc_sim,
                       mass_pred_lat = apply(o_lat_mass_pred_mat, 2, mean),
                       mass_pred_lat_LCI = apply(o_lat_mass_pred_mat, 2, 
                                                 function(x) quantile(x, probs = 0.055)),
                       mass_pred_lat_UCI = apply(o_lat_mass_pred_mat, 2, 
                                                 function(x) quantile(x, probs = 0.945)),
                       wing_pred_lat = apply(o_lat_wing_pred_mat, 2, mean),
                       wing_pred_lat_LCI = apply(o_lat_wing_pred_mat, 2, 
                                                 function(x) quantile(x, probs = 0.055)),
                       wing_pred_lat_UCI = apply(o_lat_wing_pred_mat, 2, 
                                                 function(x) quantile(x, probs = 0.945)),
                       ldre = o_ldre_sim,
                       ldre_sc = o_ldre_sc_sim,
                       mass_pred_ldre = o_ldre_mass_pred,
                       wing_pred_ldre = o_ldre_wing_pred,
                       mass_pred_ldre = apply(o_ldre_mass_pred_mat, 2, mean),
                       mass_pred_ldre_LCI = apply(o_ldre_mass_pred_mat, 2, 
                                                  function(x) quantile(x, probs = 0.055)),
                       mass_pred_ldre_UCI = apply(o_ldre_mass_pred_mat, 2, 
                                                  function(x) quantile(x, probs = 0.945)),
                       wing_pred_ldre = apply(o_ldre_wing_pred_mat, 2, mean),
                       wing_pred_ldre_LCI = apply(o_ldre_wing_pred_mat, 2, 
                                                  function(x) quantile(x, probs = 0.055)),
                       wing_pred_ldre_UCI = apply(o_ldre_wing_pred_mat, 2, 
                                                  function(x) quantile(x, probs = 0.945)),
                       DHI_spat = o_spat_sim,
                       DHI_spat_sc = o_spat_sc_sim,
                       mass_pred_DHI_spat = apply(o_DHI_spat_mass_pred_mat, 2, mean),
                       mass_pred_DHI_spat_LCI = apply(o_DHI_spat_mass_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.055)),
                       mass_pred_DHI_spat_UCI = apply(o_DHI_spat_mass_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.945)),
                       wing_pred_DHI_spat = apply(o_DHI_spat_wing_pred_mat, 2, mean),
                       wing_pred_DHI_spat_LCI = apply(o_DHI_spat_wing_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.055)),
                       wing_pred_DHI_spat_UCI = apply(o_DHI_spat_wing_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.945)),
                       DHI_temp = o_temp_sim,
                       DHI_temp_sc = o_temp_sc_sim,
                       mass_pred_DHI_temp = apply(o_DHI_temp_mass_pred_mat, 2, mean),
                       mass_pred_DHI_temp_LCI = apply(o_DHI_temp_mass_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.055)),
                       mass_pred_DHI_temp_UCI = apply(o_DHI_temp_mass_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.945)),
                       wing_pred_DHI_temp = apply(o_DHI_temp_wing_pred_mat, 2, mean),
                       wing_pred_DHI_temp_LCI = apply(o_DHI_temp_wing_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.055)),
                       wing_pred_DHI_temp_UCI = apply(o_DHI_temp_wing_pred_mat, 2, 
                                                      function(x) quantile(x, probs = 0.945)))


#get ylim
# str(pred_df2)
# ps <- dplyr::select(pred_df2,
#               -lat, -lat_sc,
#               -ldre, -ldre_sc,
#               -DHI_spat, -DHI_spat_sc,
#               -DHI_temp, - DHI_temp_sc) 
# rng_mass <- range(dplyr::select(ps, dplyr::starts_with('mass'))) 
# rng_wing <- range(dplyr::select(ps, dplyr::starts_with('wing')))

#LAT
lat_mass_plt <- ggplot() + 
  #mass predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = lat, y = mass_pred_lat, group = sp_id), 
  #           color = "#482677FF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  #mass predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = lat, ymin = mass_pred_lat_LCI, ymax = mass_pred_lat_UCI), 
              fill = "#482677FF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = lat, y = mass_pred_lat), 
            color = "#482677FF", 
            lwd = 2) +
  # #wing predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = lat, y = wing_pred_lat, group = sp_id), 
  #           color = "#29AF7FFF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  # #wing predictions - cross-species
  # geom_ribbon(data = pred_df2, 
  #           aes(x = lat, ymin = wing_pred_lat_LCI, ymax = wing_pred_lat_UCI), 
  #           fill = "#29AF7FFF", alpha = 0.2) +
  # geom_line(data = pred_df2, 
  #           aes(x = lat, y = wing_pred_lat), 
  #           color = "#29AF7FFF", 
  #           lwd = 2) +
  theme_bw() + 
  xlab('Latitude') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_y_continuous(breaks = c(0.020, 0.023, 0.026),
                     limits = rng_mass) +
  ggtitle('lat_mass')

lat_wing_plt <- ggplot() + 
  #wing predictions - species-specific
  geom_ribbon(data = pred_df2,
              aes(x = lat, ymin = wing_pred_lat_LCI, ymax = wing_pred_lat_UCI),
              fill = "#29AF7FFF", alpha = 0.2) +
  geom_line(data = pred_df2,
            aes(x = lat, y = wing_pred_lat),
            color = "#29AF7FFF",
            lwd = 2) +
  theme_bw() + 
  xlab('Latitude') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_y_continuous(breaks = c(0.0066, 0.0070, 0.0074),
                     limits = rng_wing) +
  ggtitle('lat wing')


#LDRE
ldre_mass_plt <- ggplot() + 
  # #mass predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = ldre, y = mass_pred_ldre, group = sp_id), 
  #           color = "#482677FF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  #mass predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = ldre, ymin = mass_pred_ldre_LCI, ymax = mass_pred_ldre_UCI), 
              fill = "#482677FF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = ldre, y = mass_pred_ldre), 
            color = "#482677FF", 
            lwd = 2) +
  # #wing predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = ldre, y = wing_pred_ldre, group = sp_id), 
  #           color = "#29AF7FFF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  # #wing predictions - cross-species
  # geom_line(data = pred_df2, 
  #           aes(x = ldre, y = wing_pred_ldre), 
  #           color = "#29AF7FFF", 
  #           lwd = 2) +
  theme_bw() + 
  xlab('LDRE') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0, log(7), log(55), log(400)), 
                     labels = c(0, 7, 55, 400)) +
  scale_y_continuous(breaks = c(0.020, 0.023, 0.026),
                     limits = rng_mass) +
  ggtitle('ldre mass')


ldre_wing_plt <- ggplot() + 
  #wing predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = ldre, ymin = wing_pred_ldre_LCI, ymax = wing_pred_ldre_UCI), 
              fill = "#29AF7FFF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = ldre, y = wing_pred_ldre), 
            color = "#29AF7FFF", 
            lwd = 2) +
  theme_bw() + 
  xlab('LDRE') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(0, log(7), log(55), log(400)), 
                     labels = c(0, 7, 55, 400)) +
  scale_y_continuous(breaks = c(0.0066, 0.0070, 0.0074),
                     limits = rng_wing) +
  ggtitle('ldre wing')



#SPAT
spat_mass_plt <- ggplot() + 
  #mass predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = DHI_spat, y = mass_pred_DHI_spat, group = sp_id), 
  #           color = "#482677FF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  #mass predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = DHI_spat, ymin = mass_pred_DHI_spat_LCI, ymax = mass_pred_DHI_spat_UCI), 
              fill = "#482677FF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = DHI_spat, y = mass_pred_DHI_spat), 
            color = "#482677FF", 
            lwd = 2) +
  # #wing predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = DHI_spat, y = wing_pred_DHI_spat, group = sp_id), 
  #           color = "#29AF7FFF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  # #wing predictions - cross-species
  # geom_line(data = pred_df2, 
  #           aes(x = DHI_spat, y = wing_pred_DHI_spat), 
  #           color = "#29AF7FFF", 
  #           lwd = 2) +
  theme_bw() + 
  xlab('SPAT') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(log(0.02), log(0.09), log(0.37)), 
                     labels = c(0.02, 0.09, 0.37)) +
  scale_y_continuous(breaks = c(0.020, 0.023, 0.026),
                     limits = rng_mass) +
  ggtitle('spat mass')

spat_wing_plt <- ggplot() + 
  #wing predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = DHI_spat, ymin = wing_pred_DHI_spat_LCI, ymax = wing_pred_DHI_spat_UCI), 
              fill = "#29AF7FFF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = DHI_spat, y = wing_pred_DHI_spat), 
            color = "#29AF7FFF", 
            lwd = 2) +
  theme_bw() + 
  xlab('SPAT') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(log(0.02), log(0.09), log(0.37)), 
                     labels = c(0.02, 0.09, 0.37)) +
  scale_y_continuous(breaks = c(0.0066, 0.0070, 0.0074),
                     limits = rng_wing) +
  ggtitle('spat wing')


#TEMP
temp_mass_plt <- ggplot() + 
  #mass predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = DHI_temp, y = mass_pred_DHI_temp, group = sp_id), 
  #           color = "#482677FF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  #mass predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = DHI_temp, ymin = mass_pred_DHI_temp_LCI, ymax = mass_pred_DHI_temp_UCI), 
              fill = "#482677FF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = DHI_temp, y = mass_pred_DHI_temp), 
            color = "#482677FF", 
            lwd = 2) +
  # #wing predictions - species-specific
  # geom_line(data = pred_df, 
  #           aes(x = DHI_temp, y = wing_pred_DHI_temp, group = sp_id), 
  #           color = "#29AF7FFF", 
  #           lwd = 1.2,
  #           alpha = 0.2) + 
  # #wing predictions - cross-species
  # geom_line(data = pred_df2, 
  #           aes(x = DHI_temp, y = wing_pred_DHI_temp), 
  #           color = "#29AF7FFF", 
  #           lwd = 2) +
  theme_bw() + 
  xlab('TEMP') + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(log(0.02), log(0.09), log(0.37)), 
                     labels = c(0.02, 0.09, 0.37)) +
  scale_y_continuous(breaks = c(0.020, 0.023, 0.026),
                     limits = rng_mass) +
  ggtitle('temp mass')

temp_wing_plt <- ggplot() + 
  #wing predictions - cross-species
  geom_ribbon(data = pred_df2, 
              aes(x = DHI_temp, ymin = wing_pred_DHI_temp_LCI, ymax = wing_pred_DHI_temp_UCI), 
              fill = "#29AF7FFF", alpha = 0.2) +
  geom_line(data = pred_df2, 
            aes(x = DHI_temp, y = wing_pred_DHI_temp), 
            color = "#29AF7FFF", 
            lwd = 2) +
  theme_bw() + 
  xlab('TEMP') + 
  ylab(NULL) +
  ylim(rng_wing) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_continuous(breaks = c(log(0.02), log(0.09), log(0.37)), 
                     labels = c(0.02, 0.09, 0.37)) +
  scale_y_continuous(breaks = c(0.0066, 0.0070, 0.0074),
                     limits = rng_wing) +
  ggtitle('temp wing')

H <- 1.8
W <- 2

ggsave(filename = paste0(dir, 'Results/lat_mass.pdf'),
       width = W,
       height = H,
       lat_mass_plt)

ggsave(filename = paste0(dir, 'Results/lat_wing.pdf'),
       width = W,
       height = H,
       lat_wing_plt)

ggsave(filename = paste0(dir, 'Results/ldre_mass.pdf'),
       width = W,
       height = H,
       ldre_mass_plt)

ggsave(filename = paste0(dir, 'Results/ldre_wing.pdf'),
       width = W,
       height = H,
       ldre_wing_plt)

ggsave(filename = paste0(dir, 'Results/spat_mass.pdf'),
       width = W,
       height = H,
       spat_mass_plt)

ggsave(filename = paste0(dir, 'Results/spat_wing.pdf'),
       width = W,
       height = H,
       spat_wing_plt)

ggsave(filename = paste0(dir, 'Results/temp_mass.pdf'),
       width = W,
       height = H,
       temp_mass_plt)

ggsave(filename = paste0(dir, 'Results/temp_wing.pdf'),
       width = W,
       height = H,
       temp_wing_plt)


#cat plots
#LAT
lat_mass_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Masslong,
                                   parameter == 'mu_theta1'),
              aes(x = factor(parameter,
                             levels = c("mu_theta1")),
                  y = as.numeric(value)),
              color = "#482677FF",
              alpha = 0.10,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanMass,
                                       parameter == 'mu_theta1'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta1")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#482677FF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

lat_wing_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Winglong,
                                   parameter == 'mu_theta1'),
              aes(x = factor(parameter,
                             levels = c("mu_theta1")),
                  y = as.numeric(value)),
              color = "#29AF7FFF",
              alpha = 0.15,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanWing,
                                       parameter == 'mu_theta1'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta1")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#29AF7FFF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

# LDRE
ldre_mass_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Masslong,
                                   parameter == 'mu_theta2'),
              aes(x = factor(parameter,
                             levels = c("mu_theta2")),
                  y = as.numeric(value)),
              color = "#482677FF",
              alpha = 0.10,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanMass,
                                       parameter == 'mu_theta2'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta2")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#482677FF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

ldre_wing_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Winglong,
                                   parameter == 'mu_theta2'),
              aes(x = factor(parameter,
                             levels = c("mu_theta2")),
                  y = as.numeric(value)),
              color = "#29AF7FFF",
              alpha = 0.15,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanWing,
                                       parameter == 'mu_theta2'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta2")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#29AF7FFF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

#SPAT
spat_mass_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Masslong,
                                   parameter == 'mu_theta3'),
              aes(x = factor(parameter,
                             levels = c("mu_theta3")),
                  y = as.numeric(value)),
              color = "#482677FF",
              alpha = 0.10,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanMass,
                                       parameter == 'mu_theta3'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta3")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#482677FF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

spat_wing_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Winglong,
                                   parameter == 'mu_theta3'),
              aes(x = factor(parameter,
                             levels = c("mu_theta3")),
                  y = as.numeric(value)),
              color = "#29AF7FFF",
              alpha = 0.15,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanWing,
                                       parameter == 'mu_theta3'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta3")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#29AF7FFF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()


# TEMP
temp_mass_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Masslong,
                                   parameter == 'mu_theta4'),
              aes(x = factor(parameter,
                             levels = c("mu_theta4")),
                  y = as.numeric(value)),
              color = "#482677FF",
              alpha = 0.10,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanMass,
                                       parameter == 'mu_theta4'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta4")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#482677FF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

temp_wing_cat <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray60",
             linewidth = 1) +
  geom_jitter(data = dplyr::filter(species_means_Winglong,
                                   parameter == 'mu_theta4'),
              aes(x = factor(parameter,
                             levels = c("mu_theta4")),
                  y = as.numeric(value)),
              color = "#29AF7FFF",
              alpha = 0.15,
              size = 3,
              stroke = 0,
              width = 0.20) +
  geom_pointrange(data = dplyr::filter(overall_meanWing,
                                       parameter == 'mu_theta4'),
                  aes(x = factor(parameter,
                                 levels = c("mu_theta4")),
                      y = mean,
                      ymin = q055,
                      ymax = q945),
                  color = "#29AF7FFF",
                  linewidth = 1.2,
                  fatten = 5) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  coord_flip()

H2 <- 1
W2 <- 1.8

ggsave(filename = paste0(dir, 'Results/lat_mass_cat.pdf'),
       width = W2,
       height = H2,
       lat_mass_cat)

ggsave(filename = paste0(dir, 'Results/lat_wing_cat.pdf'),
       width = W2,
       height = H2,
       lat_wing_cat)


ggsave(filename = paste0(dir, 'Results/ldre_mass_cat.pdf'),
       width = W2,
       height = H2,
       ldre_mass_cat)

ggsave(filename = paste0(dir, 'Results/ldre_wing_cat.pdf'),
       width = W2,
       height = H2,
       ldre_wing_cat)


ggsave(filename = paste0(dir, 'Results/spat_mass_cat.pdf'),
       width = W2,
       height = H2,
       spat_mass_cat)

ggsave(filename = paste0(dir, 'Results/spat_wing_cat.pdf'),
       width = W2,
       height = H2,
       spat_wing_cat)


ggsave(filename = paste0(dir, 'Results/temp_mass_cat.pdf'),
       width = W2,
       height = H2,
       temp_mass_cat)

ggsave(filename = paste0(dir, 'Results/temp_wing_cat.pdf'),
       width = W2,
       height = H2,
       temp_wing_cat)


# Supplemental information --------------------------------------------


#Figure S1. PPC (multipanel - 2 panels)

# Plot the Posterior Predictive Check (PPC) --------------------------------------------

# PPC mass 
y_val <- DATA_mass$CVobs
y_rep <- MCMCvis::MCMCchains(fitMass, params = 'CVobs_rep')

pdf(paste0(dir, 'Results/', 'PPC-within-pops-', run_date, '.pdf'), height = 5, width = 8)

par(mfrow=c(1,2))

plot(density(y_val), lwd = 2, main = 'PPC - Mass CV within pops', xlab = 'Value') 
for (i in 1:100){
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.1))
}
# Add legend
legend("topright", 
       legend = c("CV observed", "CV predicted"), 
       col = c("black", rgb(1, 0, 0, 0.1)), 
       lwd = c(2, 2), 
       bty = "n")  # 'bty = "n"' no box around the legend


#PPC wing
y_val <- DATA_wing$CVobs
y_rep <- MCMCvis::MCMCchains(fitWing, params = 'CVobs_rep')

plot(density(y_val), lwd = 2, main = 'PPC - Wing CV within pops', xlab = 'Value', ylim =c(0, 0.4))
for (i in 1:100){
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.1))
}
# Add legend
legend("topright", 
       legend = c("CV observed", "CV predicted"), 
       col = c("black", rgb(1, 0, 0, 0.1)), 
       lwd = c(2, 2), 
       bty = "n")  # 'bty = "n"' no box around the legend

dev.off()



# Plot the species-specific effects --------------------------------------------

#Figure S3. caterpillar plots from model #1 for mass

spNamesLabel <- unique(gsub('_', ' ', raw_data_mass$sci_name))

parameters <- c('theta1', 'theta2', 'theta3', 'theta4')
parameter_name <- c('latitude', 
                    'distance to range edge', 
                    'spatial variation of productivity',
                    'temporal variation of productivity')
xlimits <- list(c(-0.8,0.4), c(-2,2), c(-2,2), c(-6,6))

pdf(paste0(dir, 'Results/CV-Mass-per-sp-', run_date, '.pdf'), width = 8, height = 15)

for(i in 1:length(parameters)){
  MCMCvis::MCMCplot(fitMass, params = parameters[i], 
                    main = bquote(atop(beta[.(i)*k], "Effect of" ~ .(parameter_name[i]) ~ "on mass CV")), 
                    ci = c(50, 89), labels = spNamesLabel, sz_labels=0.6, sz_med =1.1, guide_lines = TRUE,
                    xlim = xlimits[[i]])
}

dev.off()


#Figure S4. caterpillar plots from model #1 for wing

spNamesLabel <- unique(gsub('_', ' ', raw_data_wing$sci_name))

parameters <- c('theta1', 'theta2', 'theta3', 'theta4')
parameter_name <- c('latitude', 
                    'distance to range edge', 
                    'spatial variation of productivity',
                    'temporal variation of productivity')
xlimits <- list(c(-0.06,0.02), c(-0.1,0.1), c(-0.15,0.10), c(-0.2,0.4))

pdf(paste0(dir, 'Results/CV-Wing-per-sp-', run_date, '.pdf'), width = 8, height = 15)

for(i in 1:length(parameters)){
  MCMCvis::MCMCplot(fitWing, params = parameters[i], 
                    main = bquote(atop(beta[.(i)*k], "Effect of" ~ .(parameter_name[i]) ~ "on wing CV")), 
                    ci = c(50, 89), labels = spNamesLabel, sz_labels=0.6, sz_med =1.1, guide_lines = TRUE,
                    xlim = xlimits[[i]])
}

dev.off()




# Plot the Prior Posterior Overlap (PPO) --------------------------------------------

#PPO mass

#Simulate values for the priors.
#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mu_gamma = rnorm(10000, 25, 5),
                    mu_theta1 = rnorm(10000, 0, 1),
                    mu_theta2 = rnorm(10000, 0, 1),
                    mu_theta3 = rnorm(10000, 0, 5),
                    mu_theta4 = rnorm(10000, 0, 5),
                    sigma = rnorm(10000, 2, 1))

MCMCtrace(fitMass, 
          params = c('mu_gamma', 'mu_theta1', 'mu_theta2', 
                     'mu_theta3', 'mu_theta4', 'sigma'),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE,
          pdf = TRUE, 
          filename = paste0('PPO-mass-within-pops-', run_date, '.pdf'),
          wd = paste0(dir, 'Results'))


# PPO wing

#Simulate values for the priors.
#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mu_gamma = rnorm(10000, 10, 5),
                    mu_theta1 = rnorm(10000, 0, 1),
                    mu_theta2 = rnorm(10000, 0, 1),
                    mu_theta3 = rnorm(10000, 0, 1),
                    mu_theta4 = rnorm(10000, 0, 1),
                    sigma = rnorm(10000, 0, 1))

MCMCtrace(fitWing, 
          params = c('mu_gamma', 'mu_theta1', 'mu_theta2', 
                     'mu_theta3', 'mu_theta4', 'sigma'),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE,
          pdf = TRUE, 
          filename = paste0('PPO-wing-within-pops-', run_date, '.pdf'),
          wd = paste0(dir, 'Results'))


# Calculate the % change over latitude for each species --------------------------------------------

#extract parameter estimates
gamma_spM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'gamma'))
thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'theta1'))

gamma_spW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'gamma'))
thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'theta1'))

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
lat_data <- raw_data_mass %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(meanLat = mean(as.numeric(lat)),
                   sdLat = sd(as.numeric(lat)),
                   Latminim = min(as.numeric(lat)),
                   Latmaxim = max(as.numeric(lat)),
                   minLat_str = min(as.numeric(lat_sc)),
                   maxLat_str = max(as.numeric(lat_sc)),
                   MinMaxLatSp = n())

result_list <- list()
percentage_changeM <- percentage_changeW <- as.data.frame(matrix(NA, nrow=nrow(lat_data), ncol = 5))
colnames(percentage_changeM) <- colnames(percentage_changeW) <- c('sci_name', 'mean_change', 'q055', 'q945','trait') 
percentage_changeM[, 'sci_name'] <- percentage_changeW[, 'sci_name'] <- lat_data$sci_name
percentage_changeM[, 'trait'] <- rep('mass', nrow(lat_data))
percentage_changeW[, 'trait'] <- rep('wing', nrow(lat_data))


for(i in 1:nrow(lat_data)) {
  tmpLat <- seq(lat_data$minLat_str[i], lat_data$maxLat_str[i], length.out = lat_data$MinMaxLatSp[i])
  
  tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpLat)))
  tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpLat)))
  
  for(t in 1:length(tmpLat)){
    for(j in 1:nrow(gamma_spM)){
      tmpCVmass[j,t] <- (gamma_spM[j,i] + (thetameanM[j,i] * tmpLat[t]))/1000
      tmpCVwing[j,t] <- (gamma_spW[j,i] + (thetameanW[j,i] * tmpLat[t]))/1000
      
      # what I need is only the predictions of the min lat and max lat (first and last columns)
      # Combine tmpLat and tmpCV into a data frame
      tmp_df <- data.frame(tmpCVmass = tmpCVmass[,c(1, length(tmpLat))], 
                           tmpCVwing = tmpCVwing[,c(1, length(tmpLat))])
      
      result_list[[i]] <- tmp_df #saving all values in case we need to extract other estimates
    }
  }
 
  #(predicted CV at max lat - predicted CV at min lat) / predicted CV at min lat = % change
  #changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  #changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  # % change in mass and wing cv per 10 degrees latitude
  changeMass <- ((((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])/lat_data[i,6])*100)*10
  changeWing <- ((((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])/lat_data[i,6])*100)*10
  
  percentage_changeM[i, "mean_change"] <- mean(changeMass)
  percentage_changeM[i, "q055"] <- quantile(changeMass, probs = 0.055)
  percentage_changeM[i, "q945"] <- quantile(changeMass, probs = 0.945)
  
  percentage_changeW[i, "mean_change"] <- mean(changeWing)
  percentage_changeW[i, "q055"] <- quantile(changeWing, probs = 0.055)
  percentage_changeW[i, "q945"] <- quantile(changeWing, probs = 0.945)
  
}

mean(as.numeric(percentage_changeM$mean_change))
mean(as.numeric(percentage_changeW$mean_change))
                                 

#join data frames
percentage_change <- rbind(percentage_changeM, percentage_changeW)

pdf(paste0(dir, 'Results/perc-change10degree-', run_date, '.pdf'), width = 8, height = 18)

#plot results
ggplot() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 0.5) +
  geom_pointrange(data = percentage_change, aes(y = reorder(sci_name, desc(sci_name)),
                                               x = mean_change, xmin = q055, xmax = q945, group = trait, colour = trait), 
                  linewidth = 1.2, fatten = 2, position = position_dodge(width = 0.6)) +
  labs(x = '% Change on body mass and wing length CV per 10 degree latitude increase',
       y = NULL,
       color = NULL) + 
  scale_color_manual(values = c("mass" = "#482677FF", "wing" = "#29AF7FFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, face='italic'),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, margin = margin(t = 15)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color='gray90', size=0.3),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.9, 0.06),
        legend.text = element_text(size = 12),  # Increase font size for legend text
        legend.title = element_text(size = 12))

dev.off()


#write out
write.csv(percentage_change, file = paste0(dir, 'Results/percentage_change-', run_date,'.csv'),
          row.names = FALSE)


#average cross-species response
#extract parameter estimates
mu_gamma_spM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'mu_gamma'))
mu_thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'mu_theta1'))

mu_gamma_spW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'mu_gamma'))
mu_thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'mu_theta1'))


mn_CVmass <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_lat_sc_sim))
mn_CVwing <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_lat_sc_sim))
mn_CVmass10 <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = 2)
mn_CVwing10 <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = 2)

#o_lat_sc_sim from above
#10 degrees of latitude from in middle of lat range for all species
slat10 <- c(mean(o_lat_sc_sim) - 5, mean(o_lat_sc_sim) + 5)

for (i in 1:NROW(mu_gamma_spM))
{
  #entire range of lat values
  mn_CVmass[i,] <- (mu_gamma_spM[i,1] + (mu_thetameanM[i,1] * o_lat_sc_sim))/1000
  mn_CVwing[i,] <- (mu_gamma_spW[i,1] + (mu_thetameanW[i,1] * o_lat_sc_sim))/1000

  #per 10 degrees lat
  mn_CVmass10[i,] <- (mu_gamma_spM[i,1] + (mu_thetameanM[i,1] * slat10))/1000
  mn_CVwing10[i,] <- (mu_gamma_spW[i,1] + (mu_thetameanW[i,1] * slat10))/1000
}

hist(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
hist(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))

hist(apply(mn_CVmass10, 1, function(x) (max(x) - min(x)) / min(x)))
hist(apply(mn_CVwing10, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVmass10, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVwing10, 1, function(x) (max(x) - min(x)) / min(x)))


# Calculate the % change over Distance to range edge for each species --------------------------------------------

#extract parameter estimates
gamma_spM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'gamma'))
thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'theta2'))

gamma_spW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'gamma'))
thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'theta2'))

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
dist_data <- raw_data_mass %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(meanDist = mean(as.numeric(distance_edge_km_buffer10)),
                   sdDist = sd(as.numeric(distance_edge_km_buffer10)),
                   Distminim = min(as.numeric(distance_edge_km_buffer10)),
                   Distmaxim = max(as.numeric(distance_edge_km_buffer10)),
                   RangeDist = Distmaxim-Distminim,
                   minDist_str = min(as.numeric(ldre_sc)),
                   maxDist_str = max(as.numeric(ldre_sc)),
                   MinMaxDistSp = n())
dist_data <- data.frame(dist_data)

result_list <- list()
percentage_changeM <- percentage_changeW <- as.data.frame(matrix(NA, nrow=nrow(dist_data), ncol = 5))
colnames(percentage_changeM) <- colnames(percentage_changeW) <- c('sci_name', 'mean_change', 'q055', 'q945','trait') 
percentage_changeM[, 'sci_name'] <- percentage_changeW[, 'sci_name'] <- dist_data$sci_name
percentage_changeM[, 'trait'] <- rep('mass', nrow(dist_data))
percentage_changeW[, 'trait'] <- rep('wing', nrow(dist_data))


for (i in 1:nrow(dist_data)) {

  print(paste0(i, ' of ', nrow(dist_data)))

  tmpDist <- seq(dist_data$minDist_str[i], dist_data$maxDist_str[i], length.out = 2)

  tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpDist)))
  tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpDist)))

  for (j in 1:nrow(gamma_spM)){
    tmpCVmass[j,] <- (gamma_spM[j,i] + (thetameanM[j,i] * tmpDist))/1000
    tmpCVwing[j,] <- (gamma_spW[j,i] + (thetameanW[j,i] * tmpDist))/1000
  }

  # what I need is only the predictions of the min lat and max lat (first and last columns)
  # Combine tmpDist and tmpCV into a data frame
  tmp_df <- data.frame(tmpCVmass = tmpCVmass[,c(1, length(tmpDist))], 
                       tmpCVwing = tmpCVwing[,c(1, length(tmpDist))])

  result_list[[i]] <- tmp_df #saving all values in case we need to extract other estimates

  #(predicted CV at max lat - predicted CV at min lat) / predicted CV at min lat = % change
  #changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  #changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  # % change in mass and wing cv per 10 degrees latitude
  changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  percentage_changeM[i, "mean_change"] <- mean(changeMass)
  percentage_changeM[i, "q055"] <- quantile(changeMass, probs = 0.055)
  percentage_changeM[i, "q945"] <- quantile(changeMass, probs = 0.945)

  percentage_changeW[i, "mean_change"] <- mean(changeWing)
  percentage_changeW[i, "q055"] <- quantile(changeWing, probs = 0.055)
  percentage_changeW[i, "q945"] <- quantile(changeWing, probs = 0.945)
}

mean(as.numeric(percentage_changeM$mean_change))
mean(as.numeric(percentage_changeW$mean_change))

#average cross-species response
#extract parameter estimates
mu_thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'mu_theta2'))
mu_thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'mu_theta2'))

mn_CVmass <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_ldre_sc_sim))
mn_CVwing <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_ldre_sc_sim))

#o_ldre_sc_sim from above
for (i in 1:NROW(mu_gamma_spM))
{
  print(paste0(i, ' of ', NROW(mu_gamma_spM)))
  #entire range of lat values
  mn_CVmass[i,] <- (mu_gamma_spM[i,1] + (mu_thetameanM[i,1] * o_ldre_sc_sim))/1000
  mn_CVwing[i,] <- (mu_gamma_spW[i,1] + (mu_thetameanW[i,1] * o_ldre_sc_sim))/1000
}

hist(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
hist(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))


# Calculate the % change over Spatial variation in DHI for each species --------------------------------------------

#extract parameter estimates
gamma_spM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'gamma'))
thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'theta3'))

gamma_spW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'gamma'))
thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'theta3'))

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
spat_data <- raw_data_mass %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(meanSpat = mean(as.numeric(DHI_spatial_var_10km)),
                   sdSpat = sd(as.numeric(DHI_spatial_var_10km)),
                   Spatminim = min(as.numeric(DHI_spatial_var_10km)),
                   Spatmaxim = max(as.numeric(DHI_spatial_var_10km)),
                   RangeSpat = Spatmaxim-Spatminim,
                   minSpat_str = min(as.numeric(DHI_spat_sc)),
                   maxSpat_str = max(as.numeric(DHI_spat_sc)),
                   MinMaxSpatSp = n())
spat_data <- data.frame(spat_data)

result_list <- list()
percentage_changeM <- percentage_changeW <- as.data.frame(matrix(NA, nrow=nrow(spat_data), ncol = 5))
colnames(percentage_changeM) <- colnames(percentage_changeW) <- c('sci_name', 'mean_change', 'q055', 'q945','trait') 
percentage_changeM[, 'sci_name'] <- percentage_changeW[, 'sci_name'] <- spat_data$sci_name
percentage_changeM[, 'trait'] <- rep('mass', nrow(spat_data))
percentage_changeW[, 'trait'] <- rep('wing', nrow(spat_data))


for(i in 1:nrow(spat_data)) {

  print(paste0(i, ' of ', nrow(dist_data)))

  tmpSpat <- seq(spat_data$minSpat_str[i], spat_data$maxSpat_str[i], length.out = 2)

  tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpSpat)))
  tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpSpat)))

  for(j in 1:nrow(gamma_spM)){
    tmpCVmass[j,] <- (gamma_spM[j,i] + (thetameanM[j,i] * tmpSpat))/1000
    tmpCVwing[j,] <- (gamma_spW[j,i] + (thetameanW[j,i] * tmpSpat))/1000
  }

  # what I need is only the predictions of the min lat and max lat (first and last columns)
  # Combine tmpDist and tmpCV into a data frame
  tmp_df <- data.frame(tmpCVmass = tmpCVmass[,c(1, length(tmpSpat))], 
                       tmpCVwing = tmpCVwing[,c(1, length(tmpSpat))])

  result_list[[i]] <- tmp_df #saving all values in case we need to extract other estimates

  #(predicted CV at max lat - predicted CV at min lat) / predicted CV at min lat = % change
  #changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  #changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  # % change in mass and wing cv per 10 degrees latitude
  changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  percentage_changeM[i, "mean_change"] <- mean(changeMass)
  percentage_changeM[i, "q055"] <- quantile(changeMass, probs = 0.055)
  percentage_changeM[i, "q945"] <- quantile(changeMass, probs = 0.945)

  percentage_changeW[i, "mean_change"] <- mean(changeWing)
  percentage_changeW[i, "q055"] <- quantile(changeWing, probs = 0.055)
  percentage_changeW[i, "q945"] <- quantile(changeWing, probs = 0.945)

}

mean(as.numeric(percentage_changeM$mean_change))
mean(as.numeric(percentage_changeW$mean_change))


#average cross-species response
#extract parameter estimates
mu_thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'mu_theta3'))
mu_thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'mu_theta3'))

mn_CVmass <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_spat_sc_sim))
mn_CVwing <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_spat_sc_sim))

#o_spat_sc_sim from above
for (i in 1:NROW(mu_gamma_spM))
{
  print(paste0(i, ' of ', NROW(mu_gamma_spM)))
  #entire range of lat values
  mn_CVmass[i,] <- (mu_gamma_spM[i,1] + (mu_thetameanM[i,1] * o_spat_sc_sim))/1000
  mn_CVwing[i,] <- (mu_gamma_spW[i,1] + (mu_thetameanW[i,1] * o_spat_sc_sim))/1000
}

hist(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
hist(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))


# Calculate the % change over Temporal variation in DHI for each species --------------------------------------------

#extract parameter estimates
gamma_spM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'gamma'))
thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'theta4'))

gamma_spW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'gamma'))
thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'theta4'))

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
temp_data <- raw_data_mass %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(meanTemp = mean(as.numeric(DHI_temporal_var_10km)),
                   sdTemp = sd(as.numeric(DHI_temporal_var_10km)),
                   Tempminim = min(as.numeric(DHI_temporal_var_10km)),
                   Tempmaxim = max(as.numeric(DHI_temporal_var_10km)),
                   RangeTemp = Tempmaxim-Tempminim,
                   minTemp_str = min(as.numeric(DHI_temp_sc)),
                   maxTemp_str = max(as.numeric(DHI_temp_sc)),
                   MinMaxTempSp = n())
temp_data <- data.frame(temp_data)

result_list <- list()
percentage_changeM <- percentage_changeW <- as.data.frame(matrix(NA, nrow=nrow(temp_data), ncol = 5))
colnames(percentage_changeM) <- colnames(percentage_changeW) <- c('sci_name', 'mean_change', 'q055', 'q945','trait') 
percentage_changeM[, 'sci_name'] <- percentage_changeW[, 'sci_name'] <- temp_data$sci_name
percentage_changeM[, 'trait'] <- rep('mass', nrow(temp_data))
percentage_changeW[, 'trait'] <- rep('wing', nrow(temp_data))


for(i in 1:nrow(temp_data)) {

  print(paste0(i, ' of ', nrow(dist_data)))

  tmpTemp <- seq(temp_data$minTemp_str[i], temp_data$maxTemp_str[i], length.out = 2)

  tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpTemp)))
  tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpTemp)))

  for(j in 1:nrow(gamma_spM)){
    tmpCVmass[j,] <- (gamma_spM[j,i] + (thetameanM[j,i] * tmpTemp))/1000
    tmpCVwing[j,] <- (gamma_spW[j,i] + (thetameanW[j,i] * tmpTemp))/1000
  }

  # what I need is only the predictions of the min lat and max lat (first and last columns)
  # Combine tmpTemp and tmpCV into a data frame
  tmp_df <- data.frame(tmpCVmass = tmpCVmass[,c(1, length(tmpTemp))], 
                       tmpCVwing = tmpCVwing[,c(1, length(tmpTemp))])

  result_list[[i]] <- tmp_df #saving all values in case we need to extract other estimates

  #(predicted CV at max lat - predicted CV at min lat) / predicted CV at min lat = % change
  #changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  #changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  # % change in mass and wing cv 
  changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100

  percentage_changeM[i, "mean_change"] <- mean(changeMass)
  percentage_changeM[i, "q055"] <- quantile(changeMass, probs = 0.055)
  percentage_changeM[i, "q945"] <- quantile(changeMass, probs = 0.945)

  percentage_changeW[i, "mean_change"] <- mean(changeWing)
  percentage_changeW[i, "q055"] <- quantile(changeWing, probs = 0.055)
  percentage_changeW[i, "q945"] <- quantile(changeWing, probs = 0.945)

}

mean(as.numeric(percentage_changeM$mean_change))
mean(as.numeric(percentage_changeW$mean_change))

#average cross-species response
#extract parameter estimates
mu_thetameanM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = 'mu_theta4'))
mu_thetameanW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = 'mu_theta4'))

mn_CVmass <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_temp_sc_sim))
mn_CVwing <- matrix(NA, nrow = NROW(mu_gamma_spM), ncol = length(o_temp_sc_sim))

#o_temp_sc_sim from above
for (i in 1:NROW(mu_gamma_spM))
{
  print(paste0(i, ' of ', NROW(mu_gamma_spM)))
  #entire range of lat values
  mn_CVmass[i,] <- (mu_gamma_spM[i,1] + (mu_thetameanM[i,1] * o_temp_sc_sim))/1000
  mn_CVwing[i,] <- (mu_gamma_spW[i,1] + (mu_thetameanW[i,1] * o_temp_sc_sim))/1000
}

hist(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
hist(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVmass, 1, function(x) (max(x) - min(x)) / min(x)))
mean(apply(mn_CVwing, 1, function(x) (max(x) - min(x)) / min(x)))

                                                      

# Calculate the % change of CV for each species --------------------------------------------

#extract the min/max of observed CV
#(max CV - min CV) / min CV = % change

#within species
obs_cv_mass <- raw_data_mass %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(mean = mean(cv_within_post_mean),
                   q055 = quantile(cv_within_post_mean, probs = 0.055),
                   q945 = quantile(cv_within_post_mean, probs = 0.945),
                   minCVmass = min(cv_within_post_mean),
                   maxCVmass = max(cv_within_post_mean),
                   percChangeMass = ((maxCVmass - minCVmass) / minCVmass)*100)

obs_cv_wing <- raw_data_wing %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(mean = mean(cv_within_post_mean),
                   q055 = quantile(cv_within_post_mean, probs = 0.055),
                   q945 = quantile(cv_within_post_mean, probs = 0.945), 
                   minCVwing = min(cv_within_post_mean),
                   maxCVwing = max(cv_within_post_mean),
                   percChangeWing = ((maxCVwing - minCVwing) / minCVwing)*100)

mean(obs_cv_mass$percChangeMass)
mean(obs_cv_wing$percChangeWing)


## Amount of variation in body mass and wing length 
# a) within populations within species 
# b) among populations within species 
# c) across species

# a) For within-pop variation, take the mean of cv_within_among_sp_post_mean 
#(this is species-specific CV within pops) across all species.

obs_within_pop_var <- cv_data_covs %>%
  dplyr::group_by(trait) %>%
  dplyr::distinct(sp_id, .keep_all = TRUE) %>%
  dplyr::summarize(meanCV = mean(cv_within_among_sp_post_mean),
                   sdCV = sd(cv_within_among_sp_post_mean))
# trait  meanCV     sdCV
# mass   0.0241   0.00586 
# wing   0.00707  0.000743


# b) Among-pop variation, take the mean of cv_across_post_mean 
#(this is species-specific CV among pops [using the mean and sd across pops of each species]) 
#across all species.

obs_among_pop_var <- cv_data_covs %>%
  dplyr::group_by(trait) %>%
  dplyr::distinct(sp_id, .keep_all = TRUE) %>%
  dplyr::summarize(meanCV = mean(cv_across_post_mean),
                   sdCV = sd(cv_across_post_mean))
# trait  meanCV    sdCV
# mass   0.0109   0.00497
# wing   0.00364  0.00177


# c) Among-species variation, take the mean of cv_among_sp_post_mean 
#(this is station-specific CV among species) across all stations.

obs_among_sp_var <- cv_data_covs %>%
  dplyr::group_by(trait) %>%
  dplyr::distinct(station_id, .keep_all = TRUE) %>%
  dplyr::summarize(meanCV = mean(cv_among_sp_post_mean, na.rm = T),
                   sdCV = sd(cv_among_sp_post_mean, na.rm = T))
# trait  meanCV   sdCV
# mass   0.208   0.0666
# wing   0.0514  0.0192


           
# Density plot for the variation among species --------------------------------------------

mass <- ggplot(raw_data_mass, aes(x = cv_within_post_mean, group = as.factor(sp_id)))  +
  geom_density(color = "#482677FF", linetype = "solid", linewidth = 0.3) +
  labs(x = 'CV within populations',
       y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))

wing <- ggplot(raw_data_wing, aes(x = cv_within_post_mean, group = as.factor(sp_id)))  +
  geom_density(color = "#29AF7FFF", linetype = "solid", linewidth = 0.3) +
  labs(x = 'CV within populations',
       y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))


pdf(paste0('Results/ExtendedDataS3-', run_date, '.pdf'), width = 6, height = 10)
gridExtra::grid.arrange(mass, wing, nrow=2)
dev.off()


