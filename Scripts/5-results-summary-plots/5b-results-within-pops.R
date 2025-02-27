#######################
# Read files and plot results from model #4 for mass and wing

# Model #4: Within pops variation on CV

# mu_cv = gamma[sp_id] + 
#   theta1[sp_id] * lat +           #latitude
#   theta2[sp_id] * distrange +     #distance to range edge (buffer 10km)
#   theta3[sp_id] * spatVar +       #spatial variation on DHI (productivity, buffer of 10km)
#   theta4[sp_id] * tempVar         #temporal variation on DHI (productivity, buffer of 10km)

#######################


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
model_run_date <- '2025-01-21'
cv_data_date <- '2025-01-08'
dir <- '~/Documents/MorphCVLatitude/' #change as needed



# load packages -----------------------------------------------------------

library(ggplot2)
library(MCMCvis)
library(tidyverse)
library(tidyr)



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


pdf(paste0('Results/Fig3-results-within-pops-', run_date, '.pdf'), width = 12, height = 6)
gridExtra::grid.arrange(mass, wing, nrow=1)
dev.off()




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
  changeMass <- ((tmp_df[,2]-tmp_df[,1])/tmp_df[,1])*100
  changeWing <- ((tmp_df[,4]-tmp_df[,3])/tmp_df[,3])*100
  
  percentage_changeM[i, "mean_change"] <- mean(changeMass)
  percentage_changeM[i, "q055"] <- quantile(changeMass, probs = 0.055)
  percentage_changeM[i, "q945"] <- quantile(changeMass, probs = 0.945)
  
  percentage_changeW[i, "mean_change"] <- mean(changeWing)
  percentage_changeW[i, "q055"] <- quantile(changeWing, probs = 0.055)
  percentage_changeW[i, "q945"] <- quantile(changeWing, probs = 0.945)
  
}


# Look at the results:
mean(as.numeric(percentage_changeM$mean_change))
# -8.009595
# ~8% lower CV mass in higher latitudes on average

mean(as.numeric(percentage_changeW$mean_change))
# -2.676539
# 2.7% lower wing CV in higher latitudes on average


#join data frames
percentage_change <- rbind(percentage_changeM, percentage_changeW)


pdf(paste0(dir, 'Results/perc-change-', run_date, '.pdf'), width = 8, height = 18)

#plot results
ggplot() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 0.5) +
  geom_pointrange(data = percentage_change, aes(y = reorder(sci_name, desc(sci_name)),
                                               x = mean_change, xmin = q055, xmax = q945, group = trait, colour = trait), 
                  linewidth = 1.2, fatten = 2, position = position_dodge(width = 0.6)) +
  labs(x = '% Change on trait CV',
    y = NULL,
    title = "Percentage Change in Mass and Wing Length",
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
# 107.3725 % change on body mass across species ranges

mean(obs_cv_wing$percChangeWing)
# 37.52299 % change on wing length across species ranges




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


