#######################
# Read files and plot results from model #5 for mass and wing

# Model #5: Within pops among species variation on CV

# mu_cv = mu_gamma + 
#   alpha +                        #species-specific random effect to account for the phylo
#   theta1 * rangesize +           #range size (in km2)
#   theta2 * genlength +           #generation length
#   theta3 * HWI +                 #hand-wing index
#   theta4 * migstatus             #migratory (1) or non-migratory (0) species

#######################


# load packages -----------------------------------------------------------

library(ggplot2)
library(MCMCvis)
library(tidyverse)
library(cowplot)
library(ape)
library(phytools)
library(scales)
library(here)


# set dirs ----------------------------------------------------------------

run_date <- Sys.Date()
model_run_date <- '2025-01-10'
tree_date <- '2025-01-10'
dir <- paste0(here::here(), '/') #location of working directory


# read in data -----------------

#mass
fitMass <- readRDS(paste0(dir, 'Results/fitMassCV-within-among-sp-', model_run_date, '/fitMass-fit-', model_run_date, '.rds'))
DATA_mass <- readRDS(paste0(dir, 'Results/fitMassCV-within-among-sp-', model_run_date, '/fitMass-data-', model_run_date, '.rds'))
raw_data_mass <- readRDS(paste0(dir, 'Results/fitMassCV-within-among-sp-', model_run_date, '/raw-mass-data-', model_run_date, '.rds'))


#wing
fitWing <- readRDS(paste0(dir, 'Results/fitWingCV-within-among-sp-', model_run_date, '/fitWing-fit-', model_run_date, '.rds'))
DATA_wing <- readRDS(paste0(dir, 'Results/fitWingCV-within-among-sp-', model_run_date, '/fitWing-data-', model_run_date, '.rds'))
raw_data_wing <- readRDS(paste0(dir, 'Results/fitWingCV-within-among-sp-', model_run_date, '/raw-wing-data-', model_run_date, '.rds'))


#phylo tree
tree <- readRDS(paste0(dir, 'Data/L1/consensus_tree-', tree_date, '.rds'))


#check result summary

#mass
MCMCsummary(fitMass, c("mu_gamma", "theta1", "theta2", "theta3", "theta4"), pg0=T)

#                 mean         sd        2.5%        50%       97.5% Rhat n.eff_bulk n.eff_tail  p>0
# mu_gamma 22.1504774 2.5880491 17.1123875 22.1495500 27.2555750    1       5578       6800 1.00
# theta1    0.2611775 0.5468504 -0.8320906  0.2666160  1.3242040    1       9561       8060 0.69
# theta2   -2.5538044 0.5367993 -3.5795190 -2.5569550 -1.4703725    1       7955       7844 0.00
# theta3   -0.2732020 0.5769200 -1.4217150 -0.2675500  0.8495275    1       8157       7938 0.31
# theta4   -0.7644900 2.0702553 -4.7489473 -0.7873445  3.3122495    1       6486       6518 0.35

# pg0=T -> proportion of the posterior that is greater than 0


#wing
MCMCsummary(fitWing, c("mu_gamma", "theta1", "theta2", "theta3", "theta4"), pg0=T)

#                 mean          sd        2.5%          50%        97.5% Rhat n.eff_bulk n.eff_tail  p>0
# mu_gamma  6.42169996 0.29450806  5.82090950  6.4298900  6.980614250    1       4908       5372 1.00
# theta1   -0.01393761 0.08561498 -0.18005518 -0.0136686  0.151937450    1       7868       8124 0.44
# theta2   -0.17837259 0.08961864 -0.35093212 -0.1798025 -0.003178685    1       7884       7265 0.02
# theta3   -0.44069719 0.08765283 -0.61127610 -0.4410320 -0.264878275    1       7759       7034 0.00
# theta4    0.50045267 0.29157786 -0.06987071  0.5024930  1.059262000    1       5131       6152 0.96


# Extracting values of each parameter -----------------

#mass
thetasM <- as.data.frame(MCMCvis::MCMCchains(fitMass, params = c('theta1', 'theta2',
                                                   'theta3', 'theta4'))) %>%
  dplyr::summarise(across(everything(), list(
    mean = ~ round(mean(.x), 3),
    q055 = ~ round(quantile(.x, probs = 0.055), 3),
    q945 = ~ round(quantile(.x, probs = 0.945), 3)))) %>%
  matrix(nrow = 4, byrow = TRUE, dimnames = list(c('theta1', 'theta2', 'theta3', 'theta4'), 
                                                 c('mean', 'q005', 'q945'))) %>%
  as.data.frame() %>%
  dplyr::mutate('p>0' = cbind(MCMCvis::MCMCsummary(fitMass, params = c('theta1', 'theta2', 
                                                               'theta3', 'theta4'), pg0=T)$'p>0'))

#          mean   q005   q945  p>0
# theta1  0.261 -0.615  1.126 0.69  #rangesize
# theta2 -2.554 -3.402 -1.689 0.00  #genlength
# theta3 -0.273 -1.187  0.632 0.31  #HWI
# theta4 -0.764 -4.034  2.531 0.35  #migration


#wing
thetasW <- as.data.frame(MCMCvis::MCMCchains(fitWing, params = c('theta1', 'theta2',
                                                                 'theta3', 'theta4'))) %>%
  dplyr::summarise(across(everything(), list(
    mean = ~ round(mean(.x), 3),
    q055 = ~ round(quantile(.x, probs = 0.055), 3),
    q945 = ~ round(quantile(.x, probs = 0.945), 3)))) %>%
  matrix(nrow = 4, byrow = TRUE, dimnames = list(c('theta1', 'theta2', 'theta3', 'theta4'), 
                                                 c('mean', 'q005', 'q945'))) %>%
  as.data.frame() %>%
  dplyr::mutate('p>0' = cbind(MCMCvis::MCMCsummary(fitWing, params = c('theta1', 'theta2', 
                                                                       'theta3', 'theta4'), pg0=T)$'p>0'))

#         mean   q005   q945  p>0
# theta1 -0.014 -0.152  0.121 0.44  #rangesize
# theta2 -0.178  -0.32 -0.035 0.02  #genlength
# theta3 -0.441  -0.58 -0.301 0.00  #HWI
# theta4    0.5  0.034  0.959 0.96  #migration



# Plot effects of each variable -----------------

# Range Size

#extract parameter estimates
gamma_spM <- MCMCvis::MCMCchains(fitMass, params = 'mu_gamma')
thetameanM <- MCMCvis::MCMCchains(fitMass, params = 'theta1')

gamma_spW <- MCMCvis::MCMCchains(fitWing, params = 'mu_gamma')
thetameanW <- MCMCvis::MCMCchains(fitWing, params = 'theta1')

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
range_data <- raw_data_mass %>%
  dplyr::summarize(meanRange = mean(as.numeric(range_size_log)),
                   sdRange = sd(as.numeric(range_size_log)),
                   Rangeminim = min(as.numeric(range_size_log)),
                   Rangemaxim = max(as.numeric(range_size_log)),
                   minRange_str = min(as.numeric(range_size_str)),
                   maxRange_str = max(as.numeric(range_size_str)))

#create a vector of range values
tmpRange <- seq(range_data$minRange_str, range_data$maxRange_str, length.out = nrow(raw_data_mass))

#loop through each iteration to predict new values of CV mass and wing
tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpRange)))
tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpRange)))

for(i in 1:length(tmpRange)){
  for(j in 1:nrow(gamma_spM)){
    tmpCVmass[j,i] <- (gamma_spM[j,1] + (thetameanM[j,1] * tmpRange[i]))/1000
    tmpCVwing[j,i] <- (gamma_spW[j,1] + (thetameanW[j,1] * tmpRange[i]))/1000
  }
}

#extract means and quantiles from predicted values
pred_values_range <- data.frame(tmpRange = tmpRange, 
                          tmpCVmass = apply(tmpCVmass, 2, mean),
                          tmpCVmass_upper = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.945)),
                          tmpCVmass_lower = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.055)),
                          tmpCVwing = apply(tmpCVwing, 2, mean),
                          tmpCVwing_upper = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.945)),
                          tmpCVwing_lower = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.055)),
                          sci_name = raw_data_mass$sci_name)

#backscale the range variable
pred_values_range$tmpRangeBackscaled <- exp((pred_values_range$tmpRange * range_data$sdRange) + range_data$meanRange)
pred_values_range$tmpRangeBackscaledLog <- (pred_values_range$tmpRange * range_data$sdRange) + range_data$meanRange


PercChange <- pred_values_range %>%
  dplyr::filter(tmpRangeBackscaled == min(tmpRangeBackscaled, na.rm = TRUE) | 
                  tmpRangeBackscaled == max(tmpRangeBackscaled, na.rm = TRUE)) %>%
  dplyr::select(tmpRangeBackscaled, tmpCVmass, tmpCVwing) %>%
  dplyr::mutate(change = c((max(tmpCVmass) - min(tmpCVmass)) / min(tmpCVmass) * 100,
                           ((max(tmpCVwing) - min(tmpCVwing)) / min(tmpCVwing) * 100)))

                                                  
# Plot
RangeSize <- ggplot() + 
  #mass predictions
  geom_line(data = pred_values_range, aes(x = tmpRangeBackscaled, y = tmpCVmass), 
            color = "#482677FF", lwd = 1.2) + 
  geom_ribbon(data = pred_values_range, aes(x = tmpRangeBackscaled, ymin = tmpCVmass_lower, ymax = tmpCVmass_upper), 
              fill = "#482677FF", alpha = 0.2) +
  
  #wing predictions
  geom_line(data = pred_values_range, aes(x = tmpRangeBackscaled, y = tmpCVwing), 
            color = "#29AF7FFF", lwd = 1.2) + 
  geom_ribbon(data = pred_values_range, aes(x = tmpRangeBackscaled, ymin = tmpCVwing_lower, ymax = tmpCVwing_upper), 
              fill = "#29AF7FFF", alpha = 0.2) +
  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab(expression("Range size (km"^2*")")) + 
  ylab(NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  scale_x_log10(labels = scales::label_number(), breaks = c(300000, 3000000, 30000000)) +
  #scale_y_continuous(breaks = c(0.01, 0.02, 0.03)) +
  scale_y_log10(limits = c(0.005, 0.032), breaks = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03)) +
  #coord_cartesian(ylim = c(0, 0.032)) +
  ggtitle('C')


# Generation Length

#extract parameter estimates
gamma_spM <- MCMCvis::MCMCchains(fitMass, params = 'mu_gamma')
thetameanM <- MCMCvis::MCMCchains(fitMass, params = 'theta2')

gamma_spW <- MCMCvis::MCMCchains(fitWing, params = 'mu_gamma')
thetameanW <- MCMCvis::MCMCchains(fitWing, params = 'theta2')

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
genL_data <- raw_data_mass %>%
  dplyr::summarize(meanGenL = mean(as.numeric(GenLength_log)),
                   sdGenL = sd(as.numeric(GenLength_log)),
                   GenLminim = min(as.numeric(GenLength_log)),
                   GenLmaxim = max(as.numeric(GenLength_log)),
                   minGenL_str = min(as.numeric(GenLength_str)),
                   maxGenL_str = max(as.numeric(GenLength_str)))

#create a vector of range values
tmpGenL <- seq(genL_data$minGenL_str, genL_data$maxGenL_str, length.out = nrow(raw_data_mass))

#loop through each iteration to predict new values of CV mass and wing
tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpGenL)))
tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpGenL)))

for(i in 1:length(tmpGenL)){
  for(j in 1:nrow(gamma_spM)){
    tmpCVmass[j,i] <- (gamma_spM[j,1] + (thetameanM[j,1] * tmpGenL[i]))/1000
    tmpCVwing[j,i] <- (gamma_spW[j,1] + (thetameanW[j,1] * tmpGenL[i]))/1000
  }
}

#extract means and quantiles from predicted values
pred_values_gen <- data.frame(tmpGenL = tmpGenL, 
                          tmpCVmass = apply(tmpCVmass, 2, mean),
                          tmpCVmass_upper = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.945)),
                          tmpCVmass_lower = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.055)),
                          tmpCVwing = apply(tmpCVwing, 2, mean),
                          tmpCVwing_upper = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.945)),
                          tmpCVwing_lower = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.055)),
                          sci_name = raw_data_mass$sci_name)

#backscale the generation length variable
pred_values_gen$tmpGenBackscaled <- exp((pred_values_gen$tmpGenL * genL_data$sdGenL) + genL_data$meanGenL)

# % change in mass and wing cv for different generation times
PercChange <- pred_values_gen %>%
  dplyr::filter(tmpGenBackscaled == min(tmpGenBackscaled, na.rm = TRUE) | 
                  tmpGenBackscaled == max(tmpGenBackscaled, na.rm = TRUE)) %>%
  dplyr::select(tmpGenBackscaled, tmpCVmass, tmpCVwing) %>%
  dplyr::mutate(change = c(((max(tmpCVmass) - min(tmpCVmass)) / min(tmpCVmass))*100,
                           ((max(tmpCVwing) - min(tmpCVwing)) / min(tmpCVwing))*100))
                                                  

# Plot
GenLength <- ggplot() +
  #mass predictions
  geom_line(data = pred_values_gen, aes(x = tmpGenBackscaled, y = tmpCVmass), 
            lwd = 1.2, color = "#482677FF") +
  geom_ribbon(data = pred_values_gen, aes(x = tmpGenBackscaled, ymin = tmpCVmass_lower, ymax = tmpCVmass_upper), 
              fill = "#482677FF", alpha = 0.2) +
  
  #wing predictions
  geom_line(data = pred_values_gen, aes(x = tmpGenBackscaled, y = tmpCVwing), 
            lwd = 1.2, color = "#29AF7FFF") +
  geom_ribbon(data = pred_values_gen, aes(x = tmpGenBackscaled, ymin = tmpCVwing_lower, ymax = tmpCVwing_upper), 
              fill = "#29AF7FFF", alpha = 0.2) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Generation Time (years)") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  #scale_y_continuous(breaks = c(0.01, 0.02, 0.03)) +
  #coord_cartesian(ylim = c(0, 0.032)) +
  scale_y_log10(limits = c(0.005, 0.032), breaks = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03)) +
  ggtitle('A')


# Hand wing index

#extract parameter estimates
gamma_spM <- MCMCvis::MCMCchains(fitMass, params = 'mu_gamma')
thetameanM <- MCMCvis::MCMCchains(fitMass, params = 'theta3')

gamma_spW <- MCMCvis::MCMCchains(fitWing, params = 'mu_gamma')
thetameanW <- MCMCvis::MCMCchains(fitWing, params = 'theta3')

#extract the mean/sd of original data (log) + min and max of scaled data - scale data ends with '_str'
HWI_data <- raw_data_mass %>%
  dplyr::summarize(meanHWI = mean(as.numeric(HWI_log)),
                   sdHWI = sd(as.numeric(HWI_log)),
                   HWIminim = min(as.numeric(HWI_log)),
                   HWImaxim = max(as.numeric(HWI_log)),
                   minHWI_str = min(as.numeric(HWI_str)),
                   maxHWI_str = max(as.numeric(HWI_str)))

#create a vector of range values
tmpHWI <- seq(HWI_data$minHWI_str, HWI_data$maxHWI_str, length.out = nrow(raw_data_mass))

#loop through each iteration to predict new values of CV mass and wing
tmpCVmass <- as.data.frame(matrix(NA, nrow = nrow(gamma_spM), ncol = length(tmpHWI)))
tmpCVwing <- as.data.frame(matrix(NA, nrow = nrow(gamma_spW), ncol = length(tmpHWI)))

for(i in 1:length(tmpHWI)){
  for(j in 1:nrow(gamma_spM)){
    tmpCVmass[j,i] <- (gamma_spM[j,1] + (thetameanM[j,1] * tmpHWI[i]))/1000
    tmpCVwing[j,i] <- (gamma_spW[j,1] + (thetameanW[j,1] * tmpHWI[i]))/1000
  }
}

#extract means and quantiles from predicted values
pred_values_hwi <- data.frame(tmpHWI = tmpHWI, 
                              tmpCVmass = apply(tmpCVmass, 2, mean),
                              tmpCVmass_upper = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.945)),
                              tmpCVmass_lower = apply(tmpCVmass, 2, function(x) quantile(x, probs=0.055)),
                              tmpCVwing = apply(tmpCVwing, 2, mean),
                              tmpCVwing_upper = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.945)),
                              tmpCVwing_lower = apply(tmpCVwing, 2, function(x) quantile(x, probs=0.055)),
                              sci_name = raw_data_mass$sci_name)

#backscale the generation length variable
pred_values_hwi$tmpHWIBackscaled <- exp((pred_values_hwi$tmpHWI * HWI_data$sdHWI) + HWI_data$meanHWI)

# % change in mass and wing cv for different generation times
PercChange <- pred_values_hwi %>%
  dplyr::filter(tmpHWIBackscaled == min(tmpHWIBackscaled, na.rm = TRUE) | 
                  tmpHWIBackscaled == max(tmpHWIBackscaled, na.rm = TRUE)) %>%
  dplyr::select(tmpHWIBackscaled, tmpCVmass, tmpCVwing) %>%
  dplyr::mutate(change = c(((max(tmpCVmass) - min(tmpCVmass)) / min(tmpCVmass))*100,
                           ((max(tmpCVwing) - min(tmpCVwing)) / min(tmpCVwing))*100))
                                                   
# Plot
HWI <- ggplot() +
  geom_line(data = pred_values_hwi, aes(x = tmpHWIBackscaled, y = tmpCVmass), 
            lwd = 1.2, color = "#482677FF") +
  geom_ribbon(data = pred_values_hwi, aes(x = tmpHWIBackscaled, ymin = tmpCVmass_lower, ymax = tmpCVmass_upper), 
              fill = "#482677FF", alpha = 0.2) +
  
  geom_line(data = pred_values_hwi, aes(x = tmpHWIBackscaled, y = tmpCVwing), 
            lwd = 1.2, color = "#29AF7FFF") +
  geom_ribbon(data = pred_values_hwi, aes(x = tmpHWIBackscaled, ymin = tmpCVwing_lower, ymax = tmpCVwing_upper), 
              fill = "#29AF7FFF", alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Hand-wing Index") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  #scale_y_continuous(breaks = c(0.01, 0.02, 0.03)) +
  #coord_cartesian(ylim = c(0, 0.032)) +
  scale_y_log10(limits = c(0.004, 0.032), breaks = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03)) +
  ggtitle('B')



# Migratory status = 1

#extract parameter estimates
mig_status_data <- data.frame(thetameanM = as.vector(MCMCvis::MCMCchains(fitMass, params = 'theta4')),
                              thetameanW = as.vector(MCMCvis::MCMCchains(fitWing, params = 'theta4')))

mig_status_results <- data.frame(variable = c("Mass", "Wing"),
  mean = c(mean(mig_status_data$thetameanM), 
           mean(mig_status_data$thetameanW)),
  q055 = c(quantile(mig_status_data$thetameanM, probs = 0.055), 
           quantile(mig_status_data$thetameanW, probs = 0.055)),
  q945 = c(quantile(mig_status_data$thetameanM, probs = 0.945), 
           quantile(mig_status_data$thetameanW, probs = 0.945)))

# Reorder the variable column to have "Wing" at the bottom
mig_status_results$variable <- factor(mig_status_results$variable, levels = c("Wing", "Mass"))

# For mig_status = 0 
pred_mass_res <- mean(gamma_spM/1000)
pred_wing_res <- mean(gamma_spW/1000)

# For mig_status = 1 (e.g. migrant)
pred_mass_mig <- (gamma_spM + mig_status_data$thetameanM)/1000
pred_wing_mig <- (gamma_spW + mig_status_data$thetameanW)/1000

mean(((pred_mass_res - pred_mass_mig) / pred_mass_mig )*100)
# 4.665147
                                                      
#plot
MigratStatus <- ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "gray80", linewidth = 1) +
  
  geom_pointrange(data = mig_status_results, aes(x = variable,
                                                 y = mean, 
                                                 ymin = q055,
                                                 ymax = q945,
                                                 color = variable), 
                  linewidth = 1.2, fatten = 5) +
  scale_color_manual(values = c('Mass' = '#482677FF', 'Wing' = '#29AF7FFF'),
                     labels = c("Wing length", "Body mass"),
                     name = NULL) +
  theme_minimal() +
  labs(x = NULL, y = 'Effect of migration') +
  theme(axis.text.x = element_text(size = 12),
        #axis.text.y = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(t = 10)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.2, 0.15),
        legend.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 18)) +
  coord_flip() +
  ggtitle('D')


#join all four plots and export as a pdf

pdf(paste0('Results/Fig4-results-within-pops-among-sp-', run_date, '.pdf'), width = 8, height = 8)

combined_plot <- cowplot::plot_grid(GenLength, HWI, RangeSize, MigratStatus, align = "hv", ncol = 2)
cowplot::ggdraw(combined_plot) +
  draw_label("Within-population phenotypic variation", 
             y = 0.5, x = -0.03, vjust = 1, size = 14, angle = 90) +
  theme(plot.margin = margin(20, 20, 20, 25))  # Adjust margins if needed

dev.off()




# Supplemental information --------------------------------------------


#Figure S2. PPC (multipanel - 2 panels)

# Plot the Posterior Predictive Check (PPC) --------------------------------------------

# PPC mass

y_val <- DATA_mass$CVobs
y_rep <- MCMCvis::MCMCchains(fitMass, params = 'CVobs_rep')

pdf(paste0(dir, 'Results/', 'PPC-within-pops-among-sp-', run_date, '.pdf'), height = 5, width = 8)

par(mfrow=c(1,2))

plot(density(y_val), lwd = 2, main = 'PPC - Mass CV within pops among sp', xlab = 'Value', ylim = c(0, 0.1))
for (i in 1:100){
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.1))
}
# Add legend
legend("topleft", 
       legend = c("CV observed", "CV predicted"), 
       col = c("black", rgb(1, 0, 0, 0.1)), 
       lwd = c(2, 2), 
       bty = "n")  # 'bty = "n"' no box around the legend


# PPC wing 
y_val <- DATA_wing$CVobs
y_rep <- MCMCvis::MCMCchains(fitWing, params = 'CVobs_rep')

plot(density(y_val), lwd = 2, main = 'PPC - Wing CV within pops among sp', xlab = 'Value', ylim = c(0, 0.6))
for (i in 1:100){
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.1))
}
# Add legend
legend("topleft", 
       legend = c("CV observed", "CV predicted"), 
       col = c("black", rgb(1, 0, 0, 0.1)), 
       lwd = c(2, 2), 
       bty = "n")  # 'bty = "n"' no box around the legend

dev.off()



# Plot the Prior Posterior Overlap (PPO) --------------------------------------------

#PPO mass

#Simulate values for the priors.
#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mu_gamma = rnorm(10000, 20, 20),
                    theta1 = rnorm(10000, 0, 5),
                    theta2 = rnorm(10000, 0, 5),
                    theta3 = rnorm(10000, 0, 5),
                    theta4 = rnorm(10000, 0, 10),
                    sigma = rnorm(10000, 0, 2),
                    sigma_phylo = rnorm(10000, 1, 2))

MCMCtrace(fitMass, 
          params = c("mu_gamma", "theta1", "theta2", "theta3", 
                     "theta4", "sigma", 'sigma_phylo'),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE,
          pdf = TRUE, 
          filename = paste0('PPO-mass-within-pops-among-sp-', run_date, '.pdf'),
          wd = paste0(dir, 'Results'))


# PPO wing 

#Parameters are in alphabetical order and identical to the list of parameters at the MCMCtrace function.
SimsPriors <- cbind(mu_gamma = rnorm(10000, 5, 3),
                    theta1 = rnorm(10000, 0, 1),
                    theta2 = rnorm(10000, 0, 1),
                    theta3 = rnorm(10000, 0, 1),
                    theta4 = rnorm(10000, 0, 2),
                    sigma = rnorm(10000, 0, 1),
                    sigma_phylo = rnorm(10000, 1, 2))

MCMCtrace(fitWing, 
          params = c("mu_gamma", "theta1", "theta2", "theta3", "theta4",
                     "sigma", 'sigma_phylo'),
          ISB = FALSE,
          exact = TRUE,
          priors = SimsPriors,
          Rhat = TRUE,
          n.eff = FALSE,
          post_zm = FALSE,
          pdf = TRUE, 
          filename = paste0('PPO-wing-within-pops-among-sp-', run_date, '.pdf'),
          wd = paste0(dir, 'Results'))



# Plot phylo tree with average within-pops variation --------------------------------------------

#mass
cv_data_mass <- raw_data_mass %>% 
  dplyr::select(sp_id, sci_name, cv_within_among_sp_post_mean)

#reorder cv data to match tree order
cv_data_mass$sci_name <- factor(cv_data_mass$sci_name, levels = tree$tip.label)
cv_data_mass <- cv_data_mass[order(cv_data_mass$sci_name), ]
cv_mass_scaled <- scales::rescale(cv_data_mass$cv_within_among_sp_post_mean, to = c(0,1)) + 0.1
names(cv_mass_scaled) <- tree$tip.label

pdf(paste0('Results/ExtendedDataS3Mass-', run_date, '.pdf'), width = 8, height = 8)
mass <- phytools::plotTree.wBars(tree,
                                 x= cv_mass_scaled,
                                 col = '#482677FF',
                                 type= 'fan',
                                 tip.labels = T, fsize = .6,
                                 border = F,
                                 width = 5)
dev.off()


#wing
cv_data_wing <- raw_data_wing %>% 
  dplyr::select(sp_id, sci_name, cv_within_among_sp_post_mean)

#reorder cv data to match tree order
cv_data_wing$sci_name <- factor(cv_data_wing$sci_name, levels = tree$tip.label)
cv_data_wing <- cv_data_wing[order(cv_data_wing$sci_name), ]
cv_wing_scaled <- scales::rescale(cv_data_wing$cv_within_among_sp_post_mean, to = c(0,1)) + 0.1
names(cv_wing_scaled) <- tree$tip.label

pdf(paste0('Results/ExtendedDataS3Wing-', run_date, '.pdf'), width = 8, height = 8)
wing <- phytools::plotTree.wBars(tree,
                                 x= cv_wing_scaled,
                                 col = '#29AF7FFF',
                                 type= 'fan',
                                 tip.labels = T, fsize = .6,
                                 border = F,
                                 width = 5)
dev.off()



                    


# Extra --------------------------------------------

# caterpillar plots with both mass and wing effect sizes

pdf(paste0(dir, 'Results/CV_sp.pdf'), width = 8, height = 5)

param_labels <- c('Range size', 'Generation Length', 
                  'Hand-wing index', 'Migratory status = 1')

MCMCplot(fitMass, fitWing,
         params = c('theta1', 'theta2', 'theta3', 'theta4'), 
         ci = c(50, 89),
         main = 'Variable effect on within pops among species CV', 
         labels = param_labels, col = '#482677FF', col2 = '#29AF7FFF')

legend('bottomright',  # Position of the legend
       legend = c('Mass CV', 'Wing CV'), 
       col = c('#482677FF', '#29AF7FFF'), 
       lty = 1,  # Line type (1 is solid)
       lwd = 3,  bty = 'n')

dev.off()

