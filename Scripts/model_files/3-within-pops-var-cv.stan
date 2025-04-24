// Code to estimate the effect of dfferent covariates on the within-pop variation within species (across space) 

data {
  int<lower=1> N_species;                         // number of species
  int<lower=1> N_obs;                         // number of obs
  vector[N_obs] CVobs;                          // estimated CV
  vector[N_obs] sigmaCV;                        // estimated CV error
  array[N_obs] int<lower=1> sp_id; //Species ID for each species/station (cn_id)
  vector[N_obs] lat;                              
  vector[N_obs] distrange;         
  vector[N_obs] spatVar;
  vector[N_obs] tempVar;
}

parameters {
  real<lower=0> sigma;
  real mu_gamma;
  real mu_theta1;
  real mu_theta2;
  real mu_theta3;
  real mu_theta4;

  vector<lower=0>[5] sigma_gt;                  //Change when adding other covs
  cholesky_factor_corr[5] L_Rho_gt;             //Change when adding other covs
  matrix[5, N_species] z_gt;                    //Change when adding other covs             
  vector[N_obs] CVact_raw;
}

transformed parameters {
  vector[N_species] gamma;
  vector[N_species] theta1;
  vector[N_species] theta2;
  vector[N_species] theta3;
  vector[N_species] theta4;
  
  matrix[N_species, 5] gt;             //Change when adding other covs    
  matrix[5, 5] Rho_gt;                 //Change when adding other covs
  vector[N_obs] mu_cv;
  vector[N_obs] CVact;
  
  // cholesky factor of covariance matrix multiplied by z score
  gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
  Rho_gt = multiply_lower_tri_self_transpose(L_Rho_gt);

  gamma = mu_gamma + gt[,1];   
  theta1 = mu_theta1 + gt[,2];
  theta2 = mu_theta2 + gt[,3];
  theta3 = mu_theta3 + gt[,4];
  theta4 = mu_theta4 + gt[,5];

  mu_cv = gamma[sp_id] + 
  theta1[sp_id] .* lat + 
  theta2[sp_id] .* distrange + 
  theta3[sp_id] .* spatVar + 
  theta4[sp_id] .* tempVar; 

  CVact = CVact_raw * sigma + mu_cv;
}

model {
  // Priors
  // for wing
  mu_gamma ~ normal(10, 5);
  mu_theta1 ~ std_normal();
  mu_theta2 ~ std_normal();
  mu_theta3 ~ std_normal();
  mu_theta4 ~ std_normal();
  sigma ~ std_normal();
  
  // for mass
  //mu_gamma ~ normal(25, 5);
  //mu_theta1 ~ std_normal();
  //mu_theta2 ~ std_normal();
  //mu_theta3 ~ normal(0, 5);
  //mu_theta4 ~ normal(0, 5);
  //sigma ~ normal(2, 1);
  
  sigma_gt ~ normal(0, 10);
  to_vector(z_gt) ~ std_normal();
  L_Rho_gt ~ lkj_corr_cholesky(1);         //Change when adding/deleting covs
  
  CVact_raw ~ std_normal();
  
  // Likelihood
  // observation model (mean_cv is linear predictor and sd_cv is process error, which vary by species)
  // CVact ~ normal(mu_cv[cn_id], sigma);
  CVobs ~ normal(CVact, sigmaCV);
}

generated quantities {
  array[N_obs] real CVobs_rep; 
  CVobs_rep = normal_rng(CVact, sigmaCV);
}
