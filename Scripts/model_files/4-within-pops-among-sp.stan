// Within-pops variation on CV among species
// Code to estimate the effect of dfferent covariates on the within-pop variation across species (across space) 

// Index: 
// k = species
// N_sp: total number of species
// sp_id: species identification for each observation of N_sp

// We extracted the species estimated mean cv and uncertainty from the previous analysis.
// Now we want to model cv as a function of:
// range size
// generation length
// Hand wing index (HWI)
// migratory status: 0 = species that do not migrate (n=10) and 1 = species that migrate (n=88)
// by creating an observation model for the cv as:

// CVobs_k ~ normal(CVact_k, sigmaCV_k), where CVobs_jk and sigmaCV_k are data and CVact_k is non-centered as:
// CVact_k = exp(CVact_raw_k * sigma + mu_cv_k); 
// mu_cv_k = mu_gamma + theta1 * rangesize_k + theta2 * genLeng_k + theta3 * HWI_k + 
           // theta4 * migStatus_k;  //single gamma and thetas estimates



data {
  int<lower=1> N_sp;                      // number of species
  vector[N_sp] CVobs;                     // estimated CV
  vector[N_sp] sigmaCV;                   // estimated CV error
  vector[N_sp] rangesize;                              
  vector[N_sp] genlength;         
  vector[N_sp] HWI;
  vector[N_sp] migstatus;
  matrix[N_sp, N_sp] Rho;  // known correlation matrix
}

transformed data {
  matrix[N_sp, N_sp] LRho = cholesky_decompose(Rho); // get cholesky factor of known corr matrix
}

parameters {
  real mu_gamma;
  real theta1;
  real theta2;
  real theta3;
  real theta4;
  real<lower=0> sigma_phylo;  // phylo sd
  vector[N_sp] z;  // standardized group-level effects
  vector[N_sp] CVact_raw;
  // real<lower=0> tau;  // process error
  real<lower=0> sigma;
}

transformed parameters {
  vector[N_sp] mu_cv;
  vector[N_sp] alpha;  // phylo intercepts (per species)
  vector[N_sp] CVact;
  
  // implies alpha ~ MVN(0, Rho) * sigma_phylo^2
  alpha = (sigma_phylo^2 * (LRho * z));

  mu_cv = mu_gamma + 
  alpha +
  theta1 * rangesize +
  theta2 * genlength + 
  theta3 * HWI +
  theta4 * migstatus;
  
  CVact = CVact_raw * sigma + mu_cv;
}

model {
  // Priors
  mu_gamma ~ normal(5, 3);
  theta1 ~ normal(0, 1);
  theta2 ~ normal(0, 1);
  theta3 ~ normal(0, 1);
  theta4 ~ normal(0, 2);
  sigma_phylo ~ normal(1, 2);
  z ~ std_normal();     //This should be std_normal()
  //tau ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
  CVact_raw ~ std_normal();
 
  // Likelihood
  // observation model (mu_cv is linear predictor, tau is the process error, sigmaCV is the observed uncertainty)
  // CVact ~ normal(mu_cv, tau);
  CVobs ~ normal(CVact, sigmaCV);
}

generated quantities {
  array[N_sp] real CVobs_rep; 
  CVobs_rep = normal_rng(CVact, sigmaCV);
}  
