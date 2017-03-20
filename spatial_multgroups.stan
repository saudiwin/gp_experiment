
data {
  
  //N = num obs
  
  int<lower=0> N;
  
  //J= num predictors
  
  int<lower=0> J;
  
  //x = matrix of dummy variables for discrete predictors and continuous variables
  
  matrix[N,J] X;
  
  //longitude outcome
  
  vector[N] longitude;
  
  //latitude outcome
  
  vector[N] latitude;
  
  //number of clusters in data;
  
  int G;
  
  // cluster variable
  
  int cl[N];

}


parameters {
  
  // Estimate a single parameter for each algorithm indexed by J
  // Will permit us to discriminate between algorithms
  // Higher coefficients should indicate better models. 
  // Lower standard errors also equal better-performing models
  // Per Carpenter (2015), the centered parameterization is more efficient given large data,
  // While a non-centered parameterization would be optimal for small data 
  // See http://mc-stan.org/documentation/case-studies/pool-binary-trials.html

  vector[J] theta_lat;
  vector[J] theta_lon;
  
  //Because of standardization, need to drop a cluster for identifiability
  
  vector[G] loc_lat;
  vector[G] loc_lon;
  
  real mu_lat;
  real mu_lon;
  
  real<lower=0> tau_lat;
  real<lower=0> tau_lon;

  //incidental variance parameters
  real<lower=0> sigma_lat;
  real<lower=0> sigma_lon;
}

transformed parameters {
  //Drop one cluster for identifiability
  /*
  vector[G] loc_lat_full;
  vector[G] loc_lon_full;
  
  loc_lat_full = append_row(loc_lat,0);
  loc_lon_full = append_row(loc_lon,0);
  */
}


model {
  
  //Relatively vague priors
  //We assume that the thetas come from a common distribution as they all should be 
  //measuring the same thing
    
  //On the inverse logit scale, this hyperprior on mu covers the full probability range [0,1]
  
    loc_lat ~ normal(mu_lat,tau_lon);
    loc_lon ~ normal(mu_lon,tau_lat);
    
  //tau variance (common-wisdom) prior covers a wide variance on the logit scale
    tau_lon ~ normal(0,5);
    tau_lat ~ normal(0,5);
 
    mu_lat ~ normal(0,5);
    mu_lon ~ normal(0,5); 

    sigma_lon ~ normal(0,5);
    sigma_lat ~ normal(0,5);
  //if non-centered, theta receives a standard unit normal prior to decouple the upper-level from lower-level parameters (non-centered
  
  theta_lat ~ normal(0,5);
  theta_lon ~ normal(0,5);

  //Model sampling statement -- bernoulli model with logit link function (equivalent to GLM with logit link)
  for(n in 1:N) {
  latitude[n] ~ normal(loc_lat[cl[n]] + X[n,]*theta_lat,sigma_lat);
  longitude[n] ~ normal(loc_lon[cl[n]] + X[n,]*theta_lon,sigma_lon);
  }
  // 
  // 
}


generated quantities {
  /*
  vector[N] lon_predict;
  vector[N] lat_predict;
  
  for(n in 1:N) {
    lon_predict[n] = normal_rng(X[n,]*theta_lat,sigma_lat);
    lat_predict[n] = normal_rng(X[n,]*theta_lon,sigma_lon);
  }  
  */
}