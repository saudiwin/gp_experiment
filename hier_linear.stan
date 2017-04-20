
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
  
  vector[G] loc_lat_raw;
  vector[G] loc_lon_raw;
  
  real mu_lat;
  real mu_lon;
  
  real<lower=0> tau_lat;
  real<lower=0> tau_lon;

  //incidental variance parameters
  real<lower=0> sigma_lat;
  real<lower=0> sigma_lon;
}

transformed parameters {
  vector[G] loc_lat;
  vector[G] loc_lon;
  //de-centered hyperpriors to improve sampling
  loc_lat = mu_lat + loc_lat_raw*tau_lat; 
  loc_lon = mu_lon + loc_lon_raw*tau_lon; 
}


model {
  
  //Relatively vague priors

  loc_lat_raw ~ normal(0,1);
  loc_lon_raw ~ normal(0,1);
    
    tau_lon ~ normal(0,5);
    tau_lat ~ normal(0,5);
 
    mu_lat ~ normal(0,10);
    mu_lon ~ normal(0,10); 

    sigma_lon ~ normal(0,10);
    sigma_lat ~ normal(0,10);
  
  theta_lat ~ normal(0,10);
  theta_lon ~ normal(0,10);

  //Model sampling statement -
  
  for(n in 1:N) {
  latitude[n] ~ normal(loc_lat[cl[n]] + X[n,]*theta_lat,sigma_lat);
  longitude[n] ~ normal(loc_lon[cl[n]] + X[n,]*theta_lon,sigma_lon);
  }

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
