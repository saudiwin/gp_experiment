
data {
  
  //N = num obs
  
  int<lower=0> N;
  
  //J= num predictors
  
  int<lower=0> J;
  int<lower=0> W;
  
  int wifi_id[N];
  
  int GPS;
  
  int gps_id[N];
  
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
  
  matrix[J,2] theta_int;
  
  vector[G] loc_lat_raw;
  vector[G] loc_lon_raw;
  
  vector[GPS] gps_fix_raw_lon;
  vector[GPS] gps_fix_raw_lat;
  
  matrix[W,2] wifi_est_raw; 
  
  real mu_gps_lon;
  real mu_gps_lat;
  
  real mu_lat;
  real mu_lon;
  
  real mu_wifi_lat;
  real mu_wifi_lon;
  
  real<lower=0> tau_lat;
  real<lower=0> tau_lon;
  real<lower=0> tau_gps_lon;
  real<lower=0> tau_gps_lat;
  real<lower=0> tau_wifi_lat;
  real<lower=0> tau_wifi_lon;

  //incidental variance parameters
  real<lower=0> sigma_lat;
  real<lower=0> sigma_lon;
}

transformed parameters {
  vector[G] loc_lat;
  vector[G] loc_lon;
  vector[GPS] gps_fix_lon;
  vector[GPS] gps_fix_lat;
  matrix[W,2] wifi_est;
  //de-centered hyperpriors to improve sampling
  loc_lat = mu_lat + loc_lat_raw*tau_lat; 
  loc_lon = mu_lon + loc_lon_raw*tau_lon; 
  gps_fix_lon = mu_gps_lon + gps_fix_raw_lon*tau_gps_lon;
  gps_fix_lat = mu_gps_lat + gps_fix_raw_lat*tau_gps_lat;
  wifi_est[,1] = mu_wifi_lat + wifi_est_raw[,1]*tau_wifi_lat;
  wifi_est[,2] = mu_wifi_lon + wifi_est_raw[,2]*tau_wifi_lon;
}


model {
  
  //Relatively vague priors

  loc_lat_raw ~ normal(0,1);
  loc_lon_raw ~ normal(0,1);
  gps_fix_raw_lon ~ normal(0,1);
  gps_fix_raw_lat ~ normal(0,1);
  to_vector(wifi_est_raw) ~ normal(0,1);
    
    tau_lon ~ normal(0,5);
    tau_lat ~ normal(0,5);
    tau_gps_lon ~ normal(0,5);
    tau_gps_lat ~ normal(0,5);
    tau_wifi_lon ~ normal(0,5);
    tau_wifi_lat ~ normal(0,5);
 
    mu_lat ~ normal(0,10);
    mu_lon ~ normal(0,10); 
    mu_gps_lon ~ normal(0,5);
    mu_gps_lat ~ normal(0,5);

    sigma_lon ~ normal(0,10);
    sigma_lat ~ normal(0,10);
  
  theta_lat ~ normal(0,10);
  theta_lon ~ normal(0,10);
  to_vector(theta_int) ~ normal(0,10);

  //Model sampling statement -
  
  for(n in 1:N) {
    latitude[n] ~ normal(loc_lat[cl[n]] + 
                        gps_fix_lat[gps_id[n]] + 
                        wifi_est[wifi_id[n],1] + 
                        X[n,]*theta_lat + 
                        (X[n,]*theta_int[,1])*wifi_est[wifi_id[n],1],
                        sigma_lat);
                        
    longitude[n] ~ normal(loc_lon[cl[n]] + 
                          gps_fix_lon[gps_id[n]] + 
                          wifi_est[wifi_id[n],2] + 
                          X[n,]*theta_lon +(X[n,]*theta_int[,2])*wifi_est[wifi_id[n],2],
                          sigma_lon);
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
