
data {
  
  //N = num obs
  
  int<lower=0> N;
  
  //J= num predictors
  
  int<lower=0> J;
  
  //x = matrix of dummy variables for discrete predictors and continuous variables
  
  int bssid[N];
  matrix[J,2] bssid_dim;
  matrix[J,2] bssid_err;
  
  //longitude outcome
  
  vector[N] longitude;
  
  //latitude outcome
  
  vector[N] latitude;
  
  //GPS accuracy
  
  vector[N] acc;
  
  //number of clusters in data;
  
  int G;
  
  // cluster variable
  
  int cl[N];
  
  //distances from beacons
  
  vector[N] distance_ft;

}


parameters {
  
  // Estimate a single parameter for each algorithm indexed by J
  // Will permit us to discriminate between algorithms
  // Higher coefficients should indicate better models. 
  // Lower standard errors also equal better-performing models
  // Per Carpenter (2015), the centered parameterization is more efficient given large data,
  // While a non-centered parameterization would be optimal for small data 
  // See http://mc-stan.org/documentation/case-studies/pool-binary-trials.html
  
  //parameters to estimate, not observed
  matrix[J,2] bssid_dim_par;
  
  //need scaling parameters for bssid
  real scale1_lat;
  real scale1_lon;
  real scale2_lat;
  real scale2_lon;
  //hyper-parameters for bssid_dim
  real mu_err_lon;
  real mu_err_lat;
  real<lower=0>sigma_err_lon;
  real<lower=0>sigma_err_lat;
  
  real dist_par_lon;
  real dist_par_lat;
  
  //GPS accuracy parameters
  real acc_par_lon;
  real acc_par_lat;
  
  //Because of standardization, need to drop a cluster for identifiability
  
  vector[G] loc_lat;
  vector[G] loc_lon;
  
  real mu_lat;
  real mu_lon;
  
  real<lower=0> tau_lat;
  real<lower=0> tau_lon;
  
  
  //interaction parameters
  
  real lat_int;
  real lon_int;
  
  matrix[2,2] med_int;

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
  
    loc_lat ~ normal(mu_lat,tau_lon);
    loc_lon ~ normal(mu_lon,tau_lat);
    
    lat_int ~ normal(0,5);
    lon_int ~ normal(0,5);
    
    tau_lon ~ normal(0,5);
    tau_lat ~ normal(0,5);
    
    mu_err_lon ~ normal(0,5);
    mu_err_lat ~ normal(0,5);
    mu_lat ~ normal(0,5);
    mu_lon ~ normal(0,5); 

    sigma_err_lon ~ normal(0,5);
    sigma_err_lat ~ normal(0,5);
  //if non-centered, theta receives a standard unit normal prior to decouple the upper-level from lower-level parameters (non-centered

  
  bssid_dim_par[,1] ~ normal(mu_err_lon,sigma_err_lon);
  bssid_dim_par[,2] ~ normal(mu_err_lat,sigma_err_lat);
  scale1_lat ~ normal(0,5);
  scale2_lat ~ normal(0,5);
  scale1_lon ~ normal(0,5);
  scale2_lon ~ normal(0,5);
  to_vector(med_int) ~ normal(0,5);
  
  acc_par_lon ~ normal(0,5);
  acc_par_lat ~ normal(0,5);
  
  dist_par_lat ~ normal(0,5);
  dist_par_lon ~ normal(0,5);

  //Model sampling statement -- bernoulli model with logit link function (equivalent to GLM with logit link)
  for(n in 1:N) {
  
  /* Use accuracy to measure scale variance */
  
  
  sigma_lon ~ normal(acc_par_lon * acc[n],5);
  sigma_lat ~ normal(acc_par_lat * acc[n],5);
     
  /* Set up measurement model */
  
  bssid_dim[bssid[n],1] ~ normal(bssid_dim_par[bssid[n],1],bssid_err[bssid[n],1]);
  bssid_dim[bssid[n],2] ~ normal(bssid_dim_par[bssid[n],2],bssid_err[bssid[n],2]);
  
  /* Model longitude and latitude */
  
  latitude[n] ~ normal(loc_lat[cl[n]] + 
                        scale1_lat*bssid_dim_par[bssid[n],1] + 
                        dist_par_lat*distance_ft[n] + 
                        med_int[1,1]*distance_ft[n]*bssid_dim_par[bssid[n],1] +
                        med_int[2,1]*distance_ft[n]*bssid_dim_par[bssid[n],2] +
                        scale2_lat*bssid_dim_par[bssid[n],2]+
                        lat_int*distance_ft[n]*bssid_dim_par[bssid[n],1]*bssid_dim_par[bssid[n],2],
                        sigma_lat);
  
  longitude[n] ~ normal(loc_lon[cl[n]] + 
                          scale1_lon*bssid_dim_par[bssid[n],1] + 
                          dist_par_lon*distance_ft[n] + 
                          med_int[1,2]*distance_ft[n]*bssid_dim_par[bssid[n],1] +
                          med_int[2,2]*distance_ft[n]*bssid_dim_par[bssid[n],2] +
                          scale2_lon*bssid_dim_par[bssid[n],2]+
                          lon_int*distance_ft[n]*bssid_dim_par[bssid[n],1]*bssid_dim_par[bssid[n],2],
                          sigma_lon);
  
  }
  
}


generated quantities {
  
}