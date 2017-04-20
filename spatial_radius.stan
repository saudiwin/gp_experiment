
data {
  
  //N = num obs
  
  int<lower=0> N;
  
  //J= num predictors
  
  int<lower=0> J;
  
  //x = matrix of dummy variables for discrete predictors and continuous variables
  
  int bssid[N];
  
  int num_dist;
  
  int bssid_fix_num;
  
  int bssid_id[num_dist,2];
  
  vector[num_dist] bssid_dist;
  
  matrix[bssid_fix_num,2] bssid_fix_vals;
  
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

transformed data {
  vector[N] distance_ft_sq;
  vector[J-bssid_fix_num] zero;
  
  distance_ft_sq = rows_dot_self(distance_ft);
  for(j in 1:(J-bssid_fix_num)) {
    zero[j] = 0;  
  }
   
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
  matrix[J-bssid_fix_num,2] bssid_dim_par_free;
  matrix[bssid_fix_num,2] bssid_dim_par_fixed;
  /*
  cholesky_factor_corr[J-bssid_fix_num] lcorr1;  
  cholesky_factor_corr[J-bssid_fix_num] lcorr2;
  vector<lower=0>[J-bssid_fix_num] sigma_free1; 
  vector<lower=0>[J-bssid_fix_num] sigma_free2; */
  //need scaling parameters for bssid
  real scale1_lat;
  real scale1_lon;
  
  //hyper-parameters for bssid_dim
  real mu_err_lon;
  real mu_err_lat;
  real<lower=0>sigma_err_lon;
  real<lower=0>sigma_err_lat;
  
  //GPS accuracy parameters
  real acc_par_lon;
  real acc_par_lat;
  
  vector[N] latitude_meas;
  vector[N] longitude_meas;
  
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
  real<lower=0> sigma_dist;
  real<lower=0> sigma_mds;
}

transformed parameters {
  //create combined distance coefficients/locations matrix
  matrix[J,2] bssid_full;
  
  bssid_full = append_row(bssid_dim_par_fixed,bssid_dim_par_free);
  
}


model {
  
  //Relatively vague priors
  /*
    sigma_free1 ~ cauchy(0, 5);
    sigma_free2 ~ cauchy(0, 5);
    
    lcorr1 ~ lkj_corr_cholesky(4);
    lcorr2 ~ lkj_corr_cholesky(4); */
    
    loc_lat ~ normal(mu_lat,tau_lon);
    loc_lon ~ normal(mu_lon,tau_lat);
    
    mu_err_lon ~ normal(0,5);
    mu_err_lat ~ normal(0,5);
    mu_lat ~ normal(0,5);
    mu_lon ~ normal(0,5);
    acc_par_lon ~ normal(0,5);
    acc_par_lat ~ normal(0,5);

    sigma_dist ~ normal(0,10);
    sigma_mds ~ normal(0,1  );
  //if non-centered, theta receives a standard unit normal prior to decouple the upper-level from lower-level parameters (non-centered

  
  /*bssid_dim_par_free[,1] ~ multi_normal_cholesky(zero,
    diag_pre_multiply(sigma_free1,lcorr1));
  bssid_dim_par_free[,2] ~ multi_normal_cholesky(zero,
    diag_pre_multiply(sigma_free2,lcorr2)); */
    to_vector(bssid_dim_par_free) ~ normal(0,100);
  //to_vector(bssid_dim_par_fixed) ~ normal(to_vector(bssid_fix_vals),.01);
  to_vector(bssid_dim_par_fixed) ~ normal(0,1);
  for(i in 1:num_dist) {
    bssid_dist[i] ~ lognormal(
      log(
        sqrt(
        pow(bssid_full[bssid_id[i,1],1] - 
        bssid_full[bssid_id[i,2],1],2) +
        pow(bssid_full[bssid_id[i,1],2] - 
        bssid_full[bssid_id[i,2],2  ],2)
        )),
        sigma_mds);
  }
  //Model sampling statement -- bernoulli model with logit link function (equivalent to GLM with logit link)
  /*
  for(n in 1:N) {
  
  //first step -- latent MDS measurements
  //not indexed by n, but needs to run every n times to achieve a full joint
  //distribution over the latent positions for MDS

  
  /* Use accuracy to measure scale variance */
  
  /*
  
  sigma_lon ~ normal(acc_par_lon * acc[n],5);
  sigma_lat ~ normal(acc_par_lat * acc[n],5);
  
  */
  
  /* Model longitude and latitude using a latent stress model */
  
  /*
  
  latitude[n] ~ normal(loc_lat[cl[n]] + latitude_meas[n],sigma_lon);
  
  longitude[n] ~ normal(loc_lon[cl[n]] + longitude_meas[n],sigma_lat);
  
  
  
  
  distance_ft_sq[n] ~ normal(pow(bssid_full[bssid[n],2] - latitude_meas[n],2) +
                          pow(bssid_full[bssid[n],1] - longitude_meas[n],
                          2),
                        sigma_dist);
                        
  
  }
  */
}


generated quantities {
  
}
