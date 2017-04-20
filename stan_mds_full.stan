
data {
  
  int N;
  
  //W=num wifi beacons
  int C;
  int<lower=0> W;
  
  int wifi_id[N];
  int clust_id[W];
  
  //observed longitude/latitude
  
  vector[N] lat_meas;
  vector[N] lon_meas;
  
  //vector of observed dissimilarities == Euclidean distances
  
  vector[N] points_dist;
  
  matrix[C,2] wifi_center;
  
  matrix[C,2] wifi_sd;

}

parameters {
  
  //parameters to estimate, not observed
  matrix[W,2] wifi_est;
  vector[N] latitude;
  vector[N] longitude;
  real<lower=0> sigma_mds;
  real dim1_mu;
  real dim2_mu;
  real<lower=0> dim1_sigma;
  real<lower=0> dim2_sigma;
}

model {
  
  //vague priors

  sigma_mds ~ normal(0,10);
  latitude ~ normal(lat_meas,1);
  longitude ~ normal(lon_meas,1);
  
    for(j in 1:W) {
  wifi_est[j,1] ~ normal(wifi_center[clust_id[j],1],wifi_sd[clust_id[j],1]);
  wifi_est[j,2] ~ normal(wifi_center[clust_id[j],2],wifi_sd[clust_id[j],2]);   
  }
  
  for(n in 1:N) {
    points_dist[n] ~ lognormal(
        log(sqrt(
        pow(wifi_est[wifi_id[n],1] - 
          latitude[n],2) +
        pow(wifi_est[wifi_id[n],2] - 
          longitude[n],2)
          )),
        sigma_mds);
  }
}
