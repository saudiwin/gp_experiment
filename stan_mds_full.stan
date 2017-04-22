
data {
  
  int N;
  
  //W=num wifi beacons
  int C;
  int G;
  int<lower=0> W;
  
  int wifi_id[N];
  int clust_id[N];
  
  //observed longitude/latitude
  
  int gps_id[N];
  vector[G] lon_meas;
  vector[G] lat_meas;
  //vector of observed dissimilarities == Euclidean distances
  
  vector[N] points_dist;
  
  //matrix[C,2] wifi_center;
  
  //matrix[C,2] wifi_sd;

}

parameters {
  
  //parameters to estimate, not observed
  //matrix[W,2] wifi_sd_raw;
  matrix[W,2] wifi_est;
  matrix[C,2] clust_par;
  matrix<lower=0>[C,2] wifi_sd;
  vector[G] latitude_raw;
  vector[G] longitude_raw;
  real<lower=0> sigma_mds;
  real lat_mu;
  real lon_mu;
  real wifi_lat_mu;
  real wifi_lon_mu;
  real<lower=0> lat_sigma;
  real<lower=0> lon_sigma;
  real<lower=0> lat_tau;
  real<lower=0> lon_tau;
  real<lower=0> wifi_lon_tau;
  real<lower=0> wifi_lat_tau;
    vector[G] latitude;
  vector[G] longitude;
}

transformed parameters {
  /*matrix[W,2] wifi_est;


  wifi_est[,1] = wifi_lat_mu + wifi_lat_tau*wifi_sd_raw[,1];
  wifi_est[,2] = wifi_lon_mu + wifi_lon_tau*wifi_sd_raw[,2];
  
  latitude = lat_mu + latitude_raw*lat_tau;
  longitude = lon_mu + longitude_raw*lon_tau;
  */
}

model {
  matrix[N,2] store_pts;
  vector[N] est_distance;
  //vague priors

  sigma_mds ~ normal(0,5);
  lat_sigma ~ normal(0,10);
  lon_sigma ~ normal(0,10);
  
  lat_mu ~ normal(0,10);
  lon_mu ~ normal(0,10);
  
  lat_tau ~ normal(0,10);
  lon_tau ~ normal(0,10);
  
  latitude_raw ~ normal(0,1);
  longitude_raw ~ normal(0,1);
  
  /*
  wifi_est[,1] ~ normal(0,lat_sigma);
  wifi_est[,2] ~ normal(0,lon_sigma);
  wifi_center[,1] ~ normal(clust_par[,1],wifi_sd[,1]);
  wifi_center[,2] ~ normal(clust_par[,2],wifi_sd[,2]);
  */
  /*to_vector(clust_par) ~ normal(0,10);
  to_vector(wifi_sd) ~ normal(0,10); */
  wifi_lon_tau ~ normal(0,5);
  wifi_lat_tau ~ normal(0,5);
  wifi_lon_mu ~ normal(0,10);
  wifi_lat_mu ~ normal(0,10);
  to_vector(wifi_est) ~ normal(0,10); 
  
  latitude ~ normal(lat_meas,.1);
  longitude ~ normal(lon_meas,.1);
  
  for(n in 1:N) {

      store_pts[n,1] = latitude[gps_id[n]];
      store_pts[n,2] = longitude[gps_id[n]];
      
      est_distance[n] = distance(wifi_est[wifi_id[n],] ,
      store_pts[n,]);
      
      print(est_distance[n]);
  }

  points_dist ~ exponential(est_distance);
  
}
