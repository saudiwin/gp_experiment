
data {
  
  //J= num points
  
  int<lower=0> J;
  
  //number of distances
  
  int num_dist;
  
  //points fixed to true values
  
  int points_fix_num;
  
  //point pairs
  
  int points_id[num_dist,2];
  
  //vector of observed dissimilarities == Euclidean distances
  
  vector[num_dist] points_dist;
  
  matrix[points_fix_num,2] points_fix_vals;

}

transformed data {
  vector[num_dist] dist_std;
  real dist_sd;
  real dist_mean;
  
  dist_sd = sd(points_dist);
  dist_mean = mean(points_dist);
  
  dist_std = (points_dist - dist_mean)./dist_sd;
}


parameters {
  
  //parameters to estimate, not observed
  matrix[J-points_fix_num,2] points_dim_par_free;
  matrix[points_fix_num,2] points_dim_par_fixed;

  vector<lower=0>[J-points_fix_num] sigma_dim1_raw;
  vector<lower=0>[J-points_fix_num] sigma_dim2_raw;
  real<lower=0> sigma_mds;
  real dim1_mu;
  real dim2_mu;
  real<lower=0> dim1_sigma;
  real<lower=0> dim2_sigma;
}

transformed parameters {
  //create combined distance coefficients/locations matrix
  matrix[J,2] points_full;
  
  //use a non-centered hierarchical prior for point variances
  //for some reason, point variances are difficult to compute in this model
  vector<lower=0>[J-points_fix_num] sigma_dim1_nc;
  vector<lower=0>[J-points_fix_num] sigma_dim2_nc;
  
  sigma_dim1_nc = dim1_mu + sigma_dim1_raw*dim1_sigma;
  sigma_dim2_nc = dim2_mu + sigma_dim2_raw*dim2_sigma;
  
  points_full = append_row(points_dim_par_fixed,points_dim_par_free);
  
  //center coefficients around fixed point
  //divide by magnitude of second point
  
  for(j in 1:J) {
    for(d in 1:2) {
      points_full[j,d] = (points_full[j,d] - points_full[1,d])/((points_full[1,d] - points_full[2,d])^2);
    }
  }
  
}


model {
  
  //vague priors

  sigma_mds ~ normal(0,5);
  dim1_mu ~ normal(0,5);
  dim2_mu ~ normal(0,5);
  dim1_sigma ~ normal(0,5);
  dim2_sigma ~ normal(0,5);
  sigma_dim1_raw ~ normal(0,1);
  sigma_dim2_raw ~ normal(0,1);
  
  for(j in 1:(J-points_fix_num)) {
  points_dim_par_free[j,1] ~ normal(0,sigma_dim1_nc[j]);
  points_dim_par_free[j,2] ~ normal(0,sigma_dim2_nc[j]);   
  }
  //fix to true values
  to_vector(points_dim_par_fixed) ~ normal(to_vector(points_fix_vals),.01); 
  
  for(i in 1:num_dist) {
    points_dist[i] ~ lognormal(
        log(sqrt(
        pow(points_full[points_id[i,1],1] - 
          points_full[points_id[i,2],1],2) +
        pow(points_full[points_id[i,1],2] - 
          points_full[points_id[i,2],2  ],2)
          )),
        sigma_mds);
  }
}

generated quantities {
   matrix[J,2] points_unstd;
  for(j in 1:J) {
    for(d in 1:2) {
      points_unstd[j,d] = (points_full[j,d]*dist_sd)+dist_mean;
    }
  }
}
