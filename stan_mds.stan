functions {
 vector rescale(vector X, int N) {
   vector[N] X_trans;
   for(n in 1:N) {
     X_trans[n] = (2*pi())/((max(X) - min(X)) * (X[n] - max(X))) + 2*pi();
   }
   return X_trans;
 } 
  
}

data {
  
  //J= num points
  
  int<lower=0> J;
  
  //number of distances
  
  int N;
  
  int polar;
  
  //points fixed to true values
  
  int points_fix_num;
  
  //point pairs
  
  int points_id[N,2];
  
  //vector of observed dissimilarities == Euclidean distances
  
  vector[N] points_dist;
  
  matrix[points_fix_num,2] points_fix_vals;

}

transformed data {
  matrix[points_fix_num,2] points_fix_trans;
  vector[N] points_dist_trans;
  real magnitude;
  int pos_or_neg;
  real hypot_pts;
  vector[3] triangle_dist;
  
  for(i in 1:2) {
    triangle_dist[i] = points_dist[i];
  }
  for(n in 1:N) {
    if(((points_id[n,1]==2) && (points_id[n,2]==3)) || ((points_id[n,2]==3) && (points_id[n,1]==2))) {
          triangle_dist[3] = points_dist[n];
    }
  }
  hypot_pts = (triangle_dist[1]^2 + triangle_dist[2]^2) - triangle_dist[3];
  
  points_fix_trans[1,1] = 0.0;
  points_fix_trans[1,2] = 0.0;
  
  //put distances to fixed origin at zero
  

  //center all fixed points as unit vectors around the first fixed point (the origin)
  
  if(points_fix_num>1) {
    if(points_fix_num==3) {
      
      points_fix_trans[2,1] = 0.0;
      //points_fix_trans[3,2] = 0.0;
      //now assign distances to diagonal points, need to loop to find them
      /*
      for(n in 1:N) {
        if(((points_id[n,1]==1) && (points_id[n,2]==2)) || ((points_id[n,2]==1) && (points_id[n,1]==2))) {
          points_fix_trans[2,2] = points_dist[n];
        } else if(((points_id[n,1]==1) && (points_id[n,2]==3)) || ((points_id[n,2]==1) && (points_id[n,1]==3))) {
          //points_fix_trans[3,1] = points_dist[n];
        }
      }
      */
      
    } else if(points_fix_num==2) {
      points_fix_trans[2,1] = 0.0;
      /*
        for(n in 1:N) {
          if(((points_id[n,1]==1) && (points_id[n,2]==2)) || ((points_id[n,2]==1) && (points_id[n,1]==2))) {
          points_fix_trans[2,2] = points_dist[n];
          }
        }
        */
    } else {
      //if more than 3 points are fixed, then just keep them as they are == extra identification
      
      points_fix_trans[2,1] = 0.0;
      //points_fix_trans[3,2] = 0.0;
      //now assign distances to diagonal points, need to loop to find them
      /*
      for(n in 1:N) {
        if(((points_id[n,1]==1) && (points_id[n,2]==2)) || ((points_id[n,2]==1) && (points_id[n,1]==2))) {
          points_fix_trans[2,2] = points_dist[n];
        } else if(((points_id[n,1]==1) && (points_id[n,2]==3)) || ((points_id[n,2]==1) && (points_id[n,1]==3)))
          points_fix_trans[3,1] = points_dist[n];
      } */
      for(i in 4:points_fix_num) {
        for(d in 1:2) {
          points_fix_trans[i,d] = points_fix_vals[i,d];  
        }
      }
      
    }
  }
  if(polar==1) {
    for(n in 1:N) {
      points_dist_trans = rescale(points_dist,N);
    }
  }
}


parameters {
  
  //parameters to estimate, not observed
  
  //matrix[J-points_fix_num,2] p_free;
  
  matrix[points_fix_num,2] p_fixed;

  vector<lower=0>[J-points_fix_num] sigma_dim1_raw;
  vector<lower=0>[J-points_fix_num] sigma_dim2_raw;
  real<lower=0> sigma_mds;
  real dim1_mu;
  real dim2_mu;
  real three_float;
  real<lower=0> three_constrain_pos;
  real<upper=0> three_constrain_neg;
  real two_float;
  real<lower=0> two_constrain;
  real<lower=0> dim1_sigma;
  real<lower=0> dim2_sigma;
}

transformed parameters {
  //create combined distance coefficients/locations matrix
  matrix[J,2] points_full;
  matrix[points_fix_num,2] p_fixed_trans;
  matrix[J-points_fix_num,2] p_free;
  
  p_fixed_trans[1,] = points_fix_trans[1,];
  
  if(points_fix_num==2) {
  p_fixed_trans[2,1] = 0.0;
  p_fixed_trans[2,2] = two_float;
  }
  if(points_fix_num==3) {
  p_fixed_trans[2,1] = 0.0;
  p_fixed_trans[2,2] = two_float;
  p_fixed_trans[3,1] = three_constrain;
  p_fixed_trans[3,2] = three_float;
  }
  
  //matrix[N,4] m_points_full;
  
  //use a non-centered hierarchical prior for point variances
  //for some reason, point variances are difficult to compute in this model
  
  
  p_free[,1] = dim1_mu + sigma_dim1_raw*dim1_sigma;
  p_free[,2] = dim2_mu + sigma_dim2_raw*dim2_sigma; 
  
  points_full = append_row(p_fixed,p_free);
  
}


model {
  vector[N] store_pts;
  //vague priors

  sigma_mds ~ normal(0,5);
  dim1_mu ~ normal(0,10);
  dim2_mu ~ normal(0,10);
  dim1_sigma ~ normal(0,20);
  dim2_sigma ~ normal(0,20);
  sigma_dim1_raw ~ normal(0,1);
  sigma_dim2_raw ~ normal(0,1);
  three_float ~ normal(0,20);
  three_constrain ~ normal(0,20);
  two_float ~ normal(0,20);
  
  //fix to true values
  to_vector(p_fixed[1:2,]) ~ normal(to_vector(points_fix_trans[1:2,]),.01); 
  if(polar==0) {
    
    for(n in 1:N) {
      store_pts[n] = distance(points_full[points_id[n,1],], 
                            points_full[points_id[n,2],]);
    } 
      points_dist ~ lognormal(
          log(store_pts),
          sigma_mds);
  } else {
  for(n in 1:N) {
    store_pts[n] = acos(
                        dot_product(points_full[points_id[n,1],], 
                                    points_full[points_id[n,2],])/
                            (sqrt(dot_self(points_full[points_id[n,1],]))*
                            sqrt(dot_self(points_full[points_id[n,2],]))));
  }
  points_dist_trans ~ von_mises(store_pts,sigma_mds);
}
}
