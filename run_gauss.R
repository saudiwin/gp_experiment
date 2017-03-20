# Gaussian process regression


require(readr)
require(dplyr)
require(rstan)
require(bayesplot)
require(ggplot2)
require(mclust)
gps_data <- read_csv('gps_u00.csv')
wifi_location <- read_csv('wifi_location_u00.csv')
wifi <- read_csv('wifi_u00.csv')
blue_tooth <- read_csv('bt_u00.csv')

normalize <- function(x) {
  y <- (x-min(x))/(max(x)-min(x))
  y[y==0] <- 0.001
  y[y==1] <- .999
  return(y)
}

change_time <- function(dataset) {
  dataset$time <- as.character(dataset$time) %>% substr(start=1,stop=7) %>% as.numeric
  return(dataset)
}
gps_data <- change_time(gps_data)
blue_tooth <- change_time(blue_tooth)
wifi <- change_time(wifi)
join_data <- left_join(x = gps_data,y=blue_tooth,by='time') %>% left_join(y=wifi,by='time')

# Too much data, just keep 100,000 rows
join_data <- sample_n(join_data,20000)


# run mclust on latitude/longitude to find distinct geographical areas
to_cluster <- data_frame(latitude=join_data$latitude,
                         longitude=join_data$longitude)

mclust_result <- densityMclust(to_cluster)
groups <- data_frame(clusters=mclust_result$classification) %>% left_join(data_frame(clusters=1:mclust_result$G,
                                                                                     proportion=mclust_result$parameters$pro),
                                                                          by='clusters') %>% 
  mutate(to_sample=(1/proportion))
join_data$clusters <- groups$clusters
# Sample the join_data so that it balances the internal clusters better

to_model <- sample_n(join_data,500,weight = groups$to_sample)

#basic_model <- lm(cbind(latitude,longitude)~accuracy + altitude + bearing + travelstate + MAC,data=join_data)

# Predictive performance seems OK
#apply(basic_model$residuals,2,function(x) mean(x^2))

# Let's set this up in Stan
model_frame_data <- model.frame(cbind(longitude,latitude)~scale(accuracy) +  clusters +travelstate + BSSID + scale(level.y),data=to_model)
model_data <- model.matrix(cbind(longitude,latitude)~scale(accuracy) +   travelstate + BSSID + scale(level.y),data=model_frame_data)
#Get rid of the  intercept because variables are standardized
model_data <- model_data[,-1]
model_response <-  model.response(model_frame_data)

keep_model <- apply(model_data,2,function(x) sum(x==1))[-c(1,2,ncol(model_data))]
keep_model <- keep_model>2
keep_model <- c(TRUE,TRUE,keep_model,TRUE)
model_data <- model_data[,keep_model]

stan_data <- list(N=nrow(model_data),
                  J=ncol(model_data),
                  X=model_data,
                  longitude=scale(model_response[,'longitude'])[,1],
                  latitude=scale(model_response[,'latitude'])[,1])

# To make it more interesting, drop all wifi posts that only have one location attached (i.e., user only visited one time)

compiled_stan <- stan_model(file='spatial_predict.stan',model_name='spatial_location_regression')

stan_run <- sampling(object = compiled_stan,data=stan_data,chains=2,iter=1000,warmup=800,cores=2)
posteriors <- rstan::extract(stan_run,permuted=TRUE)

# pull generated quantities and call loo

lat_pred <- matrix(nrow=400,ncol=stan_data$N)
lon_pred <- matrix(nrow=400,ncol=stan_data$N)

for(i in 1:400) {
  for(n in 1:stan_data$N) {
    lat_pred[i,n] <- posteriors$theta_lat[i,]%*%stan_data$X[n,]
    lon_pred[i,n] <- posteriors$theta_lon[i,]%*%stan_data$X[n,]
  } 
}
predictions <- data_frame(lat_mean=apply(lat_pred,2,mean),
                          lon_mean=apply(lon_pred,2,mean),
                          lat_sd=apply(lat_pred,2,sd),
                          lon_sd=apply(lon_pred,2,sd))
# Look at just greater than 0
ppc_dens_overlay(stan_data$latitude,lat_pred[1:50,])
ppc_dens_overlay(stan_data$longitude,lon_pred[1:50,])


# Add in varying group means according to clusters we have identified


stan_data <- list(N=nrow(model_data),
                  J=ncol(model_data),
                  X=model_data,
                  longitude=scale(model_response[,'longitude'])[,1],
                  latitude=scale(model_response[,'latitude'])[,1],
                  G=mclust_result$G,
                  cl=model_frame_data$clusters)

compiled_stan <- stan_model(file='spatial_multgroups.stan',model_name='spatial_multgroups')


stan_run_multgroups <- sampling(object = compiled_stan,data=stan_data,chains=2,iter=1000,warmup=800,cores=2)
posteriors_groups <- rstan::extract(stan_run_multgroups,permuted=TRUE)

posteriors <- rstan::extract(stan_run_multgroups,permuted=TRUE)

# pull generated quantities and call loo

lat_pred <- matrix(nrow=400,ncol=stan_data$N)
lon_pred <- matrix(nrow=400,ncol=stan_data$N)

for(i in 1:400) {
  for(n in 1:stan_data$N) {
    lat_pred[i,n] <- posteriors$loc_lat[i,stan_data$cl[n]]
    lon_pred[i,n] <- posteriors$loc_lon[i,stan_data$cl[n]]
  } 
}
predictions <- data_frame(lat_mean=apply(lat_pred,2,mean),
                          lon_mean=apply(lon_pred,2,mean),
                          lat_sd=apply(lat_pred,2,sd),
                          lon_sd=apply(lon_pred,2,sd))
# Look at just greater than 0
ppc_dens_overlay(stan_data$latitude,lat_pred[1:50,])
ppc_dens_overlay(stan_data$longitude,lon_pred[1:50,])

for_bayes <- as.array(stan_run_multgroups)
mcmc_trace(for_bayes,regex_pars='loc_lat.*')
mcmc_intervals(for_bayes,regex_pars='theta_lat.*')
mcmc_intervals(for_bayes,regex_pars='theta_lon.*')
mcmc_trace(for_bayes,regex_pars='tau.*')
mcmc_trace(for_bayes,regex_pars='mu.*')
