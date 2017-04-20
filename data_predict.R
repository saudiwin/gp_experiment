# Predict values from stan objects

stan_obj <- rstan::extract(readRDS('stan_circle_latest.rds'),permuted=TRUE)

posteriors <- rstan::extract(stan_run_measerr,permuted=TRUE)

# pull generated quantities and call loo

lat_pred <- matrix(nrow=niter-nwarmup,ncol=stan_data$N)
lon_pred <- matrix(nrow=niter-nwarmup,ncol=stan_data$N)

for(i in 1:(niter-nwarmup)) {
  for(n in 1:stan_data$N) {
    bssid <- stan_data$bssid
    lat_pred[i,n] <- rnorm(n=1,mean=posteriors$loc_lat[i,stan_data$cl[n]] + 
                             posteriors$scale1_lat[i]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$scale2_lat[i]*posteriors$bssid_dim_par[i,bssid[n],2] +
                             posteriors$med_int[i,1,1]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$med_int[i,2,1]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],2] +  
                             posteriors$dist_par_lat[i]*stan_data$distance_ft[n] +
                             posteriors$lat_int[i]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],1]*posteriors$bssid_dim_par[i,bssid[n],2],
                           sd=posteriors$sigma_lat[i])
    lon_pred[i,n] <- rnorm(n=1,mean=posteriors$loc_lon[i,stan_data$cl[n]] + 
                             posteriors$scale1_lon[i]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$scale2_lon[i]*posteriors$bssid_dim_par[i,bssid[n],2] +
                             posteriors$dist_par_lon[i]*stan_data$distance_ft[n] +
                             posteriors$med_int[i,1,2]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$med_int[i,2,2]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],2] + 
                             posteriors$lon_int[i]*stan_data$distance_ft[n]*posteriors$bssid_dim_par[i,bssid[n],1]*posteriors$bssid_dim_par[i,bssid[n],2],
                           sd=posteriors$sigma_lon[i])
  } 
}
predictions <- data_frame(lat_mean=apply(lat_pred,2,mean),
                          lon_mean=apply(lon_pred,2,mean),
                          lat_sd=apply(lat_pred,2,sd),
                          lon_sd=apply(lon_pred,2,sd))

# Plot the predictive density with the true outcome super-imposed
require(bayesplot)
ppc_dens_overlay(stan_data$latitude,lat_pred[1:50,])
ppc_dens_overlay(stan_data$longitude,lon_pred[1:50,])

for_bayes <- as.array(stan_run_measerr)
mcmc_trace(for_bayes,regex_pars='loc_lat.*')
mcmc_intervals(for_bayes,regex_pars='int.*')
mcmc_intervals(for_bayes,regex_pars='bssid_dim_par.*')
mcmc_trace(for_bayes,regex_pars='err.*')
mcmc_trace(for_bayes,regex_pars='acc.*')
mcmc_trace(for_bayes,regex_pars='sigma.*')
mcmc_intervals(for_bayes,regex_pars='scale.*')
mcmc_intervals(for_bayes,regex_pars='dist.*')

all_predict_lat <- as_data_frame(lat_pred) %>% gather(data_pt,latitude) %>% mutate(data_pt=stringr::str_extract(data_pt,'[0-9]+'),
                                                                                   data_pt=as.numeric(data_pt))
all_predict_lon <- as_data_frame(lon_pred) %>% gather(data_pt,longitude)
all_predict <- mutate(all_predict_lat,longitude=all_predict_lon$longitude)
join_data <- mutate(join_data,data_pt=1:n())
all_predict <- left_join(all_predict,select(join_data,data_pt,time))
predictions$time <- join_data$time

true_speech <- filter(join_data,combined_dist<1) %>% select(latitude_scale,longitude_scale,
                                                            combined_speech) %>% 
  group_by(combined_speech) %>% summarize(avg_lat=mean(latitude_scale),
                                          avg_lon=mean(longitude_scale))

ggplot(predictions) + geom_point(aes(y=lat_mean,x=lon_mean,colour=time)) +
  geom_point(data=join_data,aes(y=latitude_scale,x=longitude_scale),colour='red',alpha=0.5) +
  geom_text(data=true_speech,aes(y=avg_lat,x=avg_lon,label=combined_speech),fontface='bold') +
  theme_minimal() + ylab('latitude') + scale_color_distiller() +
  xlab('longitude')

ggsave('point_overlay.png',width=8,units='in',scale=1.1)

speech_path <- filter(join_data,combined_dist<1) %>% select(latitude_scale,longitude_scale,
                                                            combined_speech,
                                                            timestamp_geo) %>% distinct(.keep_all=TRUE) %>% 
  mutate(timestamp_geo=as.POSIXct(timestamp_geo,tz='GMT',origin='1970-01-01 00:00.00 UTC'))



speech_predict <- filter(join_data,combined_dist<1) %>% select(time,bssid_num,distance,
                                                               timestamp_geo) %>% 
  mutate(distance=as.numeric(scale(distance)),
         timestamp_geo=as.POSIXct(timestamp_geo,tz='GMT',origin='1970-01-01 00:00.00 UTC')) %>% 
  distinct(.keep_all=TRUE)

lat_pred <- matrix(nrow=100,ncol=nrow(speech_predict))
lon_pred <- matrix(nrow=100,ncol=nrow(speech_predict))

for(i in 1:100) {
  for(n in 1:nrow(speech_predict)) {
    bssid <- speech_predict$bssid_num
    distance <- speech_predict$distance
    time <- speech_predict$time
    lat_pred[i,n] <- rnorm(n=1,mean=posteriors$loc_lat[i,time[n]] + 
                             posteriors$scale1_lat[i]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$scale2_lat[i]*posteriors$bssid_dim_par[i,bssid[n],2] +
                             posteriors$med_int[i,1,1]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$med_int[i,2,1]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],2] +  
                             posteriors$dist_par_lat[i]*distance[n] +
                             posteriors$lat_int[i]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],1]*posteriors$bssid_dim_par[i,bssid[n],2],
                           sd=posteriors$sigma_lat[i])
    lon_pred[i,n] <- rnorm(n=1,mean=posteriors$loc_lon[i,time[n]] + 
                             posteriors$scale1_lon[i]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$scale2_lon[i]*posteriors$bssid_dim_par[i,bssid[n],2] +
                             posteriors$dist_par_lon[i]*distance[n] +
                             posteriors$med_int[i,1,2]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],1] +
                             posteriors$med_int[i,2,2]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],2] + 
                             posteriors$lon_int[i]*distance[n]*posteriors$bssid_dim_par[i,bssid[n],1]*posteriors$bssid_dim_par[i,bssid[n],2],
                           sd=posteriors$sigma_lon[i])
  } 
}
speech_predict <- mutate(speech_predict,lat_pred=apply(lat_pred,2,mean),
                         lon_pred=apply(lon_pred,2,mean))

tween_predict <- data_frame(timestamp_geo=tween_datetime(speech_predict$timestamp_geo,n=20,ease='quadratic-in-out')[[1]],
                            latitude_scale=tween_numeric(speech_predict$lat_pred,n=20,ease='quadratic-in-out')[[1]],
                            longitude_scale=tween_numeric(speech_predict$lon_pred,n=20,ease='quadratic-in-out')[[1]]) %>% 
  left_join(speech_predict) %>% mutate(frame=1:n())


combined_speech <- left_join(speech_path,speech_predict,by='timestamp_geo') %>% group_by(timestamp_geo) %>% 
  mutate(lat_pred_mean=mean(lat_pred),
         lon_pred_mean=mean(lon_pred)) %>% 
  distinct(timestamp_geo,.keep_all=T)


tween_speech <- data_frame(timestamp_geo=tween_datetime(combined_speech$timestamp_geo,n=20,ease='quadratic-in-out')[[1]],
                           latitude_scale=tween_numeric(combined_speech$latitude_scale,n=20,ease='quadratic-in-out')[[1]],
                           longitude_scale=tween_numeric(combined_speech$longitude_scale,n=20,ease='quadratic-in-out')[[1]],
                           latitude_mean=tween_numeric(combined_speech$lat_pred_mean,n=20,ease='quadratic-in-out')[[1]],
                           longitude_mean=tween_numeric(combined_speech$lon_pred_mean,n=20,ease='quadratic-in-out')[[1]]) %>% 
  left_join(speech_path) %>% fill(combined_speech) %>% mutate(frame=1:n())

tween_speech <- gather(tween_speech,type,estimate,longitude_scale,latitude_scale,longitude_mean,latitude_mean) %>% 
  separate(type,into = c('var','origin')) %>% spread(key=var,value=estimate)

# animate
require(gganimate)
require(tweenr)


ani_plot <- ggplot(tween_speech,aes(y=latitude,x=longitude)) + geom_path(aes(frame=frame,cumulative=TRUE,colour=origin,linetype=origin)) +
  geom_text(data=true_speech,aes(y=avg_lat,x=avg_lon,label=combined_speech),fontface='bold') +
  theme_minimal() + ylab('latitude') +
  xlab('longitude')
gganimate(ani_plot,interval=.1,'follow_path.mp4')

# Create test data using the original speech dataset

test_data <- difference_join(speech,join_data,by=c('timestamp'='timestamp_geo'))

#predict_out <- 
