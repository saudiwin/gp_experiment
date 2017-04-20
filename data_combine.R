require(fuzzyjoin)
require(dplyr)
require(readr)
require(tidyr)
require(ggplot2)
require(FactoMineR)


# Load Data ---------------------------------------------------------------



#make_data <- function(n_dist) {
load_files <- list.files('client_data/with_location/',full.names = TRUE)
all_files <- lapply(load_files,read_csv)
names(all_files) <- list.files('client_data/with_location/')
# Get rid of a row names column

all_files <- lapply(all_files, function(x)  {
  x$X1 <- NULL
  x$timestamp <- as.numeric(x$timestamp)
  return(x)})

# See how well the data overlap each other

ggplot(all_files$speech.csv,aes(x=timestamp)) + geom_density(fill='blue',alpha=0.5) + 
  geom_density(data=all_files$wifi.csv,fill='red',alpha=0.5) +
  geom_density(data=all_files$wifi.csv,fill='red',alpha=0.5) +
  geom_density(data=all_files$geo.csv,fill='green',alpha=0.5) +
  geom_density(data=all_files$accelerometer.csv,fill='pink',alpha=0.5)

check_data <- all_files$speech.csv
join_data <- all_files$speech.csv
speech <- all_files$speech

names(all_files) <- sapply(names(all_files),function(x) {
  x <- gsub('.csv',x=x,replacement='')
})

# Use a join on absolute difference between numeric timestamp values. Selects lowest distance as match.
# Seems to work pretty well

#Let's only use wifi and geo

all_files <- list(wifi=all_files$wifi,
                  geo=all_files$geo)

all_files <- lapply(all_files, function(x) {
  x <- difference_join(x,join_data,by='timestamp',mode='left',max_dist=1,
                      distance_col='distance')
  x$timestamp <- x$timestamp.x
  x <- fill(x,speech)
  return(x)
})

join_data <- difference_join(all_files$geo,all_files$wifi,by='timestamp',mode='left',
                       max_dist=3)

join_data <- select(join_data,reading_type.x.x,latitude,longitude,accuracy,
                    heading,speed,timestamp.x.x,timestamp.x.y,
                    reading_type.x.y,SSID,BSSID,
                    level,frequency,speech.x,speech.y,distance.x,distance.y) %>% rename(vartype_geo=reading_type.x.x,
                                                                  vartype_wifi=reading_type.x.y,
                                                                  speech_geo=speech.x,
                                                                  speech_wifi=speech.y,
                                                                  timestamp_geo=timestamp.x.x,
                                                                  timestamp_wifi=timestamp.x.y,
                                                                  distance_geo=distance.x,
                                                                  distance_wifi=distance.y)

join_data <- mutate(join_data,
                    combined_speech=coalesce(speech_geo,speech_wifi),
                    combined_dist=coalesce(distance_geo,distance_wifi)) %>% arrange(timestamp_geo) %>% 
  fill(combined_speech)

filter(join_data,combined_speech!='outdoors') %>% ggplot(aes(y=longitude,x=latitude,colour=combined_speech)) + geom_point() + theme_minimal() +
  stat_ellipse()

join_data <- filter(join_data,longitude<(-120.945))

for_dist <- select(join_data,BSSID,level,timestamp_wifi) %>% 
  mutate(unique_id=1:n()) %>% 
  distinct(BSSID,level,timestamp_wifi,.keep_all=TRUE) %>% 
  spread(key=BSSID,value=level)
to_dist <- lapply(select(for_dist,-timestamp_wifi,-unique_id),function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame()
to_dist$timestamp_wifi <- for_dist$timestamp_wifi
join_data <- left_join(join_data,to_dist,by='timestamp_wifi')
join_data <- filter(join_data,combined_speech!='outdoors')

# Join all then cluster observations

orig_labels <- factor(join_data$combined_speech)
#First need to separate labeled from unlabeled data

# unlabeled <- filter(join_data,is.na(orig_labels)) %>% dplyr::select(longitude,latitude,contains(':'))
labeled <- select(join_data,longitude,latitude,contains(':'))
known_labels <- orig_labels[!is.na(orig_labels)]
cluster_2 <- PCA(labeled,graph=FALSE)

# Now predict new clusters

# cluster_predict <- Rmixmod::mixmodPredict(data=select(unlabeled,longitude,latitude),cluster_2@bestResult)
# Rmixmod::plot(cluster_2)
# new_clust <- bind_rows(unlabeled,labeled)
# new_clust$clusters <- factor(c(factor(cluster_predict@partition,levels=unique(as.numeric(known_labels)),
#                                  labels=levels(known_labels)),known_labels),
#                              labels=levels(known_labels))
join_data$x <- cluster_2$ind$coord[,1]
join_data$y <- cluster_2$ind$coord[,2]
ggplot(join_data,aes(y=y,x=x,colour=combined_speech)) + geom_point() +
  stat_ellipse() + theme_minimal()

clusters <- HCPC(cluster_2,nb.clust=5)


# Stan set-up -------------------------------------------------------------



require(rstan)
require(mclust)
require(smacof)
require(proxy)

if(sample_it==TRUE) {
  join_data <- sample_n(join_data,500)
}

to_cluster <- data_frame(latitude=join_data$latitude,
                         longitude=join_data$longitude)

# Create weight matrix for missing data

fspl <- function(x_db,x_freq) {
  out_dist <- 10^((27.55 - (20 * log10(x_freq)) + abs(x_db))/20)
  return(out_dist)
}


for_dist<- select(join_data,BSSID,level,frequency,timestamp_geo) %>% filter(!is.na(level)) %>% 
  distinct(.keep_all=TRUE) %>% 
  mutate(distance=fspl(level,frequency)) %>% 
  select(-level,-frequency) %>% 
  complete(BSSID,timestamp_geo) %>% mutate(unique_id=1:n()) %>% 
  spread(key=BSSID,value=distance) %>% 
  group_by(timestamp_geo) %>% summarise_each(funs(sum(.,na.rm=TRUE))) %>% 
  mutate_each(funs(na_if(.,0))) %>% select(-unique_id)

# convert level and frequency back to feet

join_data$distance <- fspl(join_data$level,join_data$frequency)
join_data$time <- cut(join_data$timestamp_geo,100) %>% factor() %>% as.numeric()

# make a distance matrix

dist_matrix <- proxy::dist(x = select(for_dist,-timestamp_geo),
                           method='Euclidean',by_rows = F)

all_out_dist <- as.matrix(dist_matrix)

weight_dist <- apply(all_out_dist,2,function(x) {
  x <- if_else(is.na(x),0,1)
  x
})

collapsed_wifi <- smacofSym(delta=all_out_dist,weightmat=weight_dist,ndim=2)

std_err <- jackknife(collapsed_wifi)

std_err_dim1 <- apply(std_err$jackknife.conf[,1,],1,sd)
std_err_dim2 <- apply(std_err$jackknife.conf[,2,],1,sd)

# First identify what is latitude and what is longitude

# Need to transform using Procustes back to latitude/longitude

join_data <- mutate(join_data,
                    latitude_scale=as.numeric(scale(latitude)),
                    longitude_scale=as.numeric(scale(longitude)))

to_stan <- data_frame(BSSID=row.names(collapsed_wifi$conf),dim1=collapsed_wifi$conf[,1],
                       dim2=collapsed_wifi$conf[,2],
                       std_err_dim1=std_err_dim1,
                       std_err_dim2=std_err_dim2)
join_data <- left_join(join_data,to_stan,by='BSSID')
join_data <- mutate(join_data,BSSID=factor(BSSID),bssid_num=as.numeric(BSSID))

to_stan <- mutate(to_stan) %>% 
  mutate(BSSID=factor(BSSID,levels=levels(join_data$BSSID)),
         BSSID=as.numeric(BSSID))



sum_latlon <- group_by(join_data,BSSID) %>% summarize(lat=mean(latitude_scale[distance_wifi==min(distance_wifi,na.rm=TRUE)],na.rm=TRUE),
                                                      lon=mean(longitude_scale[distance_wifi==min(distance_wifi,na.rm=TRUE)],na.rm=TRUE)) %>%
  mutate(new_order=match(BSSID,row.names(collapsed_wifi$conf))) %>%
  arrange(new_order)

highest_pt <- mutate(sum_latlon,highest_pt=lat + lon) %>%
  filter(highest_pt==max(highest_pt))

lat_conv <- highest_pt$lat/collapsed_wifi$conf[highest_pt$BSSID,]



procrusted_one <- Procrustes(as.matrix(select(sum_latlon,lat,lon)),collapsed_wifi$conf)


# Merge back the collapsed values into join_data

collapsed_wifi$conf <- procrusted_one$Yhat

#check to see which dimensions correlate with which dimensions

cor(collapsed_wifi$conf[,1],sum_latlon$lat)

cor(collapsed_wifi$conf[,1],sum_latlon$lon)


# bayesian mds data prep --------------------------------------------------

long_dist <- gather(as_data_frame(all_out_dist),edge1,distance) %>% mutate(edge2=rep(row.names(all_out_dist),
                                                                                     ncol(all_out_dist)),
                                                                           edge1_fac=as.numeric(factor(edge1)),
                                                                           edge2_fac=as.numeric(factor(edge2))) %>% 
  filter(!is.na(distance))

# Need to find 3 points to fix with longitude/latitude

all_bssids <- group_by(join_data,BSSID) %>% summarize(length=sum(!is.na(longitude))+
                                                        sum(!is.na(level))) %>% arrange(desc(length)) %>% slice(1:3) %>% 
  left_join(long_dist,by=c('BSSID'='edge1')) %>% distinct(BSSID,length,edge1_fac)

fixed_bssids <- lapply(all_bssids$BSSID,function(x) {
  log_cond <- join_data$BSSID==x
  run_data <- filter(join_data,log_cond) %>% select(longitude_scale,latitude_scale,distance)
  model1 <- lm(formula=cbind(longitude_scale,latitude_scale)~distance + I(distance^2),data=run_data)
  max_vals <- predict(model1,newdata=data_frame(distance=0,`I(distance^2)`=0))
  return(max_vals)
}) %>% lapply(as_data_frame) %>% bind_rows()

fixed_bssids$BSSID <- all_bssids$BSSID

#relevel factor variables to put the fixed vals up front

join_data$BSSID <- fct_relevel(join_data$BSSID,fixed_bssids$BSSID)

#finally re-arrange distance matrix

long_dist <- mutate(long_dist,edge1=fct_relevel(factor(edge1),
                                                fixed_bssids$BSSID),
                    edge1_fac=as.numeric(edge1),
                    edge2=fct_relevel(factor(edge2),
                                      fixed_bssids$BSSID),
                    edge2_fac=as.numeric(edge2)) %>% 
  arrange(edge1,edge2) %>% 
  filter(distance!=0) %>% 
  mutate(distance=log(distance))

# Set clusters equal to BSSID networks

join_data <- mutate(join_data,
                    cluster=substr(as.character(BSSID),1,9))
  

# stan data prep ----------------------------------------------------------

model_frame_data <- model.frame(cbind(longitude_scale,latitude_scale)~scale(accuracy) + time +scale(distance) + combined_speech + bssid_num + dim1 + dim2 +std_err_dim1 +std_err_dim2,data=join_data)
#model_frame_data <- sample_n(model_frame_data,500)
bssid <- model_frame_data$bssid_num
acc <- model.matrix(cbind(longitude_scale,latitude_scale)~scale(accuracy),data=model_frame_data)[,-1]
distance <- as.numeric(model_frame_data$`scale(distance)`)
#Get rid of the  intercept because variables are standardized
model_response <-  model.response(model_frame_data)

# keep_model <- apply(model_data,2,function(x) sum(x==1))[-c(1,2,ncol(model_data))]
# keep_model <- keep_model>2
# keep_model <- c(TRUE,TRUE,keep_model,TRUE)
# model_data <- model_data[,keep_model]





stan_data <- list(N=nrow(model_frame_data),
                  J=length(unique(join_data$bssid_num)),
                  bssid=bssid,
                  acc=as.numeric(acc)*-1,
                  num_dist=nrow(long_dist),
                  bssid_dist=long_dist$distance,
                  bssid_id=select(long_dist,edge1_fac,edge2_fac),
                  bssid_fix_num=nrow(fixed_bssids),
                  bssid_fix_vals=fixed_bssids,
                  longitude=model_response[,1],
                  latitude=model_response[,2],
                  distance_ft=distance,
                  G=length(unique(join_data$cluster)),
                  cl=as.integer(factor(join_data$cluster)))

# Run the Stan sampler

niter <- 2000
nwarmup <- 1000

# Bayesian MDS model ------------------------------------------------------





stan_circle <- stan_model(file='spatial_radius.stan',model_name='spatial_radius')

stan_run_radius <- sampling(object = stan_circle,data=stan_data,chains=2,iter=niter,warmup=nwarmup,cores=2)
  
saveRDS(stan_run_radius,'stan_circle_latest.rds')

compiled_stan <- stan_model(file='spatial_measerr.stan',model_name='spatial_measerr')

stan_run_measerr <- sampling(object = compiled_stan,data=stan_data,chains=2,iter=niter,warmup=nwarmup,cores=2)

saveRDS(object = stan_run_measerr,file = 'most_recent_stan.rds')

#Pull out the posterior estimates from the model

stan_run_measerr <- readRDS('most_recent_stan.rds')

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
