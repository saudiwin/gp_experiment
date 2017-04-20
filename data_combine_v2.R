require(fuzzyjoin)
require(dplyr)
require(readr)
require(tidyr)
require(ggplot2)
require(FactoMineR)
require(forcats)


# model parameters --------------------------------------------------------

sample_it <- FALSE
niter <- 2000
nwarmup <- 1000
simulate <- TRUE

# Load Data ---------------------------------------------------------------

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

# Stan set-up -------------------------------------------------------------



require(rstan)
require(mclust)
require(smacof)
require(proxy)

if(sample_it==TRUE) {
  join_data <- sample_n(join_data,500)
}

# scale longitude / latitude before doing any calculations

join_data <- mutate(join_data,
                    longitude_scale=as.numeric(scale(longitude)),
                    latitude_scale=as.numeric(scale(latitude)))

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

distance_with_nas <- function(x,y) {
  # only compute on values without NAs
  all_x <- !is.na(x)
  all_y <- !is.na(y)
  both <- all_x & all_y
  total <- sum(both)
  eu_dist <- sqrt(sum((x[both] - y[both])^2))/total
  return(eu_dist)
}

if(simulate==TRUE) {
  # for_dist <- matrix(nrow=length(unique(join_data$BSSID)),
  #                    ncol=length(unique(join_data$BSSID)))
  # for_dist <- apply(for_dist,2,function(x) {
  #   x <- rlnorm(n=length(x),runif(n = 1,min = 3,max=6),1)
  # })
  # row.names(for_dist) <- levels(join_data$BSSID)
  # colnames(for_dist) <- levels(join_data$BSSID)
  # dist_matrix <- proxy::dist(x = for_dist,
  #                            method=distance_with_nas,by_rows = F)
  n_sim <- length(unique(join_data$BSSID))
  true_points <- matrix(nrow=n_sim,c(rnorm(n_sim,0,10),rnorm(n_sim,0,10)))
  
  dist_matrix <- proxy::dist(true_points) %>% as.matrix() %>% as_data_frame()
  
  all_out_dist <- as.matrix(dist_matrix)
  row.names(all_out_dist) <- levels(factor(join_data$BSSID))
  colnames(all_out_dist) <- levels(factor(join_data$BSSID))
} else {

dist_matrix <- proxy::dist(x = select(for_dist,-timestamp_geo),
                           method=distance_with_nas,by_rows = F)
all_out_dist <- as.matrix(dist_matrix)
}


# Compute starting values 

dist_start <- smacof::smacofSym(all_out_dist)
dist_start <- dist_start$conf

# bayesian mds data prep --------------------------------------------------

long_dist <- gather(as_data_frame(all_out_dist),edge1,distance) %>% mutate(edge2=rep(row.names(all_out_dist),
                                                                                     ncol(all_out_dist)),
                                                                           edge1_fac=as.numeric(factor(edge1)),
                                                                           edge2_fac=as.numeric(factor(edge2))) %>% 
  filter(!is.na(distance)) %>% distinct(edge1,edge2,distance,.keep_all=T)

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

join_data$bssid_num <- as.integer(join_data$BSSID)

#finally re-arrange distance matrix

long_dist <- mutate(long_dist,edge1=fct_relevel(factor(edge1),
                                                fixed_bssids$BSSID),
                    edge1_fac=as.numeric(edge1),
                    edge2=fct_relevel(factor(edge2),
                                      fixed_bssids$BSSID),
                    edge2_fac=as.numeric(edge2)) %>% 
  arrange(edge1,edge2) %>% filter(distance!=0)

dist_start <- dist_start[levels(join_data$BSSID),]

# Need to remove duplicate pairs -- only need to calculate distance once
max_edge <- max(as.numeric(long_dist$edge1))
long_dist <- group_by(long_dist,edge1_fac) %>% mutate(unnecessary=if_else((edge1_fac>edge2_fac),T,F)) %>% 
  filter(!unnecessary) %>% ungroup

# Set clusters equal to BSSID networks

join_data <- mutate(join_data,
                    cluster=substr(as.character(BSSID),1,9))
  

# stan data prep ----------------------------------------------------------

model_frame_data <- model.frame(cbind(longitude_scale,latitude_scale)~scale(accuracy) + time +scale(distance) + combined_speech + bssid_num,data=join_data)
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
                  bssid_fix_vals=select(fixed_bssids,longitude_scale,latitude_scale),
                  longitude=model_response[,1],
                  latitude=model_response[,2],
                  distance_ft=distance,
                  G=length(unique(join_data$cluster)),
                  cl=as.integer(factor(join_data$cluster)))

# Bayesian MDS model ------------------------------------------------------


inits <- list(bssid_dim_par_free=dist_start[(stan_data$bssid_fix_num+1):stan_data$J,],
                       bssid_dim_par_fixed=select(fixed_bssids,longitude_scale,latitude_scale))
inits <- list(inits,inits)


stan_circle <- stan_model(file='spatial_radius.stan',model_name='spatial_radius')

stan_run_radius <- sampling(object = stan_circle,data=stan_data,chains=2,iter=niter,warmup=nwarmup,cores=2,
                            control=list(adapt_delta=0.95),
                            init=inits)
  
saveRDS(stan_run_radius,'stan_circle_latest.rds')

# compiled_stan <- stan_model(file='spatial_measerr.stan',model_name='spatial_measerr')
# 
# stan_run_measerr <- sampling(object = compiled_stan,data=stan_data,chains=2,iter=niter,warmup=nwarmup,cores=2)
# 
# saveRDS(object = stan_run_measerr,file = 'most_recent_stan.rds')

#Pull out the posterior estimates from the model



