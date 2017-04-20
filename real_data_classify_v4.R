
require(readr)
require(dplyr)
require(ggplot2)
require(magrittr)
require(mirt)
require(tidyr)
require(FactoMineR)


# Let's redo the data 

source('data_combine.R')

final_data <- mutate(join_data,
                     timestamp=timestamp_geo,
                     speech=labels)

# Remove BSSID's with very few values

wifi_data <- select(final_data,level,BSSID,timestamp,speech) %>% group_by(BSSID)

# Try MDS
# Create distance matrix
# As with IRT, missing values are allowed

for_dist <- select(wifi_data,speech,BSSID,level,timestamp) %>% filter(!is.na(level)) %>% 
  distinct(.keep_all=TRUE) %>% 
  spread(key=BSSID,value=level)
to_dist <- lapply(select(for_dist,-timestamp),function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame()
to_dist$timestamp <- for_dist$timestamp
to_dist <- select(to_dist,timestamp,everything()) %>% dist()
mds_model <- cmdscale(to_dist)

x <- mds_model[,1]
y <- mds_model[,2]

data_frame(x=x,y=y,labels=for_dist$timestamp) %>% ggplot(aes(x=x,y=y)) + geom_text(aes(label=labels),check_overlap = TRUE) +
  theme_minimal()

#Now PCA
# for_dist <- select(wifi_data,BSSID,level,timestamp,speech) %>% filter(!is.na(level)) %>% group_by(BSSID,timestamp) %>% 
#   summarize(avg_level=mean(level,na.rm=TRUE)) %>% filter(!is.na(avg_level)) %>% spread(key=BSSID,value=avg_level)

to_dist <- lapply(select(for_dist,-timestamp),function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame()

to_dist$speech <- for_dist$speech
to_dist$timestamp <- for_dist$timestamp
to_dist <- select(to_dist,timestamp,speech,everything())

# Join with GPS data

gps_data <- select(final_data,timestamp,latitude,longitude,speech)
combined <- left_join(gps_data,to_dist,by='timestamp') %>% 
  mutate(speech=speech.x) %>% select(-speech.y,-speech.x)
combined <- select(combined,speech,everything())
pca_model <- PCA(combined,quali.sup=1)
plot(pca_model,choix='ind',label='quali')

to_plot <- data_frame(speech=combined$speech,x=pca_model$ind$coord[,1],
                      y=pca_model$ind$coord[,2])

# to_plot_sum <- group_by(to_plot,speech) %>%  summarize(mean_val_x=weighted.mean(x, 1/(distance_geo + distance_wifi)),
#                                                 mean_val_y=weighted.mean(y, 1/(distance_geo + distance_wifi)))
# to_plot <- left_join(to_plot,to_plot_sum,by='speech')

filter(to_plot,x<50,y<50) %>% 
ggplot(aes(y=y,x=x,colour=speech)) + geom_point(alpha=0.1) + geom_text(aes(label=speech)) +
  theme_minimal()

# Try it without clustering on labels first



# Let's cluster this sucker

# cluster_pca <- HCPC(pca_model_trans,nb.clust=5)
# 
# coords <- data_frame(x=pca_model$ind$coord[,1],y=pca_model$ind$coord[,2],labels=for_dist_spread$location)
# full_plot <- ggplot(coords,aes(y=y,x=x,label=labels)) + geom_text(alpha=0.5)+ theme_minimal()
# # Test with a person prediction
# sample_locs <- sample(unique(sample_beacon$location),3)
# 
# sample_people <- filter(sample_beacon, location %in% sample_locs) %>% select(beacon_id,location,rssi,files_id)
# 
# for(i in sample_locs) {
#   
#   this_sample <- filter(sample_people,location==i) %>% distinct(beacon_id,.keep_all=TRUE) 
#   to_join <- data_frame(all_beacon=unique(for_dist$beacon_id))
#   this_sample <- left_join(to_join,this_sample,c('all_beacon'='beacon_id')) %>% mutate(rssi=coalesce(rssi,-999L)) %>% 
#     select(all_beacon,rssi) %>% spread(all_beacon,rssi) %>% mutate(location=i) %>% select(location,everything())
#   
#   outcome <- predict(pca_model,newdata = this_sample)
#   full_plot <- full_plot + geom_label(data=data_frame(x=outcome$coord[,1],y=outcome$coord[,2],labels=i),
#                                      aes(x=x,y=y,label=labels),colour='red')
#   # Need to add in -999 for all missing beacon data
#   
#   
#   
# }
# print(full_plot)
# 
# 
# ggsave('full_plot.png',width=16,height=10,units='in')
    # Let's try rpart

#rpart_model <- rpart(location~ rssi + new_beacon_id,data=sample_beacon)
  
  
# Let's try random Forests

 # forest_model <- randomForest(x=feature_vector,y=as.factor(code_vector$location))
 # print(forest_model)
