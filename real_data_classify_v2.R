require(mnlogit)
require(maxent)
require(readr)
require(dplyr)
require(ggplot2)
require(randomForest)
require(rpart)
require(magrittr)
require(mirt)
require(tidyr)
require(FactoMineR)

# Load beacon/location data and merge

wifi_data <- read_csv('client_data/wifi.csv')

# Remove BSSID's with very few values

wifi_data <- group_by(wifi_data,BSSID) %>% filter(n()>3)

# Try MDS
# Create distance matrix
# As with IRT, missing values are allowed

for_dist <- select(wifi_data,BSSID,level,timestamp) %>% filter(!is.na(level)) %>% group_by(BSSID,timestamp) %>% 
  summarize(avg_level=mean(level,na.rm=TRUE)) %>% filter(!is.na(avg_level)) %>% spread(key=BSSID,value=avg_level)
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
for_dist <- select(wifi_data,BSSID,level,timestamp) %>% filter(!is.na(level)) %>% group_by(BSSID,timestamp) %>% 
  summarize(avg_level=mean(level,na.rm=TRUE)) %>% filter(!is.na(avg_level)) %>% spread(key=BSSID,value=avg_level)
to_dist <- lapply(select(for_dist,-timestamp),function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame()

to_dist_trans <- t(as.matrix(to_dist))
to_dist$timestamp <- for_dist$timestamp
to_dist <- select(to_dist,timestamp,everything())

pca_model_trans <- PCA(to_dist_trans)
pca_model <- PCA(to_dist,quali.sup=1)
plot(pca_model,choix='ind',label='quali')

# Let's cluster this sucker

cluster_pca <- HCPC(pca_model_trans,nb.clust=5)

coords <- data_frame(x=pca_model$ind$coord[,1],y=pca_model$ind$coord[,2],labels=for_dist_spread$location)
full_plot <- ggplot(coords,aes(y=y,x=x,label=labels)) + geom_text(alpha=0.5)+ theme_minimal()
# Test with a person prediction
sample_locs <- sample(unique(sample_beacon$location),3)

sample_people <- filter(sample_beacon, location %in% sample_locs) %>% select(beacon_id,location,rssi,files_id)

for(i in sample_locs) {
  
  this_sample <- filter(sample_people,location==i) %>% distinct(beacon_id,.keep_all=TRUE) 
  to_join <- data_frame(all_beacon=unique(for_dist$beacon_id))
  this_sample <- left_join(to_join,this_sample,c('all_beacon'='beacon_id')) %>% mutate(rssi=coalesce(rssi,-999L)) %>% 
    select(all_beacon,rssi) %>% spread(all_beacon,rssi) %>% mutate(location=i) %>% select(location,everything())
  
  outcome <- predict(pca_model,newdata = this_sample)
  full_plot <- full_plot + geom_label(data=data_frame(x=outcome$coord[,1],y=outcome$coord[,2],labels=i),
                                     aes(x=x,y=y,label=labels),colour='red')
  # Need to add in -999 for all missing beacon data
  
  
  
}
print(full_plot)


ggsave('full_plot.png',width=16,height=10,units='in')
    # Let's try rpart

#rpart_model <- rpart(location~ rssi + new_beacon_id,data=sample_beacon)
  
  
# Let's try random Forests

 # forest_model <- randomForest(x=feature_vector,y=as.factor(code_vector$location))
 # print(forest_model)
