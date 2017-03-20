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

all_files <- list.files('0306csv/BeaconData/',pattern='.csv',full.names=TRUE)

all_files <- lapply(all_files,read_csv)

all_files <- lapply(1:length(all_files), function(x) mutate(all_files[[x]],files_id=paste0('beacon_',x)))

beacon_data <- bind_rows(all_files)

#sample_beacon <- sample_n(beacon_data,10000) %>% mutate(beacon_id = paste(major,minor,sep='_'))
sample_beacon <- mutate(beacon_data,beacon_id = paste(major,minor,sep='_'))

#Make the beacon IDs easier to use

new_ids <- data_frame(beacon_id=unique(sample_beacon$beacon_id),
                      new_beacon_id=paste0('ID_',1:length(unique(sample_beacon$beacon_id))))
sample_beacon %<>% left_join(new_ids) %>% mutate(new_beacon_id=as.factor(new_beacon_id))


#sample_beacon <- mutate(beacon_data,beacon_id = as.factor(paste(major,minor,sep='_')))

# group_by(sample_beacon,location,beacon_id) %>% filter(!is.na(location)) %>% summarize(avg_rssi=mean(rssi)) %>%
#   ggplot(aes(x=beacon_id,y=avg_rssi)) + geom_point() + facet_wrap(~location) +
#   theme(axis.text.x=element_blank())

x <- seq(0.1, 5, by = 0.05)
y <- log(x) + rnorm(x, sd = 0.2)

m <- maxent(x, y)
new <- predict(m, x)

code_vector <- model.frame(location~ rssi + beacon_id + rssi*beacon_id,data=sample_beacon)
feature_vector <- model.matrix(location~ rssi + beacon_id + rssi*beacon_id,data=code_vector)

first_model <- maxent(feature_matrix = feature_vector,code_vector = as.factor(code_vector$location))

prediction <- predict(first_model,feature_vector)

predict_vals <- prediction[,1]

lookat <- data_frame(predict_vals,orig=code_vector$location)

sum(lookat$predict_vals==lookat$orig)/nrow(lookat)

results <- group_by(lookat,orig) %>% summarize(accuracy=sum(predict_vals==orig)/length(orig)) %>% arrange(accuracy)

sample_beacon %>% filter(location %in% c('1641','2625','2541','5587')) %>% group_by(location,beacon_id) %>%  summarize(avg_rssi=mean(rssi)) %>%
  ggplot(aes(x=beacon_id,y=avg_rssi)) + geom_point() + facet_wrap(~location) + theme_minimal() +
  theme(axis.text.x=element_blank())


# Try MDS
# Create distance matrix
# As with IRT, missing values are allowed

for_dist <- select(sample_beacon,beacon_id,location,rssi) %>% filter(!is.na(location)) %>% group_by(beacon_id,location) %>% 
  summarize(avg_rssi=mean(rssi,na.rm=TRUE)) %>% filter(!is.na(avg_rssi)) %>% spread(key=beacon_id,value=avg_rssi)
to_dist <- lapply(for_dist,function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame() %>% dist()
mds_model <- cmdscale(to_dist)

x <- mds_model[,1]
y <- mds_model[,2]

data_frame(x=x,y=y,labels=for_dist$location) %>% ggplot(aes(x=x,y=y)) + geom_text(aes(label=labels),check_overlap = TRUE) +
  theme_minimal()

#Now PCA
for_dist <- select(sample_beacon,beacon_id,location,rssi) %>% filter(!is.na(location)) %>% group_by(beacon_id,location) %>% 
  summarize(avg_rssi=mean(rssi,na.rm=TRUE))
for_dist_spread <- spread(for_dist,key=beacon_id,value=avg_rssi)
to_dist <- lapply(for_dist_spread,function(x) {
  x[is.na(x)] <- -999
  x
}) %>% as_data_frame() %>% as.matrix()

pca_model <- PCA(to_dist,quali.sup=1)
plot(pca_model,choix='ind',label='quali')
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
