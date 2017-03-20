# Generate test data that is highly skewed

require(dplyr)
require(mclust)

latitude <- c(rnorm(1000,-10,1),
              rnorm(8000,0,1),
              rnorm(1000,10,1))

longitude <- c(rnorm(1000,5,1),
               rnorm(8000,20,1),
               rnorm(1000,50,1))

test_data <- data_frame(latitude,longitude)


#This simple random sample doesn't return many records from locations with small density

simple_sample <- sample_n(test_data,2000)

plot(density(simple_sample$latitude))
plot(density(simple_sample$longitude))

# We can fix this by first clustering the data around longitude/latitude and identifying the distinct clusters
# where the individual traveled to

mclust_result <- densityMclust(test_data)
groups <- data_frame(clusters=mclust_result$classification) %>% left_join(data_frame(clusters=1:mclust_result$G,
                                                                                     proportion=mclust_result$parameters$pro),
                                                                          by='clusters') %>% 
  mutate(to_sample=(1/proportion))

# Sample the join_data so that it balances the geographical clusters

balanced_sample <- sample_n(test_data,2000,weight = groups$to_sample)

# Sample is now more balanced

plot(density(balanced_sample$latitude))
plot(density(balanced_sample$longitude))

