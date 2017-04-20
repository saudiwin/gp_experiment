require(dplyr)
require(tidyr)
require(rstan)
require(bayesplot)
require(abind)
require(trelliscopejs)

# model parameters --------------------------------------------------------

niter <- 2000
nwarmup <- 1000
num_points <- 50
test_fix <- 3

true_points <- matrix(nrow=num_points,c(rnorm(num_points,0,10),rnorm(num_points,0,10)))

distance_true <- proxy::dist(true_points) %>% as.matrix() 

# Add lognormal errors

distance_obs <- apply(distance_true,c(1,2),function(x) {
  x <- x + rlnorm(1)
}) %>% as_data_frame()



# put data in long form so that each point is matched with each other point exactly once
# total of n(n-1)/2 pairs

long_dist <- gather(distance_obs,edge1,distance) %>% mutate(edge2=rep(1:num_points,num_points)) %>% 
  group_by(edge1) %>% 
  mutate(unnecessary=if_else((edge1>edge2) | (edge1==edge2),T,F)) %>% 
  filter(!unnecessary) %>% ungroup

compiled <- stan_model(model_name = 'MDSBayes',file = 'stan_mds.stan')

# Check to see how many fixed points we need for identification

all_sims <- lapply(1:test_fix,function(i) { 
fixed_points <- matrix(true_points[1:i,],ncol=2)

stan_data <- list(J=num_points,
                  num_dist=nrow(long_dist),
                  points_fix_num=i,
                  points_id=select(long_dist,edge1,edge2),
                  points_dist=long_dist$distance,
                  points_fix_vals=fixed_points)
output <- sampling(compiled,chains=4,cores=4,iter=niter,warmup=nwarmup,data=stan_data)
})

  

to_plot_params <- lapply(all_sims,function(x) {
  all_out <- as.array(x)
  all_pars <- attributes(all_out)$dimnames$parameters
  all_out <- all_out[,,grepl('unstd\\[',all_pars)]  
  return(all_out)
})

to_plot_data <- lapply(to_plot_params,function(x) {
  outplot <- mcmc_recover_intervals(x,true=c(true_points))
  return(outplot$data)
})

combine_data <- bind_rows(to_plot_data)

combine_data$Batch <- factor(paste0(rep(1:test_fix,each=(num_points*2)),'_Fixed'))
combine_data <- gather(combine_data,type,estimate,True,Point)

ggplot(combine_data,aes(y=estimate,x=Parameter,colour=type)) + geom_point(aes(shape=type),size=1) + 
  geom_linerange(aes(ymin=l,ymax=u)) +
  theme_minimal() +
  facet_trelliscope(~Batch) +
  theme(panel.grid = element_blank()) +
  scale_color_brewer()
