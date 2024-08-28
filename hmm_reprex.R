### HMM example - problem with negative variances ### 
### Reproducible example ###

# Load packages ####

library(momentuHMM)

# Load data ####

test_data <- readRDS("output/reprex_data.rds")

# Prepare data ####

hmm_data <- test_data %>%  
  mutate(ID = burst_) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Fit model ####

# 3 states, no covariates 

mod <- fitHMM(data = hmm_data,
              nbStates = 3,
              dist = list(step = "gamma", angle = "vm"),
              estAngleMean = list(angle = TRUE),
              Par0 = list(step = c(mean_1 = 100, 
                                   mean_2 = 500,
                                   mean_3 = 4000, 
                                   sd_1 = 50, 
                                   sd_2 = 500,
                                   sd_3 = 2000),
                                   # zeromass_1 = 0.5, 
                                   # zeromass_2 = 0.001, 
                                   # zeromass_3 = 0.001),
                          angle = c(mean_1 = pi,
                                    mean_2 = pi,
                                    mean_3 = 0,
                                    concentration_1 = 0.1,
                                    concentration_2 = 0.5,
                                    concentration_3 = 0.99)))

# Look at var-covar matrix ####

# Negative or zero values
mod$mod$Sigma

# Not only in covariances, but variances too
diag(mod$mod$Sigma)
