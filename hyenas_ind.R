# # # # # # # # # # # # # # # # # # # # # # # # # 
# Individual datasets # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Oct 27, 2020 # # # # # # # # # # # # #
# Last updated Oct 27, 2020 # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # #

# In our last meeting on October 15, we talked about possible ways to control
# for individual identity in our hyena models. One option was to fit models to
# individual datasets separately. To check if this is possible, I need to verify
# that the model can converge with a limited amount of data. I'll find the 
# individual that has the smallest dataset as a test.

# Another option is to balance the design by individual, sex, and possibly also
# age, and month of the year (if there are any seasonal effects)

# Load packages ####

library(tidyverse)
library(momentuHMM)

# Load CTMM nightly data ####

hcd_night <- readRDS("output/hcd_1hCTMM_night.rds")

# Find individual with the smallest dataset ####

test_ind <- hcd_night %>% 
  group_by(id) %>% 
  tally() %>% 
  arrange(n) %>% 
  slice(1) %>% 
  pull(id)

# Try to fit the model on him/her ####

source("calc_time_from_midnight.R")

hmm_data <- hcd_night %>% 
  filter(id == test_ind) %>% 
  calc_tfm() %>%
  mutate(ID = burst_,
         tfm = as.numeric(time_from_midnight)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

form <- ~ tfm

hmm_tfm <- fitHMM(data = hmm_data,
                  nbStates = 3,
                  stationary = TRUE,
                  formula = form,
                  dist = list(step = "gamma", angle = "vm"),
                  estAngleMean = list(angle = TRUE),
                  Par0 = list(step = c(mean_1 = 100, 
                                       mean_2 = 500,
                                       mean_3 = 4000, 
                                       sd_1 = 50, 
                                       sd_2 = 500,
                                       sd_3 = 2000),
                              angle = c(mean_1 = pi,
                                        mean_2 = pi,
                                        mean_3 = 0,
                                        concentration_1 = 0.1,
                                        concentration_2 = 0.5,
                                        concentration_3 = 0.99)))

plot(hmm_tfm)

# ind 1, 3, 4 do not fit state 2
# ind 2, 5 do not converge
# ind 6 stalls for 30 min
# ind 7 and 8 work
# some of them need zeromass, some don't 