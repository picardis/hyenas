### HMM example - problem with negative variances ### 

# Load packages ####

library(tidyverse)
library(momentuHMM)

# Load data ####

hcd_night <- readRDS("output/hcd_1hCTMM_night.rds")

# Find individual with the smallest dataset ####

test_ind <- hcd_night %>% 
  group_by(id) %>% 
  tally() %>% 
  arrange(n) %>% 
  slice(3) %>% 
  pull(id)

test_data <- hcd_night %>% 
  filter(id == test_ind) %>% 
  as.data.frame()

saveRDS(test_data, "output/reprex_data.rds")
