# # # # # # # # # # # # # # # # # # # # # # # # # 
# Prep data for multiple imputation HMM # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Nov 16, 2020 # # # # # # # # # # # # #
# Last updated Nov 17, 2020 # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# The approach I used in hyenas_ctmm.R to regularize the hyena data to 1h 
# resolution results in linearly interpolated tracks instead of a random 
# realization of locations based on the underlying movement model. Here,
# I use multiple imputation to draw several realization of missing data to 
# fill gaps and then fit an HMM to each of those. 

# Load packages ####

library(tidyverse)
library(momentuHMM)
library(lubridate)

# Set seed ####

set.seed(1)

# Load cleaned data ####

hyenas <- readRDS("output/hyenas_data_cleaned.rds")

# Format data ####

# Round timestamps and add ID
hyenas <- hyenas %>% 
  filter(minute(timestamp) %in% c(0:3, 57:59)) %>% 
  mutate(timestamp_original = timestamp) %>% 
  mutate(timestamp = round(timestamp, units = "hours")) %>% 
  mutate(ID = id) %>% 
  mutate(timestamp = ymd_hms(timestamp),
         timestamp_original = ymd_hms(timestamp_original)) 

# Remove duplicates
hyenas <- hyenas[!duplicated(paste0(hyenas$timestamp, hyenas$ID)),]

# Fit CTMM ####

# If I associate covariates with the data before I do crawlWrap, the predicted
# locations will have NAs. 

ct <- crawlWrap(obsData = hyenas,
                timeStep = "1 hours",
                coord = c("x", "y"), 
                Time.name = "timestamp")

# Draw multiple imputations ####

mi <- MIfitHMM(ct, nSims = 10, fit = FALSE)

# Process multiple imputation data ####

# Do the following:
# 1. Filter night only
# 2. Change the ID to be the individual-night, not the individual
# 3. Toss incomplete individual-nights (< 15 data points)
# 4. Calculate time from midnight
source("calc_time_from_midnight.R")
# 5. Associate individual attributes

info <- readxl::read_xlsx("input/ID.xlsx") %>% 
  janitor::clean_names() %>% 
  # Combine sex and reproductive status into one
  mutate(sex_rs = case_when(
    sex == "Female" & comments == "Mother with pups" ~ "RF",
    sex == "Female" & is.na(comments) ~ "F", 
    sex == "Male" ~ "M"
  ), 
  legs = case_when(
    comments == "3 legs" ~ "3", 
    TRUE ~ "4"
  )) %>% 
  dplyr::select(id, district = location, sex_rs, legs)

# Do it:
mi$miData <- lapply(mi$miData, function(x) {
  
  x <- x %>% 
    filter(hour(timestamp) < 8 | hour(timestamp) > 16) %>% 
    group_by(ID) %>% 
    mutate(lag = as.numeric((timestamp) - dplyr::lag(timestamp))) %>% 
    mutate(lag = case_when(
      is.na(lag) ~ 1,
      TRUE ~ lag
    )) %>% 
    mutate(burst_ = NA_real_)
  
  # Add a burst column that designates the individual-night
  # Calculate when one night starts based on lag time between locations
  
  x$burst_[1] <- 1
  
  for (i in 2:nrow(x)) {
    if (x$ID[i-1] != x$ID[i] | x$lag[i] > 1) {
      x$burst_[i] <- x$burst_[i-1] + 1
    } else if (x$lag[i] == 1) {
      x$burst_[i] <- x$burst_[i-1] 
    }
  }
  
  toss <- x %>% 
    group_by(burst_) %>% 
    tally() %>% 
    filter(n < 15) %>% 
    pull(burst_)
  
  x <- x %>% 
    filter(!burst_ %in% toss) %>% 
    calc_tfm() %>% 
    mutate(tfm = as.numeric(time_from_midnight)) %>% 
    left_join(info, by = c("ID" = "id")) %>% 
    as.data.frame()
  
  class(x) <- c("momentuHMMData", "data.frame")
  
  return(x)
  
})

# Fit multiple imputation HMM ####

form <- ~ tfm

hmm_tfm <- MIfitHMM(miData = mi$miData,
                    nSims = 10,
                  nbStates = 3,
#                  formula = form,
                  dist = list(step = "gamma", angle = "vm"),
                  estAngleMean = list(angle = TRUE),
                  Par0 = list(step = c(mean_1 = 100, 
                                       mean_2 = 500,
                                       mean_3 = 4000, 
                                       sd_1 = 50, 
                                       sd_2 = 500,
                                       sd_3 = 2000,  
                                       zeromass_1 = 0.5, 
                                       zeromass_2 = 0.001, 
                                       zeromass_3 = 0.001),
                              angle = c(mean_1 = pi,
                                        mean_2 = pi,
                                        mean_3 = 0,
                                        concentration_1 = 0.1,
                                        concentration_2 = 0.5,
                                        concentration_3 = 0.99)))

saveRDS(hmm_tfm, "output/hmm_10mi_test_null.rds")

hmm_tfm <- MIfitHMM(miData = mi$miData,
                    nSims = 10,
                    nbStates = 2,
                    #                  formula = form,
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    Par0 = list(step = c(mean_1 = 200, 
                                         mean_2 = 3000,
                                         sd_1 = 400, 
                                         sd_2 = 2000,
                                         zeromass_1 = 0.5, 
                                         zeromass_2 = 0.001),
                                angle = c(mean_1 = pi,
                                          mean_2 = 0,
                                          concentration_1 = 0.1,
                                          concentration_2 = 0.9)))

saveRDS(hmm_tfm, "output/hmm_10mi_test_null_2states.rds")
