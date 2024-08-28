# Try to keep NAs instead of interpolating missing data in the tracks

library(momentuHMM)
library(tidyverse)
library(lubridate)

set.seed(7)

hyenas <- readRDS("output/hyenas_data_cleaned.rds")

# Round timestamps and add ID
hyenas <- hyenas %>% 
  filter(minute(timestamp) %in% c(0:3, 57:59)) %>% 
  mutate(timestamp_original = timestamp) %>% 
  mutate(timestamp = round(timestamp, units = "hours")) %>% 
  mutate(ID = id) %>% 
  mutate(timestamp = ymd_hms(timestamp, tz = "Israel"),
         timestamp_original = ymd_hms(timestamp_original, tz = "Israel")) 

# Remove duplicates
hyenas <- hyenas[!duplicated(paste0(hyenas$timestamp, hyenas$ID)),]

hyenas_night <- hyenas %>% 
  # filter nighttime data only (6PM to 8AM)
  filter(hour(timestamp) > 17 | hour(timestamp) < 9) %>% 
  group_by(ID) %>% 
  mutate(lag = as.numeric((timestamp) - lag(timestamp))) %>% 
  mutate(lag = case_when(
    is.na(lag) ~ 1,
    TRUE ~ lag
  )) %>% 
  mutate(burst_ = NA_real_)

hyenas_night$burst_[1] <- 1

for (i in 2:nrow(hyenas_night)) {
  if (hyenas_night$ID[i-1] != hyenas_night$ID[i] | hyenas_night$lag[i] > 1) {
    hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] + 1
  } else if (hyenas_night$lag[i] == 1) {
    hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] 
  }
}

# Get rid of bursts that are too short to fit a model
drop <- hyenas_night %>% 
  group_by(burst_, ID) %>% 
  tally() %>% 
  filter(n < 7) 
# 15 locations expected per night, we'll take nights that have at least 7

hyenas_night <- hyenas_night %>% 
  filter(!burst_ %in% drop$burst_)

# burst_ becomes ID
hyenas_night <- hyenas_night %>%
  ungroup() %>% 
  dplyr::select(id, ID = burst_, timestamp, x, y) 

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
  dplyr::select(id, district = location, sex, sex_rs, legs)
# I left age out because the age is NA for the hyenas in the South District 
# and all the others are adult except for one young adult. 

# In the data, the IDs of the South District hyenas follow the pattern 
# sex-number, in the info file they follow number-sex. Change these manually:
info$id[9:13] <- c("Female16", "Female32", "Female40", "Male32", "Male33")

hyenas_night <- hyenas_night %>% 
  left_join(info, by = "id")

# Add NAs where a location was expected but not taken ####

ranges <- hyenas_night %>% 
  group_by(ID) %>% 
  summarize(start = ymd_hms(paste0(min(as_date(timestamp)), " 18:00:00"), tz = "Israel"),
            end = ymd_hms(paste0(min(as_date(timestamp)) + days(1), " 08:00:00"), tz = "Israel"))

hyenas_night_nas <- data.frame(id = NA,
                               ID = NA, 
                               timestamp = ymd_hms(NA, tz = "Israel"),
                               x = NA,
                               y = NA, 
                               district = NA,
                               sex = NA,
                               sex_rs = NA, 
                               legs = NA)

for (i in 1:length(unique(hyenas_night$ID))) {
  
  who <- unique(hyenas_night$ID)[i]
  
  ts <- seq(ranges[i,]$start, ranges[i, ]$end, by = 3600)
  id <- unique(hyenas_night[hyenas_night$ID == who, ]$id)
  
  coord <- hyenas_night %>% 
    filter(ID == who) %>% 
    dplyr::select(id, ID, timestamp, x, y)
  
  info <- hyenas_night %>% 
    filter(ID == who) %>% 
    dplyr::select(id, ID, district, sex, sex_rs, legs) %>% 
    distinct()
  
  res <- data.frame(id = id,
                    ID = who,
                    timestamp = ts) %>% 
    left_join(coord, by = c("id", "ID", "timestamp")) %>% 
    left_join(info, by = c("id", "ID"))
  
  hyenas_night_nas <- rbind(hyenas_night_nas, res)
  
}

hyenas_night_nas <- hyenas_night_nas[-1, ]

# Calculate time from midnight ####

source("calc_time_from_midnight.R")

hyenas_night_tfm <- hyenas_night_nas %>% 
  mutate(x_original = x,
         y_original = y) %>% 
  # when coordinates are NA, we can't compute tfm because sunrise and 
  # sunset depend on the location. So we can use the same x and y as the 
  # previous location that wasn't NA.
  mutate(x = case_when(
    is.na(x_original) ~ zoo::na.locf(x_original, na.rm = FALSE),
    TRUE ~ x_original
  ),
  y = case_when(
    is.na(y_original) ~ zoo::na.locf(y_original, na.rm = FALSE),
    TRUE ~ y_original
  )) %>% 
  # There are still two NAs at the beginning of the data frame, where there's no
  # prior available location to use. Assign these manually:
  mutate(x = case_when(
    is.na(x) ~ 693385.4,
    TRUE ~ x
  ),
  y = case_when(
    is.na(y) ~ 3525480,
    TRUE ~ y
  )) %>% 
  calc_tfm() %>% 
  mutate(tfm = as.numeric(time_from_midnight)) %>% 
  dplyr::select(-time_from_midnight) %>% 
  dplyr::select(-x, -y) %>% 
  rename(x = x_original, y = y_original)

hmm_data <- hyenas_night_tfm %>% 
  mutate(sex = factor(sex), 
         sex_rs = factor(sex_rs),
         district = factor(district),
         legs = factor(legs)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Write a function to compute shape and scale of Gamma given mean and sd
# Otherwise I'm just picking initial values with my eyes closed
gamma_pars <- function(mean = NULL, sd = NULL, 
                       scale = NULL, shape = NULL) {
  if(exists("mean") & exists("sd")) {
    
    scale_theta <- sqrt(sd^2/mean)
    shape_kappa <- mean/scale_theta
  
    pars <- c(scale_theta = scale_theta,
              shape_kappa = shape_kappa)
  } else if(exists("scale") & exists("shape")) {
    
    mean <- shape * scale
    sd <- sqrt(scale^2 * shape)
    
    pars <- c(mean = mean,
              sd = sd)
    
  }
  
  return(pars)
}

test1 <- gamma_pars(mean = 50, sd = 50)
test2 <- gamma_pars(mean = 1000, sd = 2500)
test3 <- gamma_pars(mean = 3000, sd = 7000)

hist(hmm_data$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(1:12000), shape = test1[2], scale = test1[1]), col = "red")
lines(dgamma(x = c(200:12000), shape = test2[2], scale = test2[1]), col = "blue")
lines(dgamma(x = c(1:12000), shape = test3[2], scale = test3[1]), col = "green")

hmm_null3 <- fitHMM(data = hmm_data,
                   nbStates = 3,
                   dist = list(step = "gamma", angle = "vm"),
                   estAngleMean = list(angle = TRUE),
                   Par0 = list(step = c(mean_1 = 50, 
                                        mean_2 = 1000,
                                        mean_3 = 3500, 
                                        sd_1 = 50, 
                                        sd_2 = 2500,
                                        sd_3 = 8000,  
                                        zeromass_1 = 0.5, 
                                        zeromass_2 = 0.001, 
                                        zeromass_3 = 0.001),
                               angle = c(mean_1 = pi,
                                         mean_2 = pi,
                                         mean_3 = 0,
                                         concentration_1 = 0.1,
                                         concentration_2 = 0.5,
                                         concentration_3 = 0.99)))

diag(hmm_null3$mod$Sigma)
plot(hmm_null3)

results_hmm3_null <- hmm_data %>% 
  mutate(state = viterbi(hmm_null3)) %>% 
  dplyr::select(id, timestamp, x, y, state)

write.csv(results_hmm3_null, "output/results_hmm_3states_null.csv",
          row.names = FALSE)
  
form <- ~ sex_rs + district + tfm

hmm_full3_rs <- fitHMM(data = hmm_data,
                    nbStates = 3,
                    formula = form, 
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    Par0 = list(step = c(mean_1 = 50, 
                                         mean_2 = 1000,
                                         mean_3 = 3000, 
                                         sd_1 = 50, 
                                         sd_2 = 2500,
                                         sd_3 = 7000,  
                                         zeromass_1 = 0.5, 
                                         zeromass_2 = 0.001, 
                                         zeromass_3 = 0.001),
                                angle = c(mean_1 = pi,
                                          mean_2 = pi,
                                          mean_3 = 0,
                                          concentration_1 = 0.1,
                                          concentration_2 = 0.5,
                                          concentration_3 = 0.99)))

diag(hmm_full3$mod$Sigma)
diag(hmm_full3_rs$mod$Sigma)
plot(hmm_full3)

results_hmm3_full_rs <- hmm_data %>% 
  mutate(state = viterbi(hmm_full3_rs)) %>% 
  dplyr::select(id, timestamp, x, y, state)

write.csv(results_hmm3_full_rs, "output/results_hmm_3states_full.csv",
          row.names = FALSE)

test1 <- gamma_pars(mean = 100, sd = 300)
test2 <- gamma_pars(mean = 1500, sd = 5000)

hist(hmm_data$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(1:12000), shape = test1[2], scale = test1[1]), col = "red")
lines(dgamma(x = c(200:12000), shape = test2[2], scale = test2[1]), col = "blue")

hmm_full2 <- fitHMM(data = hmm_data,
                    nbStates = 2,
                    formula = form, 
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    Par0 = list(step = c(mean_1 = 100, 
                                         mean_2 = 1500,
                                         sd_1 = 300, 
                                         sd_2 = 5000,
                                         zeromass_1 = 0.5, 
                                         zeromass_2 = 0.001),
                                angle = c(mean_1 = pi,
                                          mean_2 = 0,
                                          concentration_1 = 0.1,
                                          concentration_2 = 0.9)))

diag(hmm_full2$mod$Sigma)
plot(hmm_full2)
