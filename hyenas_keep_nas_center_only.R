# # # # # # # # # # # # # # # # # # # # # # # # # 
# Behavioral segmentation of striped hyena tracks
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Feb 17, 2021 # # # # # # # # # # # # #
# Last updated Feb 17, 2021 # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# After talking with Einat, we decided to analyze data from the two different
# districts separately. We also decided to try to analyze the 20-min resolution
# data as well as the 1-h. 

library(momentuHMM)
library(tidyverse)
library(lubridate)
library(CircStats)

set.seed(7)

# 1-hour data from Central district ####

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
  # 6168 is a male, not a female!
  mutate(sex = case_when(
    id == "6168" ~ "Male",
    TRUE ~ sex
  )) %>% 
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

# Check if this introduced any bursts with < 7 locations
keepers <- hyenas_night_nas %>% 
  mutate(location_exists = case_when(
    is.na(x) | is.na(y) ~ 0,
    TRUE ~ 1
  )) %>% 
  group_by(ID) %>% 
  summarize(n = sum(location_exists)) %>% 
  arrange(n) %>% 
  filter(n > 7) %>% 
  pull(ID)

hyenas_night_nas <- hyenas_night_nas %>% 
  filter(ID %in% keepers)

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

# Keep only hyenas from Center district
hyenas_night_tfm <- hyenas_night_tfm %>% 
  filter(district == "Center") %>% 
  dplyr::select(-district) %>% 
  # and only those with 4 legs
  filter(legs == 4)

hmm_data <- hyenas_night_tfm %>% 
  mutate(sex = factor(sex), 
         sex_rs = factor(sex_rs)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Summaries ####

# Distance travelled (per night)
hmm_data %>% 
  group_by(ID, sex) %>% 
  summarize(dist_trav = sum(step, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  summarize(min_dist_trav = min(dist_trav),
            mean_dist_trav = mean(dist_trav),
            max_dist_trav = max(dist_trav),
            sd_dist_trav = sd(dist_trav))
  
# Speed (per hour)
hmm_data %>% 
  group_by(sex) %>% 
  summarize(min_dist_hour = min(step, na.rm = T),
            mean_dist_hour = mean(step, na.rm = T),
            max_dist_hour = max(step, na.rm = T),
            sd_dist_hour = sd(step, na.rm = T))

# Net squared displacement
traj <- adehabitatLT::as.ltraj(hmm_data[!is.na(hmm_data$x), c("x", "y")],
                               date = hmm_data[!is.na(hmm_data$x),]$timestamp,
                               id = hmm_data[!is.na(hmm_data$x),]$id,
                               burst = hmm_data[!is.na(hmm_data$x),]$ID)
hmm_data %>% 
  left_join(adehabitatLT::ld(traj)) %>% 
  group_by(ID, sex) %>% 
  summarize(nsd = max(R2n, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  summarize(min_net_displ = sqrt(min(nsd, na.rm = T)),
            mean_net_displ = sqrt(mean(nsd, na.rm = T)),
            max_net_displ = sqrt(max(nsd, na.rm = T)),
            sd_net_displ = sqrt(sd(nsd, na.rm = T)))

# T-test for male vs. female step length
stepf <- hmm_data[hmm_data$sex == "Female", ]$step
stepm <- hmm_data[hmm_data$sex == "Male", ]$step
t.test(stepf, stepm, alternative = "less")
t.test(stepf, stepm, alternative = "two.sided")

# T-test for male vs. female nightly distance
nightdf <- hmm_data %>% 
  group_by(ID, sex) %>% 
  summarize(dist_trav = sum(step, na.rm = T)) %>% 
  pivot_wider(names_from = sex, values_from = dist_trav) %>% 
  pull(Female) %>% 
  na.omit()
nightdm <- hmm_data %>% 
  group_by(ID, sex) %>% 
  summarize(dist_trav = sum(step, na.rm = T)) %>% 
  pivot_wider(names_from = sex, values_from = dist_trav) %>% 
  pull(Male) %>% 
  na.omit()
t.test(nightdf, nightdm, alternative = "less")

# T-test for male vs. female NSD
nsdf <- hmm_data %>% 
  left_join(adehabitatLT::ld(traj)) %>% 
  group_by(ID, sex) %>%
  pivot_wider(names_from = sex, values_from = R2n) %>% 
  pull(Female) %>% 
  na.omit()
nsdm <- hmm_data %>% 
  left_join(adehabitatLT::ld(traj)) %>% 
  group_by(ID, sex) %>%
  pivot_wider(names_from = sex, values_from = R2n) %>% 
  pull(Male) %>% 
  na.omit()
t.test(nsdf, nsdm, alternative = "less")
  
# Make plot for movement summaries
mov_summaries <- hmm_data %>% 
  left_join(adehabitatLT::ld(traj)) %>% 
  group_by(ID, sex) %>% 
  summarize(dist_trav = sum(step, na.rm = T),
            net_disp = sqrt(max(R2n, na.rm = T)))

dpn <- ggplot(mov_summaries, aes(x = sex, y = dist_trav/1000, fill = sex)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = " ", y = "Distance travelled per night (km)") +
  theme(legend.position = "none")

nsd <- ggplot(mov_summaries, aes(x = sex, y = net_disp/1000, fill = sex)) +
  geom_boxplot() +
  ylim(0, 35) +
  theme_bw() +
  labs(x = " ", y = "Net displacement per night (km)") +
  theme(legend.position = "none")

dph <- ggplot(hmm_data, aes(x = sex, y = step/1000, fill = sex)) +
  geom_boxplot() +
  ylim(0, 35) +
  theme_bw() +
  labs(x = " ", y = "Distance travelled per hour (km)") +
  theme(legend.position = "none")

ggpubr::ggarrange(dph, dpn, nsd, nrow = 1)

ggsave("output/movement-summaries_2021-08-24.tiff",
       compression = "lzw", width = 10, height = 5, dpi = 300, scale = 0.7)

# Choose initial values

source("gamma_pars.R")

test1 <- gamma_pars(mean = 50, sd = 50)
test2 <- gamma_pars(mean = 1000, sd = 2500)
test3 <- gamma_pars(mean = 3000, sd = 5000)

hist(hmm_data$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(1:12000), shape = test1[2], scale = test1[1]), col = "red")
lines(dgamma(x = c(1:12000), shape = test2[2], scale = test2[1]), col = "blue")
lines(dgamma(x = c(1:12000), shape = test3[2], scale = test3[1]), col = "green")

# Null model
hmm_null3 <- fitHMM(data = hmm_data,
                    nbStates = 3,
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    Par0 = list(step = c(mean_1 = 50, 
                                         mean_2 = 1000,
                                         mean_3 = 3500, 
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

plot(hmm_null3)
# No negative values
diag(hmm_null3$mod$Sigma)
any(diag(hmm_null3$mod$Sigma) < 0)

results_hmm3_null <- hmm_data %>% 
  mutate(state = viterbi(hmm_null3)) %>% 
  dplyr::select(id, timestamp, x, y, state)

# Model with tfm only

form <- ~ tfm

hmm_tfm3 <- fitHMM(data = hmm_data,
                    nbStates = 3,
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    formula = form,
                    Par0 = list(step = c(mean_1 = 50, 
                                         mean_2 = 1000,
                                         mean_3 = 3000, 
                                         sd_1 = 50, 
                                         sd_2 = 2500,
                                         sd_3 = 5000,  
                                         zeromass_1 = 0.5, 
                                         zeromass_2 = 0.001, 
                                         zeromass_3 = 0.001),
                                angle = c(mean_1 = pi,
                                          mean_2 = pi,
                                          mean_3 = 0,
                                          concentration_1 = 0.1,
                                          concentration_2 = 0.5,
                                          concentration_3 = 0.99)))

saveRDS(hmm_tfm3, "output/hmm_tfm3.rds")

plot(hmm_tfm3)

diag(hmm_tfm3$mod$Sigma)
any(diag(hmm_tfm3$mod$Sigma) < 0)
sum(diag(hmm_tfm3$mod$Sigma) < 0)

results_hmm3_tfm <- hmm_data %>% 
  mutate(state = viterbi(hmm_tfm3)) %>% 
  dplyr::select(ind_night = ID, id, timestamp, x, y, state)

write.csv(results_hmm3_tfm, "output/hmm_3states_time-from-midnight_2021-07-19.csv")

# Full model (sex/rs + tfm)
form <- ~ sex_rs + tfm

hmm3_rs_tfm <- fitHMM(data = hmm_data,
                       nbStates = 3,
                       formula = form, 
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

plot(hmm3_rs_tfm)
# One negative value
diag(hmm3_rs_tfm$mod$Sigma)
any(diag(hmm3_rs_tfm$mod$Sigma) < 0)

results_hmm3_rs_tfm <- hmm_data %>% 
  mutate(state = viterbi(hmm3_rs_tfm)) %>% 
  dplyr::select(id, timestamp, x, y, state)

# 20-min data from Central district ####

hyenas <- readRDS("output/hyenas_data_cleaned.rds")

# Filter data at 20-min resolution
hyenas20 <- hyenas %>% 
  mutate(delta_t = timestamp - lag(timestamp)) %>% 
  filter(delta_t %in% c(18:22)) %>% 
  mutate(timestamp_original = timestamp) 

# Round timestamps
minute(hyenas20[minute(hyenas20$timestamp) > 18 & 
                  minute(hyenas20$timestamp) < 22, ]$timestamp) <- 20

minute(hyenas20[minute(hyenas20$timestamp) > 38 & 
                  minute(hyenas20$timestamp) < 42, ]$timestamp) <- 40

minute(hyenas20[minute(hyenas20$timestamp) > 58 | 
                  minute(hyenas20$timestamp) < 2, ]$timestamp) <- 0

# Check
table(minute(hyenas20$timestamp))

# Add ID
hyenas20 <- hyenas20 %>% 
  mutate(ID = id) %>% 
  mutate(timestamp = ymd_hms(timestamp, tz = "Israel"),
         timestamp_original = ymd_hms(timestamp_original, tz = "Israel")) 

# Remove duplicates
hyenas20 <- hyenas20[!duplicated(paste0(hyenas20$timestamp, hyenas20$ID)),]

hyenas20_night <- hyenas20 %>% 
  # filter nighttime data only (6PM to 8AM)
  filter(hour(timestamp) > 17 | hour(timestamp) < 9) %>% 
  mutate(delta_t = timestamp - lag(timestamp)) %>% 
  mutate(burst_ = NA_real_)

hyenas20_night$burst_[1] <- 1

for (i in 2:nrow(hyenas20_night)) {
  if (hyenas20_night$ID[i-1] != hyenas20_night$ID[i] | hyenas20_night$delta_t[i] > minutes(200)) {
    hyenas20_night$burst_[i] <- hyenas20_night$burst_[i-1] + 1
  } else if (hyenas20_night$delta_t[i] < 200) {
    hyenas20_night$burst_[i] <- hyenas20_night$burst_[i-1] 
  }
}

# Get rid of bursts that are too short to fit a model
drop <- hyenas20_night %>% 
  group_by(burst_, ID) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  filter(n < 10) 
# 43 locations expected per night, we'll take nights that have at least 10

hyenas20_night <- hyenas20_night %>% 
  filter(!burst_ %in% drop$burst_)

# burst_ becomes ID
hyenas20_night <- hyenas20_night %>%
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

hyenas20_night <- hyenas20_night %>% 
  left_join(info, by = "id")

# Add NAs where a location was expected but not taken ####

ranges <- hyenas20_night %>% 
  group_by(ID) %>% 
  summarize(start = ymd_hms(paste0(min(as_date(timestamp)), " 18:00:00"), tz = "Israel"),
            end = ymd_hms(paste0(min(as_date(timestamp)) + days(1), " 08:00:00"), tz = "Israel"))

hyenas20_night_nas <- data.frame(id = NA,
                               ID = NA, 
                               timestamp = ymd_hms(NA, tz = "Israel"),
                               x = NA,
                               y = NA, 
                               district = NA,
                               sex = NA,
                               sex_rs = NA, 
                               legs = NA)

for (i in 1:length(unique(hyenas20_night$ID))) {
  
  who <- unique(hyenas20_night$ID)[i]
  
  ts <- seq(ranges[i,]$start, ranges[i, ]$end, by = 1200)
  id <- unique(hyenas20_night[hyenas20_night$ID == who, ]$id)
  
  coord <- hyenas20_night %>% 
    filter(ID == who) %>% 
    dplyr::select(id, ID, timestamp, x, y)
  
  info <- hyenas20_night %>% 
    filter(ID == who) %>% 
    dplyr::select(id, ID, district, sex, sex_rs, legs) %>% 
    distinct()
  
  res <- data.frame(id = id,
                    ID = who,
                    timestamp = ts) %>% 
    left_join(coord, by = c("id", "ID", "timestamp")) %>% 
    left_join(info, by = c("id", "ID"))
  
  hyenas20_night_nas <- rbind(hyenas20_night_nas, res)
  
}

hyenas20_night_nas <- hyenas20_night_nas[-1, ]

# Drop bursts with no locations at all 
keepers <- hyenas20_night_nas %>% 
  mutate(location_exists = case_when(
    is.na(x) | is.na(y) ~ 0,
    TRUE ~ 1
  )) %>% 
  group_by(ID) %>% 
  summarize(n = sum(location_exists)) %>% 
  arrange(n) %>% 
  filter(n > 9) %>% 
  pull(ID)

hyenas20_night_nas <- hyenas20_night_nas %>% 
  filter(ID %in% keepers)

# Calculate time from midnight ####

source("calc_time_from_midnight.R")

hyenas20_night_tfm <- hyenas20_night_nas %>% 
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

# Keep only hyenas20 from Center district
hyenas20_night_tfm <- hyenas20_night_tfm %>% 
  filter(district == "Center") %>% 
  dplyr::select(-district)

hmm_data_20 <- hyenas20_night_tfm %>% 
  mutate(sex = factor(sex), 
         sex_rs = factor(sex_rs),
         legs = factor(legs)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Choose initial values

test1 <- gamma_pars(mean = 20, sd = 30)
test2 <- gamma_pars(mean = 400, sd = 200)
test3 <- gamma_pars(mean = 1000, sd = 800)

hist(hmm_data_20$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(1:12000), shape = test1[2], scale = test1[1]), col = "red")
lines(dgamma(x = c(1:12000), shape = test2[2], scale = test2[1]), col = "blue")
lines(dgamma(x = c(1:12000), shape = test3[2], scale = test3[1]), col = "green")

# Fit null HMM

hmm_null3_20 <- fitHMM(data = hmm_data_20,
                    nbStates = 3,
                    dist = list(step = "gamma", angle = "vm"),
                    estAngleMean = list(angle = TRUE),
                    Par0 = list(step = c(mean_1 = 20, 
                                         mean_2 = 400,
                                         mean_3 = 1000, 
                                         sd_1 = 30, 
                                         sd_2 = 200,
                                         sd_3 = 800,  
                                         zeromass_1 = 0.5, 
                                         zeromass_2 = 0.001, 
                                         zeromass_3 = 0.001),
                                angle = c(mean_1 = pi,
                                          mean_2 = pi,
                                          mean_3 = 0,
                                          concentration_1 = 0.1,
                                          concentration_2 = 0.5,
                                          concentration_3 = 0.99)))

plot(hmm_null3_20)
# Negative values?
diag(hmm_null3_20$mod$Sigma)
any(diag(hmm_null3_20$mod$Sigma) < 0)

# Fit full HMM

form <- ~ sex_rs + tfm

hmm_full3_20 <- fitHMM(data = hmm_data_20,
                       nbStates = 3,
                       formula = form, 
                       dist = list(step = "gamma", angle = "vm"),
                       estAngleMean = list(angle = TRUE),
                       Par0 = list(step = c(mean_1 = 20, 
                                            mean_2 = 400,
                                            mean_3 = 1000, 
                                            sd_1 = 30, 
                                            sd_2 = 200,
                                            sd_3 = 800,  
                                            zeromass_1 = 0.5, 
                                            zeromass_2 = 0.001, 
                                            zeromass_3 = 0.001),
                                   angle = c(mean_1 = pi,
                                             mean_2 = pi,
                                             mean_3 = 0,
                                             concentration_1 = 0.1,
                                             concentration_2 = 0.5,
                                             concentration_3 = 0.99)))

plot(hmm_full3_20)
# Negative values?
diag(hmm_full3_20$mod$Sigma)
any(diag(hmm_full3_20$mod$Sigma) < 0)

# Model predictions - transition probabilities ####

# Use model with tfm

source("predict_hmm.R")

pred_tfm <- data.frame(tfm = seq(-6, 8, length.out = 100))

preds <- predict_hmm(mod = hmm_tfm3, newdata = pred_tfm) %>% 
  mutate(state_from = case_when(
    word(tp, 1, 1, "t") == 1 ~ "From Resting",
    word(tp, 1, 1, "t") == 2 ~ "From Traveling",
    word(tp, 1, 1, "t") == 3 ~ "From Searching"
  ),
  state_to = case_when(
    word(tp, 2, 2, "t") == 1 ~ "To Resting",
    word(tp, 2, 2, "t") == 2 ~ "To Traveling",
    word(tp, 2, 2, "t") == 3 ~ "To Searching"
  ))

ggplot(preds, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ") +
  theme(legend.position = "bottom")

ggsave("output/hmm_tfm_tp_2021-07-19.tiff", compression= "lzw", 
       scale = 0.6, height = 5, width = 10, dpi = 400)

# Model predictions - initial state probabilities ####

isp <- data.frame(state = c("Resting", "Traveling", "Searching"),
                  est = hmm_tfm3$CIreal$delta$est[1, ],
                  lwr = hmm_tfm3$CIreal$delta$lower[1, ],
                  upr = hmm_tfm3$CIreal$delta$upper[1, ]) 

ggplot(isp, aes(x = state, y = est, color = state)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 1) +
  theme_bw() +
  xlab(" ") +
  ylab("Initial state probability") +
  scale_color_discrete(name = " ")  +
  theme(legend.position = "none")

ggsave("output/hmm_tfm_isp_2021-07-19.tiff", compression= "lzw", 
       scale = 0.6, height = 5, width = 5, dpi = 400)

# Distributions by state ####

# Step length

# Get canonical parameters of gamma distr from mean and variance
# k is shape
# theta is scale 
# mean = k*theta
# variance = k*theta^2
# theta = mean/variance
# k = mean/theta

# State 1
theta1 <- hmm_tfm3$mle$step[1,1]/hmm_tfm3$mle$step[2,1]
k1 <- hmm_tfm3$mle$step[1,1]/theta1

# State 2
theta2 <- hmm_tfm3$mle$step[1,2]/hmm_tfm3$mle$step[2,2]
k2 <- hmm_tfm3$mle$step[1,2]/theta2

# State 3
theta3 <- hmm_tfm3$mle$step[1,3]/hmm_tfm3$mle$step[2,3]
k3 <- hmm_tfm3$mle$step[1,3]/theta3

# Simulate and plot
sim_step <- data.frame(sim = c(rgamma(n = 10000, shape = k1, scale = theta1),
                               rgamma(n = 10000, shape = k2, scale = theta2),
                               rgamma(n = 10000, shape = k3, scale = theta3)),
                       state = rep(c("Resting", "Traveling", "Searching"), 
                                   each = 10000)) %>% 
  mutate(state = factor(state, levels = c("Resting", "Searching", "Traveling")))

ggplot(sim_step, aes(x = sim, fill = state)) +
  geom_histogram(binwidth = 30, color = "black") +
  facet_wrap(~ state) +
  labs(x = "Step length (m)", y = "Frequency") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/hmm_tfm3_sl_2021-07-19.tiff", compression= "lzw", 
       scale = 0.6, height = 5, width = 12, dpi = 400)

ggplot(sim_step, aes(x = sim, fill = state)) +
  geom_histogram(binwidth = 30, color = "black") +
  labs(x = "Step length (m)", y = "Frequency", fill = "State") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/hmm_tfm3_sl_2021-07-19_alt.tiff", compression= "lzw", 
       scale = 0.6, height = 5, width = 8, dpi = 400)

# Turning angles

# Simulate and plot
sim_angle <- data.frame(sim = c(rvm(n = 10000, mean = hmm_tfm3$mle$angle[1,1], 
                                    k = hmm_tfm3$mle$angle[2,1]),
                                rvm(n = 10000, mean = hmm_tfm3$mle$angle[1,2], 
                                    k = hmm_tfm3$mle$angle[2,2]),
                                rvm(n = 10000, mean = hmm_tfm3$mle$angle[1,3], 
                                    k = hmm_tfm3$mle$angle[2,3])),
                        state = rep(c("Resting", "Traveling", "Searching"), 
                                    each = 10000)) %>% 
  mutate(state = factor(state, levels = c("Resting", "Searching", "Traveling")))

ggplot(sim_angle, mapping = aes(x = sim, fill = state)) +
  stat_bin(breaks = seq(0, 2*pi, by = pi/12), color = "black") +
  facet_wrap(~ state) +
  scale_x_continuous(breaks = c(0, pi/2, pi, -pi/2), 
                     labels = c("0", "90", "180", "270")) +
  labs(y = "Frequency", x = "Turning angle (degrees)") +
  coord_polar(start = 0) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/hmm_tfm3_ta_2021-07-19.tiff", compression= "lzw", 
       scale = 0.6, height = 5, width = 10, dpi = 400)

# Mean distance travelled per night ####

hmm_data %>% 
  group_by(sex, ID) %>% 
  summarize(dist_tr = sum(step, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  summarize(mean_dist_tr = mean(dist_tr),
            sd_dist_tr = sd(dist_tr),
            min_dist_tr = min(dist_tr),
            max_dist_tr = max(dist_tr))
