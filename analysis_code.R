# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code to reproduce the analysis from Bar-Ziv, E., 
# Picardi, S., Kaplan, A., Avgar, T., & Berger-Tal, O. (2022). 
# Sex differences dictate the movement patterns of Striped Hyenas, 
# Hyaena hyaena, in a human-dominated landscape. 
# Frontiers in Ecology and Evolution, 10, 897132.
# # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # 

# Load packages ####

library(tidyverse)
library(sf)
library(lubridate)
library(momentuHMM)
library(CircStats)
library(readxl)
library(mclogit)
library(MASS)
library(patchwork)

# Set seed ####

set.seed(7)

# Load data ####

# This file includes data from the Central District
# Sent by Einat on 9/21/2020 via Google Drive
# I saved a .csv copy of it

hcd <- read.csv("input/Collars_Center - Copy.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::select(id, date, time, latitude, longitude)

# Make timestamp column ####

hcd <- hcd %>% 
  # If I try to assign timezone, I have a problem with 4 dates that fall 
  # right at time change on daylight savings. Fix these manually first:
  mutate(timestamp = mdy_hm(paste(date, time), tz = "UTC")) %>% 
  mutate(timestamp2 = case_when(
    year(timestamp) == 2020 & 
      month(timestamp) == 3 & 
      day(timestamp) == 27 &
      hour(timestamp) == 2 ~ ymd_hm("2020-03-27 01:59"),
    year(timestamp) == 2019 & 
      month(timestamp) == 3 & 
      day(timestamp) == 29 &
      hour(timestamp) == 2 ~ ymd_hm("2019-03-29 01:59"),    
    year(timestamp) == 2018 & 
      month(timestamp) == 3 & 
      day(timestamp) == 23 &
      hour(timestamp) == 2  ~ ymd_hm("2018-03-23 01:59"),
    TRUE ~ timestamp
  )) %>% 
  mutate(timestamp3 = ymd_hms(as.character(timestamp2), tz = "Israel")) %>% 
  dplyr::select(id, date, time, timestamp = timestamp3, latitude, longitude)

# Convert to UTM ####

hcd <- st_as_sf(hcd, coords = c("longitude", "latitude"))
st_crs(hcd) <- 4326
hcd <- st_transform(hcd, crs = 32636) 
hcd_xy <- st_coordinates(hcd)
hcd <- hcd %>% 
  as.data.frame() %>% 
  cbind.data.frame(hcd_xy) %>% 
  dplyr::select(-geometry) %>% 
  rename(x = X, y = Y)

# Remove duplicates ####

dups <- hcd %>% 
  group_by(id, timestamp) %>% 
  tally() %>% 
  filter(n > 1) 

hcd_pt1 <- hcd %>% 
  group_by(id, timestamp) %>% 
  slice(2)

# Non-duplicated locations:

non_dups <- hcd %>% 
  group_by(id, timestamp) %>% 
  tally() %>% 
  filter(n == 1) %>% 
  mutate(key = paste(id, timestamp))

hcd_pt2 <- hcd %>% 
  mutate(key = paste(id, timestamp)) %>% 
  filter(key %in% non_dups$key) %>% 
  dplyr::select(-key)

# Bind 

hcd <- bind_rows(hcd_pt1, hcd_pt2) %>% 
  group_by(id) %>% 
  arrange(timestamp, by_group = TRUE) %>% 
  ungroup()

# Speed filter ### 

# When plotting the raw data (see below), there's an obvious outlier that
# I want to get rid of. Use a speed filter:

adehabitatLT::as.ltraj(hcd[, c("x", "y")], 
                              date = hcd$timestamp, 
                              id = hcd$id) %>% 
  adehabitatLT::ld() %>% 
  summarize(min(dist, na.rm = TRUE), 
            mean(dist, na.rm = TRUE), 
            median(dist, na.rm = TRUE), 
            max(dist, na.rm = TRUE))

hcd <- adehabitatLT::as.ltraj(hcd[, c("x", "y")], 
                       date = hcd$timestamp, 
                       id = hcd$id) %>% 
  adehabitatLT::ld() %>% 
  filter(dist < 200000) %>% 
  dplyr::select(id, timestamp = date, x, y)

# saveRDS(hcd, "output/hyenas_data_cleaned.rds")

# Prep data for HMM ####

# Round timestamps and add ID
hyenas <- hcd %>% 
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

# hyenas_night$burst_[1] <- 1
#
# for (i in 2:nrow(hyenas_night)) {
#   if (hyenas_night$ID[i-1] != hyenas_night$ID[i] | hyenas_night$lag[i] > 1) {
#     hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] + 1
#   } else if (hyenas_night$lag[i] == 1) {
#     hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] 
#   }
# }

# Update August 2024: the code above did not seem to reproduce the old results.
# Something might have changed in the way dates are handled. 
# Use a different approach: create a fake timestamp using a timezone where the 
# time a "hyena-night" starts in Israel (6 PM) is midnight. This way I can
# assign the burst ID based on the date only. 
hyenas_night <- hyenas_night %>% 
  mutate(timestamp_fake = with_tz(timestamp, "Asia/Tokyo"),
         date_fake = as_date(timestamp_fake),
         date_id = format(date_fake, "%j")) %>% 
  mutate(burst_ = paste0(id, "_", date_id)) %>% 
  dplyr::select(-timestamp_fake, -date_fake, -date_id)

# Get rid of bursts that are too short to fit a model
drop <- hyenas_night %>% 
  group_by(burst_) %>% 
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
  
  ts <- seq(ranges[ranges$ID == who,]$start, ranges[ranges$ID == who, ]$end, by = 3600)
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

# ggsave("output/movement-summaries_2021-08-24.tiff",
#        compression = "lzw", width = 10, height = 5, dpi = 300, scale = 0.7)

# Fit HMM ####

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

#saveRDS(hmm_tfm3, "output/hmm_tfm3.rds")

plot(hmm_tfm3)

diag(hmm_tfm3$mod$Sigma)
any(diag(hmm_tfm3$mod$Sigma) < 0)
sum(diag(hmm_tfm3$mod$Sigma) < 0)

results_hmm3_tfm <- hmm_data %>% 
  mutate(state = viterbi(hmm_tfm3)) %>% 
  dplyr::select(ind_night = ID, id, timestamp, x, y, state)

#write.csv(results_hmm3_tfm, "output/hmm_3states_time-from-midnight_2021-07-19.csv")

# Convert to ITM ####

dat <- results_hmm3_tfm %>% 
  filter(!is.na(x))

dat_sp <- st_as_sf(dat, coords = c("x", "y"))
st_crs(dat_sp) <- 32636

dat_itm <- st_transform(dat_sp, crs = 2039) 
xy_itm <- st_coordinates(dat_itm)
colnames(xy_itm) <- c("itm_x", "itm_y")

dat_ll <- st_transform(dat_sp, crs = 4326) 
xy_ll <- st_coordinates(dat_ll)
colnames(xy_ll) <- c("long", "lat")

dat_3p <- dat %>% 
  rename(utm_x = x, utm_y = y) %>% 
  cbind(xy_itm) %>% 
  cbind(xy_ll) %>% 
  relocate(state, .after = lat)

# write.csv(dat_3p, "output/hmm_3states_time-from-midnight_2021-07-19_3projections.csv")

# Calculate distance to settlements ####

# Load shapefile of settlements ####

sett <- read_sf("input/Settlements.shp")

# Transform to spatial points ####

pts <- st_as_sf(dat_3p, coords = c("itm_x", "itm_y"), na.fail = FALSE, 
                crs = "+init=epsg:2039")

# Plot ####

ggplot(sett) +
  geom_sf() +
  geom_sf(aes(color = id), data = (pts))

# Calculate intersection ####

inter <- st_intersects(pts, sett, sparse = FALSE)

pts$within_sett <- apply(inter, 1, any)

ggplot() +
  geom_sf(data = sett) +
  geom_sf(aes(color = within_sett), data = (pts)) 

# Calculate distance ####

# st_distance already assigns 0 when a point is within a polygon so the chunk
# above is not actually needed. 

dista <- st_distance(pts, sett)

pts$dist_to_sett <- apply(dista, 1, min)

dat_3p$within_sett <- pts$within_sett
dat_3p$dist_to_sett <- pts$dist_to_sett

#write.csv(dat, "output/hmm_3states_time-from-midnight_2021-07-19_3projections_dist-to-settlements.csv")

# Multinomial regression ####

# Data intersected with distance to settlement (use this as the main data!)
dat_sett <- dat_3p %>% 
  drop_na(itm_x, itm_y)

# Data intersected with distance to roads (done by Einat) 
dat_roads <- read.csv("input/2021-08/HMM3States_time-from-midnight_190721_Disttoroads.csv") %>% 
  mutate(date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -FID, -X, -NEAR_FID, -ind_night) %>% 
  rename(dist_roads = NEAR_DIST)

# Data intersected with agricultural areas (done by Einat) 
dat_ag <- read_xls("input/2021-08/Agriculture.xls") %>% 
  mutate(in_ag = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with JNFF (done by Einat) 
dat_jnff <- read_xls("input/2021-08/JNFF.xls") %>% 
  # filter out points that are in olive trees or fruit trees
  filter(!COV_TYPE %in% c(2994, 2990)) %>% 
  mutate(in_jnff = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with JNFV (done by Einat) 
dat_jnfv <- read_xls("input/2021-08/JNFV.xls") %>% 
  mutate(in_jnfv = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with JNF campsites (done by Einat) 
dat_jnfc <- read_xls("input/2021-08/JNFCampsites_Buffer50.xls") %>% 
  mutate(in_jnfc = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with land cover (done by Einat) 
dat_lc <- read_xls("input/2021-08/Landcover.xls") %>% 
  mutate(date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with livestock facilities (done by Einat) 
dat_live <- read_xls("input/2021-08/Livestock_Buffer50.xls") %>% 
  mutate(livestock = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Data intersected with poultry facilities (done by Einat) 
dat_poul <- read_xls("input/2021-08/Poultry_Buffer50.xls") %>% 
  mutate(poultry = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -ind_night)

# Individual info 
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
  dplyr::select(id, district = location, sex, sex_rs, legs) %>% 
  # 6168 is a male, not a female!
  mutate(sex = case_when(
    id == "6168" ~ "Male",
    TRUE ~ sex
  )) %>% 
  mutate(sex_rs = case_when(
    id == "6168" ~ "M",
    TRUE ~ sex_rs
  ))

# Join all together ####

dat <- dat_sett %>% 
  mutate(itm_x = round(itm_x, digits = 4),
         itm_y = round(itm_y, digits = 4),
         date = as_date(timestamp)) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date) %>% 
  left_join(dat_roads, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  left_join(dat_ag, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag) %>% 
  left_join(dat_jnfc, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc) %>% 
  left_join(dat_jnff, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff) %>% 
  left_join(dat_jnfv, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv) %>% 
  left_join(dat_live, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv, livestock) %>% 
  left_join(dat_poul, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv, livestock, poultry) %>% 
  left_join(dat_lc, by = c("id", "itm_x", "itm_y", "state", "date")) %>% 
  left_join(info, by = "id")

# Reclassify land cover ####

# Add readable description for FCODE
dat <- dat %>% 
  mutate(fcode_descr = case_when(
    FCODE == 501 ~ "Artificial",
    FCODE == 502 ~ "Nat no veg",
    FCODE == 503 ~ "Ag Lawn Garden",
    FCODE == 504 ~ "Nat veg",
    FCODE == 0 ~ NA_character_,
    is.na(FCODE) ~ NA_character_
  ))

# Add readable description for FTYPE
dat <- dat %>% 
  mutate(ftype_descr = case_when(
    FTYPE == 11 ~ "Builtup",
    FTYPE == 31 ~ "Rocky surface",
    FTYPE == 37 ~ "Alluvium",
    FTYPE == 51 ~ "Cultivated",
    FTYPE == 52 ~ "Uncultivated",
    FTYPE == 53 ~ "Fruit trees",
    FTYPE == 54 ~ "Palm trees",
    FTYPE == 56 ~ "Olive trees",
    FTYPE == 57 ~ "Vineyards",
    FTYPE == 60 ~ "Lawn or garden",
    FTYPE == 78 ~ "Herbaceous",
    FTYPE == 79 ~ "Shrub",
    FTYPE == 81 ~ "Med low density short",
    FTYPE == 82 ~ "Med medium density short",
    FTYPE == 83 ~ "Med high density short",
    FTYPE == 84 ~ "Med low density tall",
    FTYPE == 85 ~ "Med medium density tall",
    FTYPE == 86 ~ "Med high density tall",
    FTYPE == 0 ~ NA_character_
  ))

# Create one single column for JNF
dat <- dat %>% 
  mutate(jnf = case_when(
    is.na(in_jnff) ~ 0,
    TRUE ~ 1
  ))

# Merge livestock and poultry into one
dat <- dat %>% 
  mutate(livestock_poultry = case_when(
    livestock == 1 ~ 1,
    poultry == 1 ~ 1,
    TRUE ~ 0))

# Modify distance to roads into binary variable 
# (1 if within 100m from a road, 0 otherwise)
dat <- dat %>% 
  mutate(roads_yn = case_when(
    dat$dist_roads <= 50 ~ 1,
    TRUE ~ 0
  ))

# Create land cover
dat <- dat %>% 
  mutate(land_cover_simple = case_when(
    # First, any points within settlements get assigned to anthropogenic
    within_sett ~ "Anthropogenic",
    # Then, any points labeled as artificial get assigned to antrhopogenic
    fcode_descr == "Artificial" ~ "Anthropogenic",
    # Then, any points in lawns/builtup get assigned to anthropogenic
    ftype_descr %in% c("Lawn or garden", "Builtup") ~ "Anthropogenic",
    # Then also any points within 50 m of a livestock and poultry facility
    # get assigned to anthropogenic
    livestock_poultry == 1 ~ "Anthropogenic",
    # Then any points within 50 m of JNF campsites also get assigned to anthropogenic
    !is.na(in_jnfc) ~ "Anthropogenic",
    # Then any points within 50 m of a road get assigned to roads
    roads_yn == 1 ~ "Roads",
    # Any cultivated surface gets assigned to agriculture
    ftype_descr %in% c("Vineyards", "Fruit trees", "Olive trees", 
                       "Palm trees", "Cultivated", "Uncultivated") ~ "Agriculture",
    # Any natural surface gets assigned to natural
    fcode_descr %in% c("Nat veg", "Nat no veg") ~ "Natural",
    # If there are remaining points in JNF that haven't been assigned yet,
    # they get assigned to natural
    jnf == 1 ~ "Natural",
    TRUE ~ NA_character_
  )) %>% 
  mutate(state = factor(case_when(
    state == 1 ~ "Resting",
    state == 2 ~ "Traveling",
    state == 3 ~ "Searching"
  ))) %>% 
  mutate(log_dist_settlement = log(dist_to_sett + 0.001),
         id = factor(id),
         land_cover_simple = factor(land_cover_simple))

#write.csv(dat, "output/multinomial-regression-data_2021-08-30.csv")

# Multinomial regression ####

mn_dat <- dat

mn_dat$state <- factor(mn_dat$state)
mn_dat$sex <- factor(mn_dat$sex)
mn_dat$land_cover_simple <- factor(mn_dat$land_cover_simple)

mod1 <- mblogit(state ~ land_cover_simple * sex,
                data = mn_dat,
                random = ~ 1|id)

mod2 <- mblogit(state ~ (land_cover_simple +
                           log_dist_settlement) * sex,
                data = mn_dat,
                random = ~ 1|id)

newd_land <- data.frame(land_cover_simple = factor(rep(c("Agriculture", "Anthropogenic", "Natural", "Roads"), 2)),
                        log_dist_settlement = log(mean(mn_dat$dist_to_sett, na.rm = T)),
                        roads_yn = factor(0, levels = c(0, 1)),
                        sex = factor(rep(c("Female", "Male"), each = 4)),
                        id =  factor("6163A"),
                        ind_night = factor(2))

newd_settlement <- data.frame(land_cover_simple = factor("Natural", levels = levels(mn_dat$land_cover_simple)),
                              log_dist_settlement = rep(log(seq(0, 2500, length.out = 100) + 0.001), 2),
                              dist_to_sett = rep(seq(0, 2500, length.out = 100), 2),
                              roads_yn = factor(0, levels = c(0, 1)),
                              sex = factor(rep(c("Female", "Male"), each = 100)),
                              id = factor("6163A"),
                              ind_night = factor(2))

preds_land <- predict(mod1, newdata = newd_land, conditional = FALSE, type = "response")
preds_settlement <- predict(mod1, newdata = newd_settlement, conditional = FALSE, type = "response")

# Bootstrap confidence intervals ####

boot_land <- array(dim = c(nrow(preds_land), 3, 1000))
boot_sett <- array(dim = c(nrow(preds_settlement), 3, 1000))

for (i in 1:1000) {
  
  print(i)
  
  samp <- mn_dat[sample(1:nrow(mn_dat), replace = TRUE), ]
  
  mod <- mblogit(state ~ (land_cover_simple +
                            log_dist_settlement) * sex,
                 data = samp,
                 random = ~ 1|id)
  
  boot_land[, , i] <- predict(mod, newdata = newd_land, conditional = FALSE, type = "response")
  boot_sett[, , i] <- predict(mod, newdata = newd_settlement, conditional = FALSE, type = "response")
  
}

means_land <- apply(boot_land, c(1, 2), FUN = mean)
means_sett <- apply(boot_sett, c(1, 2), FUN = mean)

lwr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.025)
lwr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.025)

upr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.975)
upr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.975)

dist_lwr_lands <- means_land - lwr_land
dist_lwr_sett <- means_sett - lwr_sett

dist_upr_lands <- upr_land - means_land
dist_upr_sett <- upr_sett - means_sett

lwr_lands <- preds_land - dist_lwr_lands
colnames(lwr_lands) <- paste0(colnames(lwr_lands), "_lwr")
lwr_sett <- preds_settlement - dist_lwr_sett
colnames(lwr_sett) <- paste0(colnames(lwr_sett), "_lwr")

upr_lands <- preds_land + dist_upr_lands
colnames(upr_lands) <- paste0(colnames(upr_lands), "_upr")
upr_sett <- preds_settlement + dist_upr_sett
colnames(upr_sett) <- paste0(colnames(upr_sett), "_upr")

# Plot ####

res_land <- cbind.data.frame(newd_land, preds_land) %>% 
  pivot_longer(cols = Resting:Traveling, names_to = "state", values_to = "pred")
res_settlement <- cbind.data.frame(newd_settlement, preds_settlement)  %>% 
  pivot_longer(cols = Resting:Traveling, names_to = "state", values_to = "pred")

upr_lands <- upr_lands %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_upr:Traveling_upr, names_to = "state", values_to = "upr") %>% 
  dplyr::select(upr)
upr_sett <- upr_sett %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_upr:Traveling_upr, names_to = "state", values_to = "upr") %>% 
  dplyr::select(upr)

lwr_lands <- lwr_lands %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_lwr:Traveling_lwr, names_to = "state", values_to = "lwr") %>% 
  dplyr::select(lwr)
lwr_sett <- lwr_sett %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_lwr:Traveling_lwr, names_to = "state", values_to = "lwr") %>% 
  dplyr::select(lwr)

res_land <- res_land %>% 
  cbind.data.frame(upr_lands) %>% 
  cbind.data.frame(lwr_lands)
res_settlement <- res_settlement %>% 
  cbind.data.frame(upr_sett) %>% 
  cbind.data.frame(lwr_sett) %>% 
  mutate(dist_settlement = exp(log_dist_settlement) - 0.001)

saveRDS(res_land, "output/bootstrapped_multinomial_2022-02-07_response-to-land-cover.rds")
saveRDS(res_settlement, "output/bootstrapped_multinomial_2022-02-07_response-to-settlement.rds")

ggplot(res_land, aes(x = land_cover_simple, y = pred, color = state)) +
  geom_point(size = 2, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, width = 0.2), position = position_dodge(0.2)) +
  facet_wrap(~ sex) +
  labs(x = " ", y = "Probability", color = "State") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_2022-02-07_land-cover_sex.tiff",
       compression = "lzw", width = 10, height = 5, dpi = 300)

res_land %>% 
  filter(land_cover_simple == "Anthropogenic") %>% 
  ggplot(aes(x = sex, y = pred, color = state)) +
  geom_point(size = 2, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, width = 0.2), position = position_dodge(0.2)) +
  labs(x = " ", y = "Probability", color = "State", 
       title = "Activity allocation in anthropogenic areas") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_2023-02-05_simplified.tiff",
       compression = "lzw", width = 4, height = 3, dpi = 300)

ggplot(res_settlement, aes(x = dist_to_sett, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_wrap(~ sex) +
  labs(x = "Distance to settlements (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 500))

ggsave("output/multinomial-model-predictions_2022-02-07_distance-to-settlement.tiff",
       compression = "lzw", width = 6, height = 4, dpi = 300)

# Keep x axis on log scale for distance to settlements
ggplot(res_settlement, aes(x = log_dist_settlement, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_wrap(~ sex) +
  labs(x = "Distance to settlements (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom", axis.title.x.top = element_text(margin = margin(1.5))) +
  scale_x_continuous(breaks = log(c(0, 10, 100, 1000) + 0.001), 
                     labels = c(0, 10, 100, 1000)) 

ggsave("output/multinomial-model-predictions_2022-02-07_log_distance-to-settlement.tiff",
       compression = "lzw", width = 6, height = 4, dpi = 300)

# Summary plots ####

# Proportion of steps in each land cover type for each sex 

# Broken down by individuals and then average:

steps_total <- mn_dat %>% 
  group_by(id) %>% 
  tally() %>% 
  rename(tot = n)

steps_per_cat <- mn_dat %>%
  group_by(id, land_cover_simple, sex) %>% 
  tally() %>% 
  left_join(steps_total, by = "id") %>% 
  mutate(prop = n/tot) %>% 
  group_by(land_cover_simple, sex) %>% 
  summarize(avg_prop = mean(prop),
            sd_prop = sd(prop)) %>% 
  filter(!is.na(land_cover_simple))

# Chisquare test
s_lc_mat <- mn_dat %>%
  group_by(land_cover_simple, sex) %>% 
  tally() %>% 
  filter(!is.na(land_cover_simple)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = sex, values_from = n) %>% 
  dplyr::select(Female, Male) %>% 
  as.matrix()

rownames(s_lc_mat) <- c("Agriculture",
                        "Anthropogenic",
                        "Natural",
                        "Roads")

chisq.test(s_lc_mat)

# Plot
ggplot(steps_per_cat, aes(x = sex, y = avg_prop, fill = land_cover_simple)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = " ", y = "Proportion of steps", fill = "Land cover") +
  scale_fill_brewer(palette = "Set2")

ggsave("output/steps-by-land-cover_2022-02-07.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300, scale = 0.7)

# Without taking the average across individuals, just crude numbers:

steps_total <- mn_dat %>% 
  group_by(sex) %>% 
  tally() %>% 
  rename(tot = n)

steps_per_cat <- mn_dat %>%
  group_by(land_cover_simple, sex) %>% 
  tally() %>% 
  left_join(steps_total, by = "sex") %>% 
  mutate(prop = n/tot) %>% 
  filter(!is.na(land_cover_simple))

ggplot(steps_per_cat, aes(x = sex, y = prop, fill = land_cover_simple)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = " ", y = "Proportion of steps", fill = "Land cover") +
  scale_fill_brewer(palette = "Set2")

ggsave("output/steps-by-land-cover-RAW_2022-02-07.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300, scale = 0.7)

# Proportion of steps in each behavioral state for each sex 

# Broken down by individuals and then average:

steps_total <- mn_dat %>% 
  group_by(id) %>% 
  tally() %>% 
  rename(tot = n)

steps_per_state <- mn_dat %>%
  group_by(id, state, sex) %>% 
  tally() %>% 
  left_join(steps_total, by = "id") %>% 
  mutate(prop = n/tot) %>% 
  group_by(state, sex) %>% 
  summarize(avg_prop = mean(prop))

ggplot(steps_per_state, aes(x = sex, y = avg_prop, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = " ", y = "Proportion of steps", fill = "State")

ggsave("output/steps-by-behavioral-state_2022-02-07.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300, scale = 0.7)

# Without taking the average across individuals, just crude numbers:

steps_total <- mn_dat %>% 
  group_by(sex) %>% 
  tally() %>% 
  rename(tot = n)

steps_per_state <- mn_dat %>%
  group_by(state, sex) %>% 
  tally() %>% 
  left_join(steps_total, by = "sex") %>% 
  mutate(prop = n/tot) 

ggplot(steps_per_state, aes(x = sex, y = prop, fill = factor(state))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = " ", y = "Proportion of steps", fill = "State")

ggsave("output/steps-by-behavioral-state-RAW_2022-02-07.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300, scale = 0.7)
