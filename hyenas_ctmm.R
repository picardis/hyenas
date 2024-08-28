# # # # # # # # # # # # # # # # # # # # # # # # # 
# CTMM resampling of striped hyena movement data
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Sept 28, 2020 # # # # # # # # # # # # #
# Last updated Nov 17, 2020 # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# Load packages ####

library(tidyverse)
library(sp)
library(lubridate)
library(momentuHMM)

# Set seed ####

set.seed(7)

# Load data ####

# This file includes data from the Central District
# Sent by Einat on 9/21/2020 via Google Drive
# I saved a .csv copy of it

hcd <- read.csv("input/Collars_Center - Copy.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::select(id = i_id, date, time, latitude, longitude)

# This file includes data from both the Central and the South District
# Sent by Einat on 9/30/2020 via Google Drive
# I saved a .csv copy of it
# I'm still calling this hcd so I don't have to change it downstream

hcd <- read.csv("input/South_Center - Copy.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::select(area, id, date, time, latitude, longitude)

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

coordinates(hcd) <- c("longitude", "latitude")
proj4string(hcd) <- CRS("+init=epsg:4326")
hcd <- spTransform(hcd, CRS("+init=epsg:32636")) %>% 
  as.data.frame() %>% 
  rename(x = longitude, y = latitude)

# Remove duplicates ####

dups <- hcd %>% 
  group_by(id, timestamp) %>% 
  tally() %>% 
  filter(n > 1) # 12 duplicates

# There's only 12 duplicates over 27666 locations and we're going to resample
# using the continuous time movement model (CTMM) anyway, so which one of the
# duplicated fixes I choose is not a big deal. I'll choose the second one 
# because sometimes that's the fix the tag takes after it resets the clock,
# and it has more time to take a high-quality fix. 

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

saveRDS(hcd, "output/hyenas_data_cleaned.rds")

# Plot raw data ####

ggplot(hcd, aes(x = x, y = y, color = id)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  labs(title = "Raw data", x = "UTM Easting", y = "UTM Northing") +
  scale_color_discrete(name = "ID") 

ggsave("output/full_raw.tiff", compression= "lzw", 
       scale = 1.2, height = 7, width = 7, dpi = 400)

# Process data ####

# The GPS tags were programmed to collect data every hour during the night,
# one location in the middle of the day, and to occasionally increase the
# fix rate to 20 minutes. The HMM wants data to be at regular intervals.
# I can think of three ways to deal with this:
# 1. Fit a continuous-time movement model (CTMM) to the data and resample 
#   to a resolution of interest, e.g., 1-hour. The problem with this is that
#   the locations during the day will be highly interpolated unless we choose
#   a coarser resolution which would mean throwing out a lot of good data. 
# 2. Split the track into nightly bursts, keeping only the nighttime data. 
#   This has two consequences: first, every night, we would reset the initial
#   behavioral state, so instead of having one initial state per individual
#   we would have as many initial states as days of tracking; second, we 
#   may want to fit the HMM with only two states (restricted and traveling)
#   instead of three (resting, restricted, and traveling.) Depending on whether
#   we think that hyenas ONLY sleep during the day or they could potentially
#   rest at night, or at dawn or dusk. If I choose this route, I 
#   still need to run the CTMM first, because the HMM does not handle NAs, 
#   so I need to interpolate missing locations. If I try to fit the CTMM 
#   after I split the bursts, the CTMM will struggle to converge on such 
#   short bursts.
# 3. Replicate the daytime point to all unsampled daytime hours. If the 
#   assumption is that hyenas sleep during the day, we can impose that on
#   the data by assuming that the daytime location is the resting site and 
#   repeating it throughout the day, instead of letting the CTMM "make up"
#   some movements. If, ideally, the last location in the morning and the 
#   first in the evening were in the same spot as the daytime location, then
#   the CTMM should give us the same result, but that is unlikely to happen
#   both because of actual behavior and because of GPS error. 

# Fit CTMM ####

# Option 1 above is the CTMM.

# Format like momentuHMM wants it
hcd <- hcd %>% 
  mutate(ID = id)

# 1 hour time lag
hcd_ct1 <- crawlWrap(obsData = hcd,
                     timeStep = "1 hours",
                     coord = c("x", "y"), 
                     Time.name = "timestamp")

# Bind and reformat data
hcd_ct1_pred <- hcd_ct1$crwPredict %>% 
  as_tibble() %>% 
  dplyr::select(timestamp, id = ID, x = mu.x, y = mu.y) 

# Select only locations at 1 hour intervals (the predicted ones)
hcd_1h <- hcd_ct1_pred %>% 
  group_by(id) %>% 
  arrange(timestamp) %>% 
  mutate(orig = timestamp[1]) %>% 
  mutate(timediff = difftime(timestamp, orig, units = "secs")) %>% 
  mutate(flag = timediff/1/60/60 == floor(timediff/1/60/60)) %>% 
  filter(flag == TRUE) %>% 
  dplyr::select(timestamp, id, x, y) %>% 
  ungroup() %>% 
  arrange(id, timestamp)

# Plot
ggplot(hcd_1h, aes(x = x, y = y, color = id)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  labs(title = "Resampled data (1-hour CTMM)", 
       x = "UTM Easting", y = "UTM Northing") +
  scale_color_discrete(name = "ID")

ggsave("output/full_1hCTMM.tiff", compression= "lzw", 
       scale = 1.2, height = 7, width = 7, dpi = 400)

# Save
saveRDS(hcd_1h, "output/full_1hCTMM.rds")

# Plot CTMM-1h data on top of original data ####

ggplot(hcd, aes(x = x, y = y)) + 
  geom_path(show.legend = FALSE, color = "black") +
  geom_path(data = hcd_1h, show.legend = FALSE, color = "tomato", alpha = 0.5) +
  geom_point(show.legend = FALSE, color = "black") +
  geom_point(data = hcd_1h, show.legend = FALSE, color = "tomato", alpha = 0.5) +
  facet_wrap(~ id, scales = "free") +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing")

ggsave("output/full_raw_over_1hCTMM.tiff", compression= "lzw", 
       height = 20, width = 20, dpi = 400, scale = 0.5)

# Nightly bursts ####

# Option 2 above is to have nightly bursts. 

hcd_night <- hcd_1h %>% 
  # filter nighttime data only (6PM to 8AM)
  filter(hour(timestamp) > 17 | hour(timestamp) < 9) %>% 
  group_by(id) %>% 
  mutate(lag = as.numeric((timestamp) - lag(timestamp))) %>% 
  mutate(lag = case_when(
    is.na(lag) ~ 1,
    TRUE ~ lag
  )) %>% 
  mutate(burst_ = NA_real_)

hcd_night$burst_[1] <- 1

for (i in 2:nrow(hcd_night)) {
  if (hcd_night$id[i-1] != hcd_night$id[i] | hcd_night$lag[i] > 1) {
    hcd_night$burst_[i] <- hcd_night$burst_[i-1] + 1
  } else if (hcd_night$lag[i] == 1) {
    hcd_night$burst_[i] <- hcd_night$burst_[i-1] 
  }
}

hcd_night <- hcd_night %>% 
  dplyr::select(-lag) %>% 
  ungroup()

# Get rid of bursts that are too short to fit a model
drop <- hcd_night %>% 
  group_by(burst_) %>% 
  tally() %>% 
  filter(n < 15) # 15 locations expected per night

hcd_night <- hcd_night %>% 
  filter(!burst_ %in% drop$burst_)

# Save
saveRDS(hcd_night, "output/full_1hCTMM_night.rds")

# Plot night data on top of original data ####

hcd_night <- readRDS("output/full_1hCTMM_night.rds")
hcd <- readRDS("output/hyenas_data_cleaned.rds")

hcd <- hcd %>% 
  filter(id %in% c("6163A", "6163B", "6164F", 
                   "6167", "6168", "6169", "6170"))
hcd_night <- hcd_night %>% 
  filter(id %in% c("6163A", "6163B", "6164F", 
                   "6167", "6168", "6169", "6170"))

ggplot(hcd, aes(x = x, y = y)) + 
  geom_path(show.legend = FALSE, color = "black") +
  geom_path(data = hcd_night, show.legend = FALSE, color = "steelblue", alpha = 0.5) +
  geom_point(show.legend = FALSE, color = "black") +
  geom_point(data = hcd_night, show.legend = FALSE, color = "steelblue", alpha = 0.5) +
  facet_wrap(~ id, scales = "free") +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing")

ggsave("output/full_raw_over_1hCTMM_night_2021-09-21.tiff", compression= "lzw", 
       scale = 0.5, height = 20, width = 20, dpi = 400)
