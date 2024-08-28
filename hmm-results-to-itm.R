# Convert points in output from HMM into Israel Transverse Mercator CRS #

# Load packages ####

library(sp)
library(tidyverse)

# Load data ####

dat <- read.csv("output/hmm_3states_time-from-midnight_2021-07-19.csv")

# Convert to ITM ####

dat_sp <- dat %>% 
  filter(!is.na(x))
coordinates(dat_sp) <- c("x", "y")
proj4string(dat_sp) <- CRS("+init=epsg:32636")

dat_itm <- spTransform(dat_sp, CRS("+init=epsg:2039")) %>% 
  as.data.frame() %>% 
  rename(itm_x = x, itm_y = y)

dat_ll <- spTransform(dat_sp, CRS("+init=epsg:4326")) %>% 
  as.data.frame() %>% 
  rename(long = x, lat = y)

dat_3p <- dat %>% 
  rename(utm_x = x, utm_y = y) %>% 
  left_join(dat_itm) %>% 
  left_join(dat_ll) %>% 
  relocate(state, .after = lat)

write.csv(dat_3p, "output/hmm_3states_time-from-midnight_2021-07-19_3projections.csv")
