# # # # # # # # # # # # # # # # # # # # # # # # # 
# CTMM simulation of striped hyena movement data
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Nov 13, 2020 # # # # # # # # # # # # #
# Last updated Nov 16, 2020 # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# The approach I used in hyenas_ctmm.R to regularize the hyena data to 1h 
# resolution results in linearly interpolated tracks instead of a random 
# realization of locations based on the underlying movement model. Here,
# I try a similar but slightly different approach using ctmm::simulate.ctmm
# that should yield a random realization of regular individual tracks based
# on the observed data. 

# Load packages ####

library(tidyverse)
library(move)
library(ctmm)
library(lubridate)

# Load cleaned data ####

hyenas <- readRDS("output/hyenas_data_cleaned.rds")

# Select one test individual ####

# Try on one individual
test_raw <- hyenas %>% 
  as.data.frame() %>% 
  filter(id == unique(hyenas$id)[1])

# Mimic Movebank format ####

# ctmm works with 'telemetry' objects. To coerce data into a 'telemetry' 
# object, I need it to mimic the Movebank format. This means changing 
# the names and using lat-long instead of UTMs. 

coordinates(test_raw) <- c("x", "y")
crs(test_raw) <- crs("+init=epsg:32636")
test_raw_ll <- spTransform(test_raw, crs("+init=epsg:4326")) %>% 
  as.data.frame

names(test_raw_ll) <- c("individual.local.identifier",
                   "timestamp",
                   "location.long",
                   "location.lat")

# Process timestamps ####

# Round timestamps to the hour and discard 20-min resolution locations
test_raw_ll <- test_raw_ll %>% 
  filter(minute(timestamp) %in% c(0:3, 57:59)) %>% 
  mutate(timestamp_original = timestamp) %>% 
  mutate(timestamp = round(timestamp, units = "hours"))

# Remove duplicates
test_raw_ll <- test_raw_ll[!duplicated(test_raw_ll$timestamp), ]

# Coerce into 'telemetry' object ####

h_telemetry <- as.telemetry(object = test_raw_ll, 
             timeformat = "%Y-%m-%d %H:%M:%S",
             timezone = "Israel",
             projection = crs("+init=epsg:32636"))

# Fit CTMM ####

# Make a guess for initial parameter values
ctmm.guess(data = h_telemetry, name = "guess1", interactive = TRUE)

# Fit different possible models
mod_sel <- ctmm.select(data = h_telemetry, CTMM = guess1, verbose = TRUE)

# Choose best model (the function automatically ranks them)
best_mod <- mod_sel[[1]]

# Simulate CTMM tracks ####

# I want to simulate tracks at 1h resolution. 
# There are three ways to do this (arguments t, dt, res) but the documentation
# for dt is really poor and res doesn't really do what I want. So I'll go with
# the t option. I need to make a vector of timestamps I want to simulate over.
# From the minimum timestamp to the maximum, by intervals of 1h (3600s).
ts <- seq(min(test_raw_ll$timestamp), max(test_raw_ll$timestamp), by = 3600)

# Expected number of data points:
length(ts)

# Simulate track
test1 <- ctmm::simulate(object = best_mod, data = h_telemetry, res = 1)
test2 <- ctmm::simulate(object = best_mod, data = h_telemetry, res = 1)

test3 <- ctmm::predict(object = best_mod, data = h_telemetry, t = ts)
test4 <- ctmm::predict(object = best_mod, data = h_telemetry, t = ts)

# Number of data points I get:
nrow(test1)

# Plot
ggplot(as.data.frame(test_raw), aes(x = x, y = y)) +
  geom_point(data = as.data.frame(test1), color = "tomato", alpha = 0.3) +
  geom_path(data = as.data.frame(test1), color = "tomato", alpha = 0.3) +
  geom_point(data = as.data.frame(test2), color = "steelblue", alpha = 0.3) +
  geom_path(data = as.data.frame(test2), color = "steelblue", alpha = 0.3) +
  geom_point(alpha = 0.3) +
  geom_path(alpha = 0.3) 

# Re-format timestamps ####

# The timestamp was transformed into numeric by CTMM. I want the timestamp 
# back like it was.

convert_ts <- data.frame(ts = ts, ts_num = as.numeric(ts))

sim_track_df <- as.data.frame(sim_track@.Data) 

names(sim_track_df) <- c("ts_num", "x", "y", "vx", "vy")

sim_track_df <- left_join(sim_track_df, convert_ts, by = "ts_num")

# Check that the process is stochastic ####

# Try and simulate once more. If the process is stochastic, the results of the
# two simulations should be different 
test4 <- ctmm::simulate(object = best_mod, data = h_telemetry[[1]], t = ts)
sim_track2 <- test4[test4$t %in% as.numeric(ts),]

# Create manual legend for plot

colors <- c("Raw data" = "black", 
            "Simulation 1" = "tomato", 
            "Simulation 2" = "steelblue")

# Plot 
ggplot(test_raw, aes(x = x, y = y)) +
  geom_point(aes(color = "Simulation 1"), data = as.data.frame(sim_track), alpha = 0.3) +
  geom_path(aes(color = "Simulation 1"), data = as.data.frame(sim_track), alpha = 0.3) +
  geom_point(aes(color = "Simulation 2"), data = as.data.frame(sim_track2), alpha = 0.3) +
  geom_path(aes(color = "Simulation 2"), data = as.data.frame(sim_track2), alpha = 0.3) +
  geom_point(aes(color = "Raw data"), alpha = 0.3) +
  geom_path(aes(color = "Raw data"), alpha = 0.3) +
  scale_color_manual(values = colors) +
  labs(x = "UTM Easting", y = "UTM Northing", 
       title = "CTMM simulations for one individual track",
       color = " ") +
  theme_bw()

ggsave("output/ctmm_simulate_example.tiff", compression = "lzw", dpi = 400,
       width = 8, height = 6)