# In our meeting in August 2021, we realized we had been calculating distance to 
# settlements wrong (from the boundary of a settlement, without differentiating
# between inside or outside the boundary). Here, I recalculate it making sure 
# that distance to settlement is 0 if the animal is within a settlement and the
# distance to the boundary otherwise. 

# Load packages ####

library(tidyverse)
library(sf)

# Load data ####

dat <- read.csv("output/hmm_3states_time-from-midnight_2021-07-19_3projections.csv")

# Load shapefile of settlements ####

sett <- read_sf("input/Settlements.shp")

# Transform to spatial points ####

pts <- st_as_sf(dat, coords = c("itm_x", "itm_y"), na.fail = FALSE, 
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

# Save ####

dat$dist_to_sett <- pts$dist_to_sett

write.csv(dat, "output/hmm_3states_time-from-midnight_2021-07-19_3projections_dist-to-settlements.csv")
