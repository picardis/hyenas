# Calculate time since/to midnight ####

# Make a variable for time since/to midnight. The value will be 0 in
# the middle of the night, negative between sunset and the middle of the night,
# and positive between the middle of the night and sunrise. Account for 
# variation of sunset/sunrise time through the lubridate::year; use the actual times when 
# the sun went up or down. 

# I want to use sunset for the current lubridate::day as a reference if before midnight, 
# but for the previous lubridate::day if after midnight. Conversely, I want to calculate 
# sunrise for the next lubridate::day if before midnight, or sunrise for the current lubridate::day 
# if before midnight. Essentially, I want to use noon instead of midnight as the 
# break between lubridate::days. 

calc_tfm <- function(x) {
  
  # x must be a data.frame with columns x, y (UTM coordinates), timestamp,
  # id, burst_
  # sp::coordinates(x) <- c("x", "y")
  # sp::proj4string(x) <- sp::CRS("+init=epsg:32636")
  
  x <- st_as_sf(x, coords = c("x", "y"))
  st_crs(x) <- 32636
  ll <- st_transform(x, crs = 4326) 
  ll <- st_coordinates(ll)
  x <- x %>% 
    as.data.frame() %>% 
    cbind.data.frame(ll) %>% 
    dplyr::select(-geometry) %>% 
    rename(x = X, y = Y)
  
  res <- x %>% 
    as.data.frame() %>% 
    mutate(
      # Establish reference dates for sunset and sunrise
      # Again, a problem with DST after subtracting 1 lubridate::day for set_date
      set_date = case_when(
        lubridate::year(timestamp) == 2020 & 
          lubridate::month(timestamp) == 3 & 
          lubridate::day(timestamp) == 28 &
          lubridate::hour(timestamp) == 2 ~ lubridate::ymd_hm("2020-03-27 01:59", tz = "Israel"),
        lubridate::year(timestamp) == 2019 & 
          lubridate::month(timestamp) == 3 & 
          lubridate::day(timestamp) == 30 &
          lubridate::hour(timestamp) == 2 ~ lubridate::ymd_hm("2019-03-29 01:59", tz = "Israel"),    
        lubridate::year(timestamp) == 2018 & 
          lubridate::month(timestamp) == 3 & 
          lubridate::day(timestamp) == 24 &
          lubridate::hour(timestamp) == 2  ~ lubridate::ymd_hm("2018-03-23 01:59", tz = "Israel"),
        lubridate::hour(timestamp) >= 12 ~ timestamp,
        lubridate::hour(timestamp) < 12 ~ timestamp - lubridate::days(1)),
      rise_date = case_when(
        lubridate::hour(timestamp) >= 12 ~ timestamp + lubridate::days(1),
        lubridate::hour(timestamp) < 12 ~ timestamp)) %>% 
    # Calculate sunrise and sunset times
    mutate(sunrise = maptools::sunriset(crds = as.matrix(.[, c("x", 
                                                               "y")]), 
                                        dateTime = rise_date,
                                        direction = "sunrise",
                                        POSIXct.out = TRUE)$time) %>%  
    mutate(sunset = maptools::sunriset(crds = as.matrix(.[, c("x", 
                                                              "y")]), 
                                       dateTime = set_date,
                                       direction = "sunset",
                                       POSIXct.out = TRUE)$time) %>% 
    # Calculate midnight
    mutate(midnight = sunset + (difftime(sunrise, sunset, tz = "Israel")/2)) %>% 
    # Count time since/to midnight
    mutate(time_from_midnight = difftime(timestamp, midnight, 
                                         units = "hours", tz = "Israel")) %>% 
    dplyr::select(-set_date, -rise_date, -sunset, -sunrise, -midnight)
  
  res <- st_as_sf(res, coords = c("x", "y"))
  st_crs(res) <- 4326
  res <- st_transform(res, crs = 32636) 
  xy <- st_coordinates(res)
  res <- res %>% 
    as.data.frame() %>% 
    cbind.data.frame(xy) %>% 
    dplyr::select(-geometry) %>% 
    rename(x = X, y = Y)
  # 
  # sp::coordinates(res) <- c("longitude", "latitude")
  # sp::proj4string(res) <- sp::CRS("+init=epsg:4326")
  # 
  # res <- sp::spTransform(res, sp::CRS("+init=epsg:32636")) %>% 
  #   as.data.frame() %>% 
  #   rename(x = longitude, y = latitude)
  
  return(res)
} 
