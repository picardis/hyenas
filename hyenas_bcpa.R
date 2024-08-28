library(bcpa)
library(tidyverse)
library(lubridate)

hyenas <- readRDS("output/hyenas_data_cleaned.rds")

hyenas_night <- hyenas %>% 
  mutate(timestamp = round.POSIXt(timestamp, "hours")) %>% 
  filter(hour(timestamp) > 17 | hour(timestamp) < 9) %>% 
  group_by(id) %>% 
  mutate(lag = as.numeric((timestamp) - lag(timestamp))/3600) %>% 
  mutate(lag = case_when(
    is.na(lag) ~ 1,
    TRUE ~ lag
  )) %>% 
  mutate(burst_ = NA_real_)

hyenas_night$burst_[1] <- 1

for (i in 2:nrow(hyenas_night)) {
  if (hyenas_night$id[i-1] != hyenas_night$id[i] | hyenas_night$lag[i] != 1) {
    hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] + 1
  } else if (hyenas_night$lag[i] == 1) {
    hyenas_night$burst_[i] <- hyenas_night$burst_[i-1] 
  }
}

expected_dat <- hyenas_night %>% 
  group_by(burst_, id) %>% 
  summarize(from = as_date(min(timestamp)), to = as_date(max(timestamp)))

hours_range <- c(18:23, 0:8)

expected_ts <- data.frame()

for (i in 1:length(unique(expected_dat$burst_))) {
  
  sub <- expected_dat[expected_dat$burst_ == unique(expected_dat$burst_)[i],]
  who <- unique(sub$id)
  ts <- c(ymd_hms(paste0(sub$from, " ", hours_range[1:6], ":00:00"), 
                  tz = "Israel"),
          ymd_hms(paste0(sub$from + days(1), " ", hours_range[7:15], ":00:00"),
                  tz = "Israel"))
  res <- cbind.data.frame(burst_ = unique(expected_dat$burst_)[i],
                          id = who,
                          timestamp = ts)
  
  expected_ts <- rbind(expected_ts, res)
  
}

hyenas_night_exp <- left_join(expected_ts, hyenas_night, 
                          by = c("burst_", "timestamp", "id"))

# Let's say we want to keep only nights that have at least ~50% of data
keep <- hyenas_night_exp %>% 
  filter(!is.na(x)) %>% 
  group_by(burst_) %>% 
  tally() %>% 
  arrange(n) %>% 
  filter(n >= 7) %>% 
  pull(burst_)

hyenas_night_exp <- hyenas_night_exp %>% 
  filter(burst_ %in% keep)

test <- hyenas %>% 
  filter(id == unique(hyenas$id)[1])

test <- hyenas_night_exp %>% 
  filter(burst_ == unique(hyenas_night_exp$burst_)[6])

track <- MakeTrack(X = test$x, 
                   Y = test$y, 
                   Time = test$timestamp)

vt <- GetVT(track)

GetRho(vt$V, vt$T.start, tau = TRUE)

bb <- GetBestBreak(vt$V, vt$T.start, tau = FALSE)

GetModels(vt$V, vt$T.start, bb[1], tau = FALSE)

simp.ws <- WindowSweep(vt, "V*cos(Theta)", windowsize = 4, K = 2)

plot(simp.ws, type = "smooth")
plot(simp.ws, type = "flat")

ChangePointSummary(simp.ws)

PathPlot(track, simp.ws, type = "flat")
PhasePlot(simp.ws, type = "flat")

DiagPlot(simp.ws)
