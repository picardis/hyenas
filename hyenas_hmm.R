# # # # # # # # # # # # # # # # # # # # # # # # # 
# Behavioral segmentation of striped hyena tracks
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Sept 28, 2020 # # # # # # # # # # # # #
# Last updated Sept 28, 2020 # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# Load packages ####

library(tidyverse)
library(momentuHMM)
library(sp)
library(lubridate)
library(CircStats)

# Set seed ####

set.seed(7)

# Fit simple 3-state HMM (full 1-hour dataset) ####

hcd_1h <- readRDS("output/hcd_1hCTMM.rds")

hmm_data <- hcd_1h %>% 
  mutate(ID = id) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Fit a model with 3 behavioral states: 
# Resting, Restricted, and Traveling
# Resting will have very short step lengths and uniform turning angles
# Restricted will have intermediate step lengths and turning angles ~180
# Traveling will have long step lengths and turning angles ~0

hmm_null <- fitHMM(data = hmm_data,
               nbStates = 3,
               dist = list(step = "gamma", angle = "vm"),
               estAngleMean = list(angle = TRUE),
               Par0 = list(step = c(mean_1 = 10, 
                                    mean_2 = 250, 
                                    mean_3 = 2500,
                                    sd_1 = 10, 
                                    sd_2 = 100,
                                    sd_3 = 800,
                                    zeromass_1 = 0.8, 
                                    zeromass_2 = 0.1,
                                    zeromass_3 = 0.001),
                           angle = c(mean_1 = pi, 
                                     mean_2 = pi,
                                     mean_3 = 0,
                                     concentration_1 = 0.01,
                                     concentration_2 = 0.99,
                                     concentration_3 = 0.99)))

plot(hmm_null)

# Fit simple 3-state HMM (nightly dataset) ####

hcd_night <- readRDS("output/hcd_1hCTMM_night.rds")

hmm_data <- hcd_night %>% 
  mutate(ID = burst_) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

hmm_null <- fitHMM(data = hmm_data,
                   nbStates = 3,
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

plot(hmm_null)

hmm_null_states <- hmm_data %>% 
  dplyr::select(timestamp, id, burst_, x, y) %>%
  mutate(state = factor(viterbi(hmm_null)))

ggplot(hmm_null_states, aes(x = x, y = y, color = state, group = burst_)) +
  geom_path(alpha = 0.2) +
  geom_point(alpha = 0.2) +
  facet_wrap(~ id, scales = "free") +
  theme_bw()

# Fit 3-state HMM with time of night ####

hcd_night <- readRDS("output/hcd_1hCTMM_night.rds")

hmm_data <- hcd_night %>% 
  mutate(ID = burst_,
         hon = hour(timestamp),
         ton = factor(case_when(
           hour(timestamp) >= 18 & hour(timestamp) <= 21 ~ "dusk",
           hour(timestamp) >= 22 | hour(timestamp) <= 4 ~ "night",
           hour(timestamp) >= 5 & hour(timestamp) <= 8 ~ "dawn"
         ))) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Initial state probability does not depend on anything
# Transition probabilities are a function of time of night

form <- ~ ton

hmm_ton <- fitHMM(data = hmm_data,
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

plot(hmm_ton)

saveRDS(hmm_ton, "output/hmm_ton.rds")

# Plot transition probabilities ####

hmm_ton <- readRDS("output/hmm_ton.rds")

best <- hmm_ton$CIbeta$beta$est
bse <- hmm_ton$CIbeta$beta$se

tps <- data.frame(ton = rep(c("dusk", "night", "dawn"), each = 9),
           tp = rep(c("s1ts2", "s1ts3", "s2ts1", 
                      "s2ts3", "s3ts1", "s3ts2",
                      "s1ts1", "s2ts2", "s3ts3"), 3))

tps$mle <- NA
tps$lwr <- NA
tps$upr <- NA

# Dawn (intercept)

tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$mle <- best[1, 1]
tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$mle - 1.96 * bse[1, 1]
tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$mle + 1.96 * bse[1, 1]

tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$mle <- best[1, 2]
tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$mle - 1.96 * bse[1, 2]
tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$mle + 1.96 * bse[1, 2]

tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$mle <- best[1, 3]
tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$mle - 1.96 * bse[1, 3]
tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$mle + 1.96 * bse[1, 3]

tps[tps$ton == "dawn" & tps$tp == "s2ts3",]$mle <- best[1, 4]
tps[tps$ton == "dawn" & tps$tp == "s2ts3",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s2ts3",]$mle - 1.96 * bse[1, 4]
tps[tps$ton == "dawn" & tps$tp == "s2ts3",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s2ts3",]$mle + 1.96 * bse[1, 4]

tps[tps$ton == "dawn" & tps$tp == "s3ts1",]$mle <- best[1, 5]
tps[tps$ton == "dawn" & tps$tp == "s3ts1",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s3ts1",]$mle - 1.96 * bse[1, 5]
tps[tps$ton == "dawn" & tps$tp == "s3ts1",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s3ts1",]$mle + 1.96 * bse[1, 5]

tps[tps$ton == "dawn" & tps$tp == "s3ts2",]$mle <- best[1, 6]
tps[tps$ton == "dawn" & tps$tp == "s3ts2",]$lwr <- 
  tps[tps$ton == "dawn" & tps$tp == "s3ts2",]$mle - 1.96 * bse[1, 6]
tps[tps$ton == "dawn" & tps$tp == "s3ts2",]$upr <- 
  tps[tps$ton == "dawn" & tps$tp == "s3ts2",]$mle + 1.96 * bse[1, 6]

# Dusk

tps[tps$ton == "dusk" & tps$tp == "s1ts2",]$mle <- best[1, 1] + best[2, 1]
tps[tps$ton == "dusk" & tps$tp == "s1ts2",]$lwr <- tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s1ts2",]$mle - 1.96 * sqrt(bse[1, 1]^2 + bse[2, 1]^2)
tps[tps$ton == "dusk" & tps$tp == "s1ts2",]$upr <- tps[tps$ton == "dawn" & tps$tp == "s1ts2",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s1ts2",]$mle + 1.96 * sqrt(bse[1, 1]^2 + bse[2, 1]^2)

tps[tps$ton == "dusk" & tps$tp == "s1ts3",]$mle <- best[1, 2] + best[2, 2]
tps[tps$ton == "dusk" & tps$tp == "s1ts3",]$lwr <- tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s1ts3",]$mle - 1.96 * sqrt(bse[1, 2]^2 + bse[2, 2]^2)
tps[tps$ton == "dusk" & tps$tp == "s1ts3",]$upr <- tps[tps$ton == "dawn" & tps$tp == "s1ts3",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s1ts3",]$mle + 1.96 * sqrt(bse[1, 2]^2 + bse[2, 2]^2)

tps[tps$ton == "dusk" & tps$tp == "s2ts1",]$mle <- best[1, 2] + best[2, 2]
tps[tps$ton == "dusk" & tps$tp == "s2ts1",]$lwr <- tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s2ts1",]$mle - 1.96 * sqrt(bse[1, 2]^2 + bse[2, 2]^2)
tps[tps$ton == "dusk" & tps$tp == "s2ts1",]$upr <- tps[tps$ton == "dawn" & tps$tp == "s2ts1",]$mle +
  tps[tps$ton == "dusk" & tps$tp == "s2ts1",]$mle + 1.96 * sqrt(bse[1, 2]^2 + bse[2, 2]^2)

# LEFT OFF HERE, TRYING TO CALCULATE THE CIs ####

# MLE
tp <- data.frame(ton = c("dusk", "night", "dawn")) %>% 
    mutate(s1ts2 = case_when(
    ton == "dusk" ~ best[1, 1] + best[2, 1],
    ton == "night" ~ best[1, 1] + best[3, 1],
    TRUE ~ best[1, 1]),
    s1ts3 = case_when(
      ton == "dusk" ~ best[1, 2] + best[2, 2],
      ton == "night" ~ best[1, 2] + best[3, 2],
      TRUE ~ best[1, 2]),
    s2ts1 = case_when(
      ton == "dusk" ~ best[1, 3] + best[2, 3],
      ton == "night" ~ best[1, 3] + best[3, 3],
      TRUE ~ best[1, 3]),
    s2ts3 = case_when(
      ton == "dusk" ~ best[1, 4] + best[2, 4],
      ton == "night" ~ best[1, 4] + best[3, 4],
      TRUE ~ best[1, 4]),
    s3ts1 = case_when(
      ton == "dusk" ~ best[1, 5] + best[2, 5],
      ton == "night" ~ best[1, 5] + best[3, 5],
      TRUE ~ best[1, 5]),
    s3ts2 = case_when(
      ton == "dusk" ~ best[1, 6] + best[2, 6],
      ton == "night" ~ best[1, 6] + best[3, 6],
      TRUE ~ best[1, 6])) %>% 
  mutate_at(vars(s1ts2:s3ts2), plogis) %>% 
  mutate(s1ts1 = 1 - (s1ts2 + s1ts3),
         s2ts2 = 1 - (s2ts1 + s2ts3),
         s3ts3 = 1 - (s3ts1 + s3ts2)) %>% 
  pivot_longer(cols = s1ts2:s3ts3, names_to = "tp", values_to = "mle") 

# Lower CI
tpl <- data.frame(ton = c("dusk", "night", "dawn")) %>% 
  mutate(s1ts2 = case_when(
    ton == "dusk" ~ tp[19, 3] + tp[1, 3] - 1.96 * sqrt(bse[1, 1]^2 + bse[2, 1]^2),
    ton == "night" ~ tp[19, 3] + tp[10, 3] - 1.96 * sqrt(bse[1, 1]^2 + bse[3, 1]^2),
    TRUE ~ tp[19, 3] - 1.96 * bse[1, 1]),
    s1ts3 = case_when(
      ton == "dusk" ~ tp[20, 3] + tp[2, 3] - 1.96 * sqrt(bse[1, 2]^2 + bse[2, 2]^2),
      ton == "night" ~ tp[20, 3] + tp[11, 3] - 1.96 * sqrt(bse[1, 2]^2 + bse[3, 2]^2),
      TRUE ~ tp[20, 3] - 1.96 * bse[1, 2]),
    s2ts1 = case_when(
      ton == "dusk" ~ tp[21, 3] + tp[3, 3] - 1.96 * sqrt(bse[1, 3]^2 + bse[2, 3]^2),
      ton == "night" ~ tp[21, 3] + tp[12, 3] - 1.96 * sqrt(bse[1, 3]^2 + bse[3, 3]^2),
      TRUE ~ tp[21, 3] - 1.96 * bse[1, 3]),
    s2ts3 = case_when(
      ton == "dusk" ~ tp[22, 3] + tp[4, 3] - 1.96 * sqrt(bse[1, 4]^2 + bse[2, 4]^2),
      ton == "night" ~ tp[22, 3] + tp[13, 3] - 1.96 * sqrt(bse[1, 4]^2 + bse[3, 4]^2),
      TRUE ~ tp[22, 3] - 1.96 * bse[1, 4]),
    s3ts1 = case_when(
      ton == "dusk" ~ tp[23, 3] + tp[5, 3] - 1.96 * sqrt(bse[1, 5]^2 + bse[2, 5]^2),
      ton == "night" ~ tp[23, 3] + tp[14, 3] - 1.96 * sqrt(bse[1, 5]^2 + bse[3, 5]^2),
      TRUE ~ tp[23, 3] - 1.96 * bse[1, 5]),
    s3ts2 = case_when(
      ton == "dusk" ~ tp[24, 3] + tp[6, 3] - 1.96 * sqrt(bse[1, 6]^2 + bse[2, 6]^2),
      ton == "night" ~ tp[24, 3] + tp[15, 3] - 1.96 * sqrt(bse[1, 6]^2 + bse[3, 6]^2),
      TRUE ~ tp[24, 3] - 1.96 * bse[1, 6])) %>% 
  mutate_at(vars(s1ts2:s3ts2), plogis) 

# Upper CI
tpu <- data.frame(ton = c("dusk", "night", "dawn")) %>% 
  mutate(s1ts2 = case_when(
    ton == "dusk" ~ bupr[1, 1] + bupr[2, 1],
    ton == "night" ~ bupr[1, 1] + bupr[3, 1],
    TRUE ~ bupr[1, 1]),
    s1ts3 = case_when(
      ton == "dusk" ~ bupr[1, 2] + bupr[2, 2],
      ton == "night" ~ bupr[1, 2] + bupr[3, 2],
      TRUE ~ bupr[1, 2]),
    s2ts1 = case_when(
      ton == "dusk" ~ bupr[1, 3] + bupr[2, 3],
      ton == "night" ~ bupr[1, 3] + bupr[3, 3],
      TRUE ~ bupr[1, 3]),
    s2ts3 = case_when(
      ton == "dusk" ~ bupr[1, 4] + bupr[2, 4],
      ton == "night" ~ bupr[1, 4] + bupr[3, 4],
      TRUE ~ bupr[1, 4]),
    s3ts1 = case_when(
      ton == "dusk" ~ bupr[1, 5] + bupr[2, 5],
      ton == "night" ~ bupr[1, 5] + bupr[3, 5],
      TRUE ~ bupr[1, 5]),
    s3ts2 = case_when(
      ton == "dusk" ~ bupr[1, 6] + bupr[2, 6],
      ton == "night" ~ bupr[1, 6] + bupr[3, 6],
      TRUE ~ bupr[1, 6])) %>% 
  mutate_at(vars(s1ts2:s3ts2), plogis) %>% 
  mutate(s1ts1 = 1 - (s1ts2 + s1ts3),
         s2ts2 = 1 - (s2ts1 + s2ts3),
         s3ts3 = 1 - (s3ts1 + s3ts2)) %>% 
  pivot_longer(cols = s1ts2:s3ts3, names_to = "tp", values_to = "upr") 

# Join
tp <- tp %>% 
  left_join(tpl) %>% 
  left_join(tpu) %>% 
  mutate(state_from = word(tp, 1, 1, "t"),
         state_to = word(tp, 2, 2, "t")) %>% 
  mutate(state_to_descr = case_when(
    state_to == "s1" ~ "To Resting",
    state_to == "s2" ~ "To Searching",
    state_to == "s3" ~ "To Traveling"
  ),
  ton_descr = factor(case_when(
    ton == "dawn" ~ "Dawn",
    ton == "dusk" ~ "Dusk",
    ton == "night" ~ "Night"
  ), levels = c("Dusk", "Night", "Dawn"))) 

facet_labs <- c("s1" = "From Resting",
                "s2" = "From Searching",
                "s3" = "From Traveling")

ggplot(tp, aes(x = ton_descr, y = mle, color = state_to_descr)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
  facet_wrap(~ state_from, labeller = as_labeller(facet_labs)) + 
  theme_bw() + 
  xlab(" ") +
  ylab("Transition probability") +
  scale_color_discrete(name = " ") 

# Fit 3-state HMM with time from midnight ####

# Instead of splitting the night into three categorical phases,
# I will use the time since sunset/to sunrise. This way I can model 
# behavioral changes continuously and have a single parameter in the model. 

hcd_night <- readRDS("output/hcd_1hCTMM_night.rds")

source("calc_time_from_midnight.R")

hmm_data <- hcd_night %>% 
  calc_tfm() %>%
  mutate(ID = burst_,
         tfm = as.numeric(time_from_midnight)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

form <- ~ tfm

hmm_tfm <- fitHMM(data = hmm_data,
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

plot(hmm_tfm)

saveRDS(hmm_tfm, "output/hmm_tfm.rds")

# Model predictions - transition probabilities ####

source("predict_hmm.R")

pred_tfm <- data.frame(tfm = seq(-7, 10, length.out = 100))

preds <- predict_hmm(mod = hmm_tfm, newdata = pred_tfm) %>% 
  mutate(state_from = case_when(
    word(tp, 1, 1, "t") == 1 ~ "From Resting",
    word(tp, 1, 1, "t") == 2 ~ "From Searching",
    word(tp, 1, 1, "t") == 3 ~ "From Traveling"
  ),
  state_to = case_when(
    word(tp, 2, 2, "t") == 1 ~ "To Resting",
    word(tp, 2, 2, "t") == 2 ~ "To Searching",
    word(tp, 2, 2, "t") == 3 ~ "To Traveling"
  ))

ggplot(preds, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

ggsave("output/hmm_tfm_tp.tiff", compression= "lzw", 
       scale = 0.8, height = 5, width = 12, dpi = 400)

# Model predictions - initial state probabilities ####

isp <- data.frame(state = c("Resting", "Searching", "Traveling"),
                  est = hmm_tfm$CIreal$delta$est[1, ],
                  lwr = hmm_tfm$CIreal$delta$lower[1, ],
                  upr = hmm_tfm$CIreal$delta$upper[1, ]) 

ggplot(isp, aes(x = state, y = est, color = state)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
  theme_bw() +
  xlab(" ") +
  ylab("Initial state probability") +
  scale_color_discrete(name = " ") 

ggsave("output/hmm_tfm_isp.tiff", compression= "lzw", 
       scale = 0.8, height = 5, width = 5, dpi = 400)

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
theta1 <- hmm_tfm$mle$step[1,1]/hmm_tfm$mle$step[2,1]
k1 <- hmm_tfm$mle$step[1,1]/theta1

# State 2
theta2 <- hmm_tfm$mle$step[1,2]/hmm_tfm$mle$step[2,2]
k2 <- hmm_tfm$mle$step[1,2]/theta2

# State 3
theta3 <- hmm_tfm$mle$step[1,3]/hmm_tfm$mle$step[2,3]
k3 <- hmm_tfm$mle$step[1,3]/theta3

# Simulate and plot
sim_step <- data.frame(sim = c(rgamma(n = 10000, shape = k1, scale = theta1),
                               rgamma(n = 10000, shape = k2, scale = theta2),
                               rgamma(n = 10000, shape = k3, scale = theta3)),
                       state = rep(c("Resting", "Searching", "Traveling"), 
                                   each = 10000)) %>% 
  mutate(state = factor(state, levels = c("Resting", "Searching", "Traveling")))

ggplot(sim_step, aes(x = sim, fill = state)) +
  geom_histogram(binwidth = 20, color = "black") +
  facet_wrap(~ state) +
  labs(x = "Step length (m)", y = "Frequency") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/hmm_tfm_sl.tiff", compression= "lzw", 
       scale = 0.8, height = 5, width = 12, dpi = 400)

# Turning angles

# Simulate and plot
sim_angle <- data.frame(sim = c(rvm(n = 10000, mean = hmm_tfm$mle$angle[1,1], 
                                    k = hmm_tfm$mle$angle[2,1]),
                                rvm(n = 10000, mean = hmm_tfm$mle$angle[1,2], 
                                    k = hmm_tfm$mle$angle[2,2]),
                                rvm(n = 10000, mean = hmm_tfm$mle$angle[1,3], 
                                    k = hmm_tfm$mle$angle[2,3])),
                        state = rep(c("Resting", "Searching", "Traveling"), 
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

ggsave("output/hmm_tfm_ta.tiff", compression= "lzw", 
       scale = 1.2, height = 5, width = 10, dpi = 400)

# Fit 3-state HMM with time from midnight and individual random effects ####

hmm_tfm_rind <- fitHMM(data = hmm_data,
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
                                             concentration_3 = 0.99)),
                       mixtures = 8,
                       formulaPi = ~ id)

saveRDS(hmm_tfm_rind, "output/hmm_8mixtures_random_id.rds")

vit_null <- viterbi(hmm_tfm)
vit_rf <- viterbi(hmm_tfm_rind)
table(vit_null == vit_rf)
