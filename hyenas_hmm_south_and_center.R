# # # # # # # # # # # # # # # # # # # # # # # # # 
# Behavioral segmentation of striped hyena tracks
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Code by Simona Picardi # # # # # # # # # # # # 
# Created Nov 05, 2020 # # # # # # # # # # # # #
# Last updated Oct 05, 2020 # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # 

# Load packages ####

library(tidyverse)
library(momentuHMM)
library(sp)
library(lubridate)
library(CircStats)

# Set seed ####

set.seed(7)

# Load data ####

hyenas_night <- readRDS("output/full_1hCTMM_night.rds")

# Prepare individual info ####

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
  dplyr::select(id, district = location, sex_rs, legs)
# I left age out because the age is NA for the hyenas in the South District 
# and all the others are adult except for one young adult. 

# In the data, the IDs of the South District hyenas follow the pattern 
# sex-number, in the info file they follow number-sex. Change these manually:
info$id[9:13] <- c("Female16", "Female32", "Female40", "Male32", "Male33")

# Fit full 3-state HMM ####

# This model will include as covariates:
# - District;
# - Sex/reproductive status;
# - Number of legs;
# - Time from midnight.

source("calc_time_from_midnight.R")

hmm_data <- hyenas_night %>% 
  calc_tfm() %>%
  mutate(ID = burst_,
         tfm = as.numeric(time_from_midnight)) %>% 
  left_join(info, by = "id") %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  mutate(sex_rs = factor(sex_rs),
         district = factor(district),
         legs = factor(legs)) %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

form <- ~ tfm + district + sex_rs

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

# I saved different objects here after changing the formula, see list of models 
# below.
saveRDS(hmm_tfm, "output/hmm_tfm_district_sexrs.rds")

# Model with full dataset, tfm only: works fine
# Model works but estimates very small or negative variances.

# Model with full dataset, tfm + sex_rs: works fine
# Model works but estimates very small or negative variances.

# Model with full dataset, tfm + district: works fine
# Model works but estimates very small or negative variances.

# Model with full dataset, tfm + sex_rs + district: runs, but then I get this 
# error when plotting. Model works but estimates very small or negative variances.
# Error in n2w(tmpPar, tmpp$bounds, beta, delta, nbStates, tmpInputs$estAngleMean,  : 
#                Check the parameter bounds for step (the initial parameters should be strictly between the bounds of their parameter space).
             
# Model with full dataset, tfm + sex_rs + district + legs: fails to converge
# The individual with 3 legs has 1516 locations, it's the third smallest 
# dataset in the Central District (5th smallest overall)
# Error in nlm(nLogLike, optPar, nbStates, newformula, p$bounds, p$parSize,  : 
#                non-finite value supplied by 'nlm'

# Model predictions: transition probabilities ####

hmm_tfm <- readRDS("output/hmm_tfm_district_sexrs.rds")

source("predict_hmm.R")

# Predict using each possible combination of reference factor levels 

# Center, F

pred_cf <- data.frame(tfm = seq(-7, 10, length.out = 100),
                       district = factor("Center", 
                                         levels = c("Center", "South")),
                       sex_rs = factor("F", levels = c("M", "F", "RF")))

preds_cf <- predict_hmm(mod = hmm_tfm, newdata = pred_cf) %>% 
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

plot_cf <- ggplot(preds_cf, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  labs(title = "Central + F") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

# Center, M

pred_cm <- data.frame(tfm = seq(-7, 10, length.out = 100),
                      district = factor("Center", 
                                        levels = c("Center", "South")),
                      sex_rs = factor("M", levels = c("M", "F", "RF")))

preds_cm <- predict_hmm(mod = hmm_tfm, newdata = pred_cm) %>% 
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

plot_cm <- ggplot(preds_cm, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  labs(title = "Central + M") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

# Center, RF

pred_crf <- data.frame(tfm = seq(-7, 10, length.out = 100),
                      district = factor("Center", 
                                        levels = c("Center", "South")),
                      sex_rs = factor("RF", levels = c("M", "F", "RF")))

preds_crf <- predict_hmm(mod = hmm_tfm, newdata = pred_crf) %>% 
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

plot_crf <- ggplot(preds_crf, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  labs(title = "Central + RF") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

# South, F

pred_sf <- data.frame(tfm = seq(-7, 10, length.out = 100),
                      district = factor("South", 
                                        levels = c("Center", "South")),
                      sex_rs = factor("F", levels = c("M", "F", "RF")))

preds_sf <- predict_hmm(mod = hmm_tfm, newdata = pred_sf) %>% 
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

plot_sf <- ggplot(preds_sf, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  labs(title = "South + F") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

# South, M

pred_sm <- data.frame(tfm = seq(-7, 10, length.out = 100),
                      district = factor("South", 
                                        levels = c("Center", "South")),
                      sex_rs = factor("M", levels = c("M", "F", "RF")))

preds_sm <- predict_hmm(mod = hmm_tfm, newdata = pred_sm) %>% 
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

plot_sm <- ggplot(preds_sm, aes(x = tfm, y = est, color = state_to, fill = state_to)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  facet_wrap(~ state_from) +
  theme_bw() +
  xlab("Hours from solar midnight") +
  ylab("Transition probability") +
  labs(title = "South + M") +
  scale_color_discrete(name = " ") +
  scale_fill_discrete(name = " ")

# Make multi-panel plot

ggpubr::ggarrange(plot_cf, plot_cm, plot_crf,
                  plot_sf, plot_sm, nrow = 2, ncol = 3,
                  common.legend = TRUE, legend = "bottom")

ggsave("output/hmm_tfm_district_sexrf_tp.tiff", compression= "lzw", 
       height = 5, width = 12, dpi = 400)

# Separate models for the two districts ####

# Central District

hmm_data_center <- hmm_data %>% 
  filter(district == "Center")

form <- ~ tfm + sex_rs

hmm_tfm <- fitHMM(data = hmm_data_center,
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

saveRDS(hmm_tfm, "output/hmm_center_tfm_sexrs.rds")

# South District

hmm_data_south <- hmm_data %>% 
  filter(district == "South") %>% 
  mutate(sex = factor(sex_rs, levels = c("M", "F")))

# There are no reproductive females in the South
form <- ~ tfm + sex

hmm_tfm <- fitHMM(data = hmm_data_south,
                  nbStates = 3,
                  formula = form,
                  dist = list(step = "gamma", angle = "vm"),
                  estAngleMean = list(angle = TRUE),
                  Par0 = list(step = c(mean_1 = 100, 
                                       mean_2 = 500,
                                       mean_3 = 4000, 
                                       sd_1 = 50, 
                                       sd_2 = 500,
                                       sd_3 = 2000),
                              angle = c(mean_1 = pi,
                                        mean_2 = pi,
                                        mean_3 = 0,
                                        concentration_1 = 0.1,
                                        concentration_2 = 0.5,
                                        concentration_3 = 0.99)))

# Error in nlm(nLogLike, optPar, nbStates, newformula, p$bounds, p$parSize,  : 
#                non-finite value supplied by 'nlm'

# Maybe not enough data for males?
hmm_data_south %>% group_by(sex) %>% tally()

# Predictions for model with tfm + sexrf (no district) ####

hmm_tfm <- readRDS("output/hmm_tfm_sexrs.rds")

pred <- data.frame(tfm = seq(-7, 10, length.out = 100),
                   sex_rs = factor("F", levels = c("M", "F", "RF")))

preds <- predict_hmm(mod = hmm_tfm, newdata = pred) %>% 
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
