# Multinomial regression of hyena behavioral state
# Simona Picardi
# May 25th, 2021

# Load packages ####

library(readxl)
library(tidyverse)
library(lubridate)
library(mclogit)
library(MASS)
library(patchwork)

# Load data ####

# Data intersected with land cover variables
dat <- read_excel("input/JoinHMM3StatesLandcovSettAgrivJNFFVPoultCowJNFCamp2105.xlsx")
# Data intersected with distance to roads
dat_roads <- read_excel("input/HMM_3States_NearDistance to Buffer100PavedRoad.xls")
# Data intersected with distance to settlement (use this as the main data!)
dat_sett <- read.csv("output/hmm_3states_time-from-midnight_2021-07-19_3projections_dist-to-settlements.csv")
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
# Read in HMM output
dat_hmm_new <- read.csv("output/hmm_3states_time-from-midnight_2021-07-19.csv") %>% 
  mutate(timestamp = ymd_hms(timestamp))

# Reclassify land cover ####

# Create new binary variables for livestock and poultry facilities
dat <- dat %>% 
  mutate(livestock = case_when(
    !is.na(bildingtyp) ~ "yes",
    is.na(bildingtyp) ~ "no"
  ),
  poultry = case_when(
    is.na(anafname) | anafname == "לול ריק" ~ "no",
    !is.na(anafname) ~ "yes"
  ))

# Add readable description for FCODE
dat <- dat %>% 
  mutate(fcode_descr = case_when(
    FCODE == 501 ~ "Artificial",
    FCODE == 502 ~ "Nat no veg",
    FCODE == 503 ~ "Ag Lawn Garden",
    FCODE == 504 ~ "Nat veg",
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
    FTYPE == 86 ~ "Med high density tall"
  ))

# Create column for JNF
dat <- dat %>% 
  mutate(jnf = case_when(
    !is.na(Types) ~ "yes",
    is.na(Types) ~ "no"
  ))

# Define new land cover variable
dat <- dat %>% 
  mutate(land_cover = case_when(
    ftype_descr == "Lawn or garden" ~ "Settlement",
    ftype_descr %in% c("Vineyards", "Fruit trees", "Olive trees", 
                       "Palm trees", "Cultivated", "Uncultivated") ~ "Agriculture",
    jnf == "yes" ~ "JNF",
    fcode_descr %in% c("Nat veg", "Nat no veg") ~ "Natural",
    TRUE ~ fcode_descr
  ))

# Create binary "settlement" column
dat <- dat %>% 
  mutate(settlement_yn = case_when(
    !is.na(settlement) ~ "yes",
    is.na(settlement) ~ "no"
  ))

# Add to land cover
dat <- dat %>% 
  mutate(land_cover = case_when(
    settlement_yn == "yes" ~ "Settlement",
    settlement_yn == "no" ~ land_cover
  ))

# We ended up not using GrowthCat because it's not as exhaustive as FTYPE 
# dat <- dat %>% 
#   mutate(ag_type = case_when(
#     GrowthCat %in% c("מטעים", "הדרים", "מטעים בכיסוי") ~ "Orchard",
#     GrowthCat %in% c("פרחים") ~ "Cultivated Flowers",
#     GrowthCat %in% c("לא מעובד") ~ "Uncultivated",
#     GrowthCat %in% c("גד\"ש וירקות בשטח פתוח" ,"ירקות בכיסוי") ~ "Crops and Veg",
#     is.na(GrowthCat) ~ NA_character_
#   ))

# Merge livestock and poultry into one
dat <- dat %>% 
  mutate(land_cover = case_when(
    livestock == "yes" ~ "Livestock and poultry facility",
    poultry == "yes" ~ "Livestock and poultry facility",
    TRUE ~ land_cover
  ))

# # Merge artificial and settlement
# dat <- dat %>% 
#   mutate(land_cover = case_when(
#     land_cover %in% c("Artificial", "Settlement") ~ "Anthropogenic",
#     TRUE ~ land_cover
#   ))

# Also use a simplified version of land cover
dat <- dat %>% 
  mutate(land_cover_simple = case_when(
    FTYPE %in% c(11, 60) ~ "Anthropogenic",
    FTYPE %in% c(31, 37, 78:86) ~ "Natural",
    FTYPE %in% c(51:57) ~ "Agriculture"
  ))

# Join new behavioral classification and distance variables
mn_dat <- dat %>% 
  dplyr::select(id = ID, timestamp = Time, state_old = State, land_cover, land_cover_simple) %>% 
  mutate(timestamp = mdy_hm(timestamp)) %>% 
  left_join(dat_hmm_new, by = c("id", "timestamp")) %>% 
  mutate(state = factor(case_when(
    state == 1 ~ "Resting",
    state == 2 ~ "Traveling",
    state == 3 ~ "Searching"
  ))) %>% 
  dplyr::select(-X) %>% 
  mutate(dist_roads = dat_roads$NEAR_DIST,
         dist_settlement = dat_sett$NEAR_DIST) %>% 
  left_join(info, by = "id") %>% 
  mutate(sex = factor(sex),
         id = factor(id))

# Multinomial regression ####

# Models 1 and 2 are not ideal because they have no random effect

mod1 <- mblogit(state ~ (land_cover + dist_settlement) * sex,
              data = mn_dat)

summary(mod1)

mod2 <- mblogit(state ~ (land_cover_simple + dist_settlement) * sex,
              data = mn_dat)

summary(mod2)

# Model 3 does not converge
mod3 <- mblogit(state ~ (land_cover + dist_settlement) * sex,
              data = mn_dat,
              random = ~ 1|id)

summary(mod3)

# Model 4 converges and it has random effect. This is our best option
mod4 <- mblogit(state ~ (land_cover_simple + dist_settlement) * sex,
              data = mn_dat,
              random = ~ 1|id)

summary(mod4)

AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)

# Following meeting on 06-30-2021, I'll make the following changes:
# 1. Use model with simple land cover;
# 2. Try to fit land cover as a series of three 0/1 boolean;
# 3. Include distance to roads again;
# 4. Log-transform distance variables;
# 5. Make all covariates interact;
# 6. Use individual-night as a random effect instead of individual.

mn_dat <- mn_dat %>% 
  mutate(agr = case_when(
    land_cover_simple == "Agriculture" ~ 1,
    TRUE ~ 0
  ),
  art = case_when(
    land_cover_simple == "Anthropogenic" ~ 1,
    TRUE ~ 0
  ),
  nat = case_when(
    land_cover_simple == "Natural" ~ 1,
    TRUE ~ 0
  ),
  log_dist_settlement = log(dist_settlement + 0.001),
  log_dist_roads = log(dist_roads + 0.001)) 
  
mod5 <-  mblogit(state ~ (agr + art + nat) * 
                   log_dist_settlement *
                   log_dist_roads *
                   sex,
                 data = mn_dat,
                 random = ~ 1|ind_night)

# Use individual-night nested within individual as random effect
mod6 <-  mblogit(state ~ (agr + art + nat) * 
                   log_dist_settlement *
                   log_dist_roads *
                   sex,
                 data = mn_dat,
                 random = ~ 1|(id/ind_night))

# Model predictions (mod 5) ####

log_dist_settlement <- log(rep(seq(0, 2500, length.out = 100) + 0.001, 6))
log_dist_roads <- log(rep(seq(0, 1000, length.out = 100) + 0.001, 6))
land_cover_simple <- rep(rep(c("Agriculture", "Anthropogenic", "Natural"), each = 100), 2)
agr <- ifelse(land_cover_simple == "Agriculture", 1, 0)
art <- ifelse(land_cover_simple == "Anthropogenic", 1, 0)
nat <- ifelse(land_cover_simple == "Natural", 1, 0)
sex <- rep(c("Female", "Male"), each = 300)
id <- "6163A"

newd_land <- data.frame(land_cover_simple = rep(c("Agriculture", "Anthropogenic", "Natural"), 2),
                        agr = c(1, 0, 0, 1, 0, 0), 
                        art = c(0, 1, 0, 0, 1, 0),
                        nat = c(0, 0, 1, 0, 0, 1),
                        log_dist_settlement = log(mean(mn_dat$dist_settlement, na.rm = T)),
                        log_dist_roads = log(mean(mn_dat$dist_roads, na.rm = T)),
                        sex = rep(c("Female", "Male"), each = 3),
                        id = id,
                        ind_night = 2)

newd_settlement <- data.frame(land_cover_simple = land_cover_simple,
                              agr = agr,
                              art = art,
                              nat = nat,
                              log_dist_settlement = log_dist_settlement,
                              log_dist_roads = log(mean(mn_dat$dist_roads, na.rm = T)),
                              sex = sex,
                              id = id,
                              ind_night = 2)

newd_roads <- data.frame(land_cover_simple = land_cover_simple,
                              agr = agr,
                              art = art,
                              nat = nat,
                              log_dist_roads = log_dist_roads, 
                              log_dist_settlement = log(mean(mn_dat$dist_settlement, na.rm = T)),
                              sex = sex,
                              id = id,
                              ind_night = 2)

preds_land <- predict(mod5, newdata = newd_land, conditional = FALSE, type = "response")
preds_settlement <- predict(mod5, newdata = newd_settlement, conditional = FALSE, type = "response")
preds_roads <- predict(mod5, newdata = newd_roads, conditional = FALSE, type = "response")

# Bootstrap confidence intervals (mod 5) ####

boot_land <- array(dim = c(6, 3, 1000))
boot_sett <- array(dim = c(600, 3, 1000))
boot_road <- array(dim = c(600, 3, 1000))

for (i in 1:1000) {
  
  print(i)
  
  samp <- mn_dat[sample(1:nrow(mn_dat), replace = TRUE), ]
  
  mod <- mblogit(state ~ (agr + art + nat) * 
                   log_dist_settlement *
                   log_dist_roads *
                   sex,
                  data = samp,
                  random = ~ 1|ind_night)
  
  boot_land[, , i] <- predict(mod, newdata = newd_land, conditional = FALSE, type = "response")
  boot_sett[, , i] <- predict(mod, newdata = newd_settlement, conditional = FALSE, type = "response")
  boot_road[, , i] <- predict(mod, newdata = newd_roads, conditional = FALSE, type = "response")
  
    }

means_land <- apply(boot_land, c(1, 2), FUN = mean)
means_sett <- apply(boot_sett, c(1, 2), FUN = mean)
means_road <- apply(boot_road, c(1, 2), FUN = mean)

lwr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.025)
lwr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.025)
lwr_road <- apply(boot_road, c(1, 2), FUN = quantile, 0.025)

upr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.975)
upr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.975)
upr_road <- apply(boot_road, c(1, 2), FUN = quantile, 0.975)

dist_lwr_lands <- means_land - lwr_land
dist_lwr_sett <- means_sett - lwr_sett
dist_lwr_road <- means_road - lwr_road

dist_upr_lands <- upr_land - means_land
dist_upr_sett <- upr_sett - means_sett
dist_upr_road <- upr_road - means_road

lwr_lands <- preds_land - dist_lwr_lands
colnames(lwr_lands) <- paste0(colnames(lwr_lands), "_lwr")
lwr_sett <- preds_settlement - dist_lwr_sett
colnames(lwr_sett) <- paste0(colnames(lwr_sett), "_lwr")
lwr_road <- preds_roads - dist_lwr_road
colnames(lwr_road) <- paste0(colnames(lwr_road), "_lwr")

upr_lands <- preds_land + dist_upr_lands
colnames(upr_lands) <- paste0(colnames(upr_lands), "_upr")
upr_sett <- preds_settlement + dist_upr_sett
colnames(upr_sett) <- paste0(colnames(upr_sett), "_upr")
upr_road <- preds_roads + dist_upr_road
colnames(upr_road) <- paste0(colnames(upr_road), "_upr")

# Plot (mod 5) ####

res_land <- cbind.data.frame(newd_land, preds_land) %>% 
  pivot_longer(cols = Resting:Traveling, names_to = "state", values_to = "pred")
res_settlement <- cbind.data.frame(newd_settlement, preds_settlement)  %>% 
  pivot_longer(cols = Resting:Traveling, names_to = "state", values_to = "pred")
res_roads <- cbind.data.frame(newd_roads, preds_roads)  %>% 
  pivot_longer(cols = Resting:Traveling, names_to = "state", values_to = "pred")

upr_lands <- upr_lands %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_upr:Traveling_upr, names_to = "state", values_to = "upr") %>% 
  dplyr::select(upr)
upr_sett <- upr_sett %>% 
  as_tibble %>% 
  pivot_longer(cols = Resting_upr:Traveling_upr, names_to = "state", values_to = "upr") %>% 
  dplyr::select(upr)
upr_road <- upr_road %>% 
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
lwr_road <- lwr_road %>% 
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
res_roads <- res_roads %>% 
  cbind.data.frame(upr_road) %>% 
  cbind.data.frame(lwr_road) %>% 
  mutate(dist_roads = exp(log_dist_roads) - 0.001)

saveRDS(res_land, "output/bootstrapped_multinomial_mod5_2021-07-20_response-to-land-cover.rds")
saveRDS(res_settlement, "output/bootstrapped_multinomial_mod5_2021-07-20_response-to-settlement.rds")
saveRDS(res_roads, "output/bootstrapped_multinomial_mod5_2021-07-20_response-to-roads.rds")

ggplot(res_land, aes(x = land_cover_simple, y = pred, color = state)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, width = 0.2)) +
  facet_wrap(~ sex) +
  labs(x = " ", y = "Probability", color = "State") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_mod5_2021-07-19_land-cover_sex.tiff",
       compression = "lzw", width = 6, height = 4, dpi = 300)

# Only plot curve within the actual range of distances found in the data for 
# each land cover type
mn_dat %>% 
  group_by(land_cover_simple) %>% 
  summarize(range(dist_settlement),
            range(dist_roads))
# Ranges are very similar. Don't bother.

ggplot(res_settlement, aes(x = dist_settlement, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(sex ~ land_cover_simple) +
  labs(x = "Distance to settlements (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_mod5_2021-07-19_distance-to-settlement.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300)

ggplot(res_roads, aes(x = dist_roads, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(sex ~ land_cover_simple) +
  labs(x = "Distance to roads (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_mod5_2021-07-19_distance-to-roads.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300)

# Keep x axis on log scale
ggplot(res_settlement, aes(x = log_dist_settlement, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(sex ~ land_cover_simple) +
  labs(x = "Distance to settlements (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom", axis.title.x.top = element_text(margin = margin(1.5))) +
  scale_x_continuous(breaks = log(c(0, 10, 100, 1000) + 0.001), 
                     labels = c(0, 10, 100, 1000)) 

ggsave("output/multinomial-model-predictions_mod5_2021-07-21_log_distance-to-settlement.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300)

ggplot(res_roads, aes(x = log_dist_roads, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(sex ~ land_cover_simple) +
  labs(x = "Distance to roads (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom", axis.title.x.top = element_text(margin = margin(1.5))) +
  scale_x_continuous(breaks = log(c(0, 10, 100, 1000) + 0.001), 
                     labels = c(0, 10, 100, 1000))

ggsave("output/multinomial-model-predictions_mod5_2021-07-21_log_distance-to-roads.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300)

# Model predictions (mod 1) ####

dist_settlement <- rep(seq(0, 2500, length.out = 100), 12)
land_cover <- rep(rep(unique(mn_dat$land_cover), each = 100), 2)
sex <- rep(c("Female", "Male"), each = 600)
id <- "6163A"

newd_land <- data.frame(land_cover = rep(unique(mn_dat$land_cover), 2),
                        dist_settlement = mean(mn_dat$dist_settlement, na.rm = T),
                        sex = rep(c("Female", "Male"), each = 6),
                        id = id)

newd_settlement <- data.frame(land_cover = land_cover,
                              dist_settlement = dist_settlement,
                              sex = sex,
                              id = id)

preds_land <- predict(mod1, newdata = newd_land, conditional = FALSE, type = "response")
preds_settlement <- predict(mod1, newdata = newd_settlement, conditional = FALSE, type = "response")

# Bootstrap confidence intervals (mod 1) ####

boot_land <- array(dim = c(12, 3, 1000))
boot_sett <- array(dim = c(1200, 3, 1000))

for (i in 1:1000) {
  
  print(i)
  
  try(
    
    {samp <- mn_dat[sample(1:nrow(mn_dat), replace = TRUE), ]
  
    mod <- mblogit(state ~ (land_cover + dist_settlement) * sex,
                 data = samp)
  
    boot_land[, , i] <- predict(mod, newdata = newd_land, conditional = FALSE, type = "response")
    boot_sett[, , i] <- predict(mod, newdata = newd_settlement, conditional = FALSE, type = "response")
  
  }) 
}

means_land <- apply(boot_land, c(1, 2), FUN = mean, na.rm = T)
means_sett <- apply(boot_sett, c(1, 2), FUN = mean, na.rm = T)

lwr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.025, na.rm = T)
lwr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.025, na.rm = T)

upr_land <- apply(boot_land, c(1, 2), FUN = quantile, 0.975, na.rm = T)
upr_sett <- apply(boot_sett, c(1, 2), FUN = quantile, 0.975, na.rm = T)

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

# Plot (mod 1) ####

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
  cbind.data.frame(lwr_sett)

res_land <- res_land %>% 
  mutate(land_cover = case_when(
    land_cover == "Livestock and poultry facility" ~ "Farming facility",
    TRUE ~ land_cover
  ))

ggplot(res_land, aes(x = land_cover, y = pred, color = state)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, width = 0.2), 
                position = position_dodge(0.3)) +
  facet_wrap(~ sex) +
  labs(x = " ", y = "Probability", color = "State") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_mod1_2021-06-02_land-cover_sex.tiff",
       compression = "lzw", width = 15, height = 5, dpi = 300)

res_settlement <- res_settlement %>% 
  mutate(land_cover = case_when(
    land_cover == "Livestock and poultry facility" ~ "Farming facility",
    TRUE ~ land_cover
  ))

ggplot(res_settlement, aes(x = dist_settlement, y = pred, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(sex ~ land_cover) +
  labs(x = "Distance to settlements (m)", y = "Probability", color = "State", fill = "State") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("output/multinomial-model-predictions_mod1_2021-06-02_distance-to-settlement.tiff",
       compression = "lzw", width = 12, height = 5, dpi = 300)

# Summary plots ####

# Proportion of steps in each land cover type for each sex 
# (broken down by individuals and then average?)

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
  summarize(avg_prop = mean(prop))

ggplot(steps_per_cat, aes(x = sex, y = avg_prop, fill = land_cover_simple)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = " ", y = "Proportion of steps", fill = "Land cover") +
  scale_fill_brewer(palette = "Set2")

ggsave("output/steps-by-land-cover.tiff",
       compression = "lzw", width = 7, height = 5, dpi = 300, scale = 0.7)
