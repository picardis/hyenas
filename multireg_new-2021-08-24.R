# Multinomial regression of hyena behavioral state
# Simona Picardi
# August 24th, 2021

# Load packages ####

library(readxl)
library(tidyverse)
library(lubridate)
library(mclogit)
library(MASS)
library(patchwork)

# Load data ####

# Data intersected with distance to settlement (use this as the main data!)
dat_sett <- read.csv("output/hmm_3states_time-from-midnight_2021-07-19_3projections_dist-to-settlements.csv") %>% 
  drop_na(itm_x, itm_y)

# Data intersected with distance to roads 
dat_roads <- read.csv("input/2021-08/HMM3States_time-from-midnight_190721_Disttoroads.csv") %>% 
  mutate(date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp, -FID, -X, -NEAR_FID) %>% 
  rename(dist_roads = NEAR_DIST)

# Data intersected with agricultural areas
dat_ag <- read_xls("input/2021-08/Agriculture.xls") %>% 
  mutate(in_ag = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with JNFF
dat_jnff <- read_xls("input/2021-08/JNFF.xls") %>% 
  # filter out points that are in olive trees or fruit trees
  filter(!COV_TYPE %in% c(2994, 2990)) %>% 
  mutate(in_jnff = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with JNFV
dat_jnfv <- read_xls("input/2021-08/JNFV.xls") %>% 
  mutate(in_jnfv = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with JNF campsites
dat_jnfc <- read_xls("input/2021-08/JNFCampsites_Buffer50.xls") %>% 
  mutate(in_jnfc = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with land cover
dat_lc <- read_xls("input/2021-08/Landcover.xls") %>% 
  mutate(date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with livestock facilities
dat_live <- read_xls("input/2021-08/Livestock_Buffer50.xls") %>% 
  mutate(livestock = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

# Data intersected with poultry facilities
dat_poul <- read_xls("input/2021-08/Poultry_Buffer50.xls") %>% 
  mutate(poultry = 1,
         date = as_date(timestamp)) %>% 
  dplyr::select(-timestamp)

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
  left_join(dat_roads, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  left_join(dat_ag, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag) %>% 
  left_join(dat_jnfc, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc) %>% 
  left_join(dat_jnff, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff) %>% 
  left_join(dat_jnfv, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv) %>% 
  left_join(dat_live, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv, livestock) %>% 
  left_join(dat_poul, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
  dplyr::select(ind_night, id, timestamp, itm_x, itm_y, state, within_sett, dist_to_sett, date, dist_roads,
                in_ag, in_jnfc, in_jnff, in_jnfv, livestock, poultry) %>% 
  left_join(dat_lc, by = c("ind_night", "id", "itm_x", "itm_y", "state", "date")) %>% 
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

write.csv(dat, "output/multinomial-regression-data_2021-08-30.csv")

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
