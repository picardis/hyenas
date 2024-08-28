test <- hcd[1:10, ]
test <- test %>% 
  mutate(timestamp = mdy_hm(paste(date, time))) %>% 
  select(id, timestamp, x = longitude, y = latitude)
saveRDS(test, "output/test_data.rds")
