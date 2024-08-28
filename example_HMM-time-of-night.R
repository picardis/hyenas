# Recreate a snippet of the data
hcd_night <- structure(list(timestamp = structure(c(1526918280, 1526921880, 
                                                    1526925480, 1526929080, 1526932680, 1526936280, 1526939880, 1526943480, 
                                                    1526947080, 1526950680, 1526954280, 1526957880, 1526961480, 1526965080, 
                                                    1526968680, 1527004680, 1527008280, 1527011880, 1527015480, 1527019080
), tzone = "Israel", class = c("POSIXct", "POSIXt")), id = c("6163A", 
                                                             "6163A", "6163A", "6163A", "6163A", "6163A", "6163A", "6163A", 
                                                             "6163A", "6163A", "6163A", "6163A", "6163A", "6163A", "6163A", 
                                                             "6163A", "6163A", "6163A", "6163A", "6163A"), x = c(693497.880154616, 
                                                                                                                 693398.229696095, 692718.988446969, 692738.795276243, 692714.246837338, 
                                                                                                                 692688.298803264, 692665.377454328, 692661.255382117, 692658.199156029, 
                                                                                                                 692662.843318162, 692651.092169263, 693347.407193595, 693364.225533052, 
                                                                                                                 693373.63013825, 693390.524689462, 693352.361889505, 693207.626417661, 
                                                                                                                 692723.268380071, 692650.833047312, 692707.366898032), y = c(3525471.97739037, 
                                                                                                                                                                              3525472.89402112, 3525940.62836049, 3526002.88619798, 3525992.54600388, 
                                                                                                                                                                              3525996.34898064, 3525989.45325536, 3525986.97119702, 3525989.12021905, 
                                                                                                                                                                              3525990.58914534, 3525998.46589508, 3525347.69448305, 3525350.09966774, 
                                                                                                                                                                              3525338.50761993, 3525317.51996006, 3525342.72127728, 3525519.84584497, 
                                                                                                                                                                              3526008.13164111, 3526026.6462346, 3525986.46665157), burst_ = c(2, 
                                                                                                                                                                                                                                               2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3)), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                        -20L), class = c("tbl_df", "tbl", "data.frame"))

# Source function to calculate time from midnight based on timestamp and 
# geographic location
source("calc_time_from_midnight.R")

# Prep the data for HMM fitting
hmm_data <- hcd_night %>% 
  calc_tfm() %>%
  mutate(ID = burst_,
         tfm = as.numeric(time_from_midnight)) %>% 
  arrange(ID, timestamp) %>% 
  as.data.frame() %>% 
  prepData(type = "UTM",
           coordNames = c("x", "y"))

# Formula for transition probabilities
form <- ~ tfm

# Fit model
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