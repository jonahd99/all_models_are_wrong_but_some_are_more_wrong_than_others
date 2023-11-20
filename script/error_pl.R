################################################################################
# Read me                                                                      #
################################################################################



################################################################################
# Preamble                                                                     #
################################################################################

library("here")
library("tidyverse")


################################################################################
# Blythe and Kiraly (2016) -- Running dataset                                  #
################################################################################

df_power_of_10 <- read_rds(here("data", "df_power_of_10.Rds")) %>%
  select(!type)

################################################################################
# IDs of runners that have completed at least 1 marathon                       #
################################################################################


df_power_of_10 %>%
  filter(distance == 42194.988) %>%
  pull(id) %>%
  unique() -> ID_run_marathon


################################################################################
# IDs of runners that initially enough training points                         #
################################################################################

df_power_of_10 %>%
  filter(id %in% ID_run_marathon & distance != 42194.988 & distance >= 400) %>%
  group_by(id) %>%
  summarise(n_dis = n_distinct(distance),
            n_train = n()) %>%
  filter(n_dis >= 3 & n_train >= 3) %>%
  pull(id) %>%
  unique() -> ID_enough_train_init

################################################################################
# Data frame containing the number of marathons each runner has completed      #
################################################################################

df_power_of_10 %>%
  filter(id %in% ID_enough_train_init) %>%
  filter(distance == 42194.988) %>% 
  group_by(id) %>%
  summarise(n_marathons = n()) -> df_n_marathons

################################################################################
# Data frame that only contains marathons -- i.e. Testing dataset              #
################################################################################

df_power_of_10 %>%
  filter(distance == 42194.988 & id %in% ID_enough_train_init) %>%
  mutate(date_lb = date %m-% months(4)) %>% 
  group_by(id) %>%
  mutate(marathon_n = row_number(),
         updated_id = paste0(id,"_",marathon_n)) %>%
  ungroup() %>%
  select(id,updated_id,marathon_n,date,date_lb,distance,duration,power) %>%
  rename(date_up = date) -> marathons_only

# data_lb = The date 4 months prior to the testing marathon 
# data_up = The date of the testing marathon 
# Both of these will be used to filter out data points that 
# will not be included in the training dataset.


###############################################################################
# Data frame that only contains observations from athletes with enough intial #
# training data                                                               #              
###############################################################################

df_power_of_10 %>%
  filter(id %in% ID_enough_train_init & distance >= 400)  -> training_data_init

###############################################################################
# For each athlete, create a dataset that contains duplicated training        #
# observations depending on the number of marathons completed by the athlete. #
# Each marathon will be treated essentially as it's own athlete. And the      #
# PL model will be fit to data from up to 4 months prior.                     #           
###############################################################################



duplicated_training_data <- tibble()

re_wrangle <- FALSE

if(re_wrangle == TRUE){
  for(i in 1:length(ID_enough_train_init)){
    
    ID <- ID_enough_train_init[i]
    
    n_marathons <- df_n_marathons %>%
      filter(id == ID) %>%
      pull(n_marathons)
    
    df_i <- training_data_init %>%
      filter(id == ID)
    
    aux_marathon_n <- rep(1:n_marathons, each = nrow(df_i))
    
    duplicated_training_data_aux <- df_i %>%
      slice(rep(1:n(), times = n_marathons)) %>%
      cbind(marathon_n = aux_marathon_n) %>%
      mutate(updated_id = paste0(id, "_", marathon_n))
    duplicated_training_data <- rbind(duplicated_training_data,duplicated_training_data_aux)
  }
  
  saveRDS(duplicated_training_data,here("data","duplicated_training_data.RDS"))
}else{
  duplicated_training_data <- readRDS(here("data","duplicated_training_data.RDS"))
}

###############################################################################
# Getting the IDs for the athletes that now have enough training data within  #        
# the 4 months prior to the testing marathon                                  #
###############################################################################

marathons_only %>%
  select(updated_id,date_lb,date_up) %>%
  left_join(duplicated_training_data, by = c("updated_id"), multiple = "all") %>%
  mutate(is_in = if_else(date >= date_lb & date < date_up,1,0)) %>%
  filter(is_in == 1 & distance >= 400) %>%
  group_by(updated_id) %>%
  summarise(n_dis = n_distinct(distance),
            n_train = n()) %>%
  filter(n_dis >= 3 & n_train >= 3) %>%
  pull(updated_id) %>%
  unique() -> ID_final


###############################################################################
# Creating the final training dataset                                         #
###############################################################################

marathons_only %>%
  select(updated_id,date_lb,date_up, duration, power, marathon_n) %>%
  rename(act_duration = duration,
         act_velocity = power) %>%
  left_join(duplicated_training_data, by = c("updated_id","marathon_n"), multiple = "all") %>%
  mutate(is_in = if_else(date >= date_lb & date < date_up,1,0)) %>%
  filter(is_in == 1 & updated_id %in% ID_final) -> training_data

###############################################################################
# Get the Endurance parameter for the PL model                                #
###############################################################################

get_E <- function(power,durations){
  data_pl <- tibble(power = power, durations = durations)
  pl_lr <- lm(data = data_pl, I(log(power)) ~ I(log(durations)))
  b <- 1 - as.numeric(-coef(pl_lr)[2])
  return(b)
}

###############################################################################
# Get the Sprint parameter for the PL model                                   #
###############################################################################

get_S <- function(power,durations){
  data_pl <- tibble(power = power, durations = durations)
  pl_lr <- lm(data = data_pl, I(log(power)) ~ I(log(durations)))
  a <- as.numeric(exp(coef(pl_lr)[1]))
  return(a)
}

###############################################################################
# Creating a data frame of the actual finishing time for each marathon        #
###############################################################################

marathons_only %>%
  select(updated_id,duration) %>%
  rename(mar_fin_time = duration) -> marathon_fin_time_df

###############################################################################
# Calculating the error for each marathon                                     #
###############################################################################

training_data %>%
  group_by(updated_id) %>%
  mutate(weighting = 1/as.numeric(date_up - date)) %>%
  summarise(E = get_E(power = power, durations = duration),
            S = get_S(power = power, durations = duration)) %>%
  mutate(mar_fin = (42194.988/S)^(1/(E))) %>%
  left_join(marathon_fin_time_df, by = "updated_id") %>%
  mutate(error = abs(mar_fin - mar_fin_time)/ mar_fin_time) -> error_df




###############################################################################
# Average error for N =  1697                                                 #
###############################################################################

mean(error_df$error, na.rm = TRUE)

