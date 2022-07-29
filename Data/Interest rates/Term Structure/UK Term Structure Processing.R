########################################################################################
# Code to calculate forward force of mortality based on spot rates for the UK
########################################################################################

library(readxl)
library(tidyverse)

#-----------------------------------------------------
#Read data
#-----------------------------------------------------
Data_1970_2015 <- read_excel(path ="Data/Interest rates/Term Structure/GLC Nominal month end data_1970 to 2015.xlsx", 
                             sheet = "4. spot curve",
                             range = "A4:AY557"  )

Data_2016_2020 <- read_excel(path ="Data/Interest rates/Term Structure/GLC Nominal month end data_2016 to present.xlsx", 
                             sheet = "4. spot curve",
                             range = "A4:CC62"  )

Data_1970_2020 <- bind_rows(Data_1970_2015, Data_2016_2020)

Data_1970_2020_long <- Data_1970_2020 %>% 
  pivot_longer(cols = -`years:`, names_to = "term", values_to = "spot") %>% 
  rename(date = `years:`) %>% 
  mutate(term = as.numeric(term),
         spot = spot/100) %>% 
  filter(!is.na(date)) 
  

#-----------------------------------------------------
#Fill in missing values for the short end
#-----------------------------------------------------
Data_1970_2020_long_short_end <- Data_1970_2020_long %>% 
  filter(term < 5) %>% fill(spot, .direction = "up")
  
Data_1970_2020_long_long_end <- Data_1970_2020_long %>% 
  filter(term >= 5) 

Data_1970_2020_long <- bind_rows(Data_1970_2020_long_short_end, Data_1970_2020_long_long_end) %>% 
  arrange(date,term)

#-----------------------------------------------------
#Compute forward force of mortality:
# Given spot rate at term i (s_i) the forward force i
# is given by
#
#       ff_i = i log(1+s_i) - (i-1)log(1+s_{i-1})
#-----------------------------------------------------

#Use only integer terms
Data_1970_2020_long <- Data_1970_2020_long %>% filter(term %in% 1:40)

#Compute forward force
Data_1970_2020_long <- Data_1970_2020_long %>% 
  group_by(date) %>% 
  mutate(forward_force =  term * log(1+spot) - lag(term) * log(lag(1+spot))) %>% 
  ungroup() %>% 
  mutate(forward_force = ifelse(term == 1, log(1+spot), forward_force)) %>% 
  fill(forward_force) #Fill long end of the curve

write_csv(Data_1970_2020_long, path = "Data/Interest rates/Term Structure/Processed_Term_Structure.csv")

