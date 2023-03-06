# pacman::p_load(tidyverse, lubridate, here)
# import.data <- TRUE

# Catch data
## Import catch data
if (import.data) {
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/catch_final.Rdata")
  saveRDS(catch, here::here("Data/trawl/catch_1907RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/catch_final.Rdata")
  saveRDS(catch, here::here("Data/trawl/catch_2107RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/catch_final.Rdata")
  saveRDS(catch, here::here("Data/trawl/catch_2207RL.rds"))
} 

# Load catch data
catch.2019 <- readRDS(here::here("Data/trawl/catch_1907RL.rds")) %>% 
  mutate(survey = "1907RL")
catch.2021 <- readRDS(here::here("Data/trawl/catch_2107RL.rds")) %>% 
  mutate(survey = "2107RL")
catch.2022 <- readRDS(here::here("Data/trawl/catch_2207RL.rds")) %>% 
  mutate(survey = "2207RL")

# Combine catch data
catch.trawl <- catch.2019 %>% bind_rows(catch.2021) %>% bind_rows(catch.2022)
saveRDS(catch.trawl, here::here("Data/trawl/catch_trawl.rds"))

# Specimen data
## Import specimen data
if (import.data) {
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_1907RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_2107RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_2207RL.rds"))
} 

# Load specimen data
lengths.2019 <- readRDS(here::here("Data/trawl/lengths_1907RL.rds")) %>% 
  mutate(survey = "1907RL")
lengths.2021 <- readRDS(here::here("Data/trawl/lengths_2107RL.rds")) %>% 
  mutate(survey = "2107RL")
lengths.2022 <- readRDS(here::here("Data/trawl/lengths_2207RL.rds")) %>% 
  mutate(survey = "2207RL")

# Combine specimen data
lengths.trawl <- lengths.2019 %>% bind_rows(lengths.2021) %>% bind_rows(lengths.2022)
saveRDS(lengths.trawl, here::here("Data/trawl/lengths_trawl.rds"))
