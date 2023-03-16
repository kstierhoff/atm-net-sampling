pacman::p_load(tidyverse, lubridate, here, scatterpie, patchwork, sf)
pacman::p_load_gh("kstierhoff/atmData")
pacman::p_load_gh("kstierhoff/atm")

import.data <- TRUE

source(here::here("Doc/settings.R"))

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

# Specimen data----------------------------------------------
## Import specimen data
if (import.data) {
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_1907RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_2107RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/lengths_final.Rdata")
  saveRDS(lengths, here::here("Data/trawl/lengths_2207RL.rds"))
} 

## Load specimen data
lengths.2019 <- readRDS(here::here("Data/trawl/lengths_1907RL.rds")) %>% 
  mutate(survey = "1907RL")
lengths.2021 <- readRDS(here::here("Data/trawl/lengths_2107RL.rds")) %>% 
  mutate(survey = "2107RL")
lengths.2022 <- readRDS(here::here("Data/trawl/lengths_2207RL.rds")) %>% 
  mutate(survey = "2207RL")  

# Combine specimen data
lengths.trawl <- lengths.2019 %>% bind_rows(lengths.2021) %>% bind_rows(lengths.2022)
saveRDS(lengths.trawl, here::here("Data/trawl/lengths_trawl.rds"))

## Import clf data
if (import.data) {
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/clf_ts_proportions.Rdata")
  saveRDS(clf, here::here("Data/trawl/clf_1907RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/clf_ts_proportions.Rdata")
  saveRDS(clf, here::here("Data/trawl/clf_2107RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/clf_ts_proportions.Rdata")
  saveRDS(clf, here::here("Data/trawl/clf_2207RL.rds"))
}

# Load clf data
clf.2019 <- readRDS(here::here("Data/trawl/clf_1907RL.rds")) %>% 
  mutate(survey = "1907RL")
clf.2021 <- readRDS(here::here("Data/trawl/clf_2107RL.rds")) %>% 
  mutate(survey = "2107RL")
clf.2022 <- readRDS(here::here("Data/trawl/clf_2207RL.rds")) %>% 
  mutate(survey = "2207RL")

# Combine specimen data
clf.all <- clf.2019 %>% bind_rows(clf.2021) %>% bind_rows(clf.2022)

# Replace missing values
clf.all <- clf.all %>% 
  # sample type with "Trawl"
  replace_na(list(sample.type = "Trawl")) %>%
  # round herring proportions with zeros (only present in 2021)
  replace_na(list(prop.rher = 0, prop.rher.wg = 0)) %>%
  # Remove purse seine sets
  filter(sample.type == "Trawl")

saveRDS(clf.all, here::here("Data/trawl/clf_all.rds"))

clf.zero <- filter(clf.all, CPS.num == 0)
clf.pos <- filter(clf.all, CPS.num > 0) %>% 
  arrange(X)


