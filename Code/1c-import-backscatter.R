# pacman::p_load(tidyverse, lubridate, here)
# pacman::p_load_gh("kstierhoff/atmData")

# import.data <- TRUE

# Backscatter data
## Core area
if (import.data) {
  ## Import backscatter data from {atmData}
  nasc.2019 <- nasc_density_1907RL %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance, Region, 
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "1907RL")
  saveRDS(nasc.2019, here::here("Data/backscatter/nasc_core_1907RL.rds"))
  
  nasc.2021 <- nasc_density_2107RL %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance, Region, 
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "2107RL")
  saveRDS(nasc.2021, here::here("Data/backscatter/nasc_core_2107RL.rds"))
  
  nasc.2022 <- nasc_density_2207RL %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance, Region, 
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "2207RL")
  saveRDS(nasc.2022, here::here("Data/backscatter/nasc_core_2207RL.rds"))
  
  # Save data
  nasc_core <- nasc.2019 %>% bind_rows(nasc.2021) %>% bind_rows(nasc.2022) %>% 
    select(survey, transect:cps.nasc, starts_with("prop."), ends_with(".dens"))
  saveRDS(nasc_core, here::here("Data/backscatter/nasc_core.rds"))
} else {
  # Load processed data
  nasc_core <- readRDS(here::here("Data/backscatter/nasc_core.rds"))
}

## Nearshore
if (import.data) {
  ## Import backscatter data from {atmData}
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/nasc_nearshore_final.Rdata")
  nasc.ns.2019 <- nasc.nearshore %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance,  
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "1907RL")
  
  saveRDS(nasc.ns.2019, here::here("Data/backscatter/nasc_nearshore_1907RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/nasc_nearshore_final.Rdata")
  nasc.ns.2021 <- nasc.nearshore %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance,  
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "2107RL")
  
  saveRDS(nasc.ns.2021, here::here("Data/backscatter/nasc_nearshore_2107RL.rds"))
  
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/nasc_nearshore_final.Rdata")
  nasc.ns.2022 <- nasc.nearshore %>% 
    select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
           lat, long, cluster, cluster.distance,  
           NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
    mutate(survey = "2207RL")
  
  saveRDS(nasc.ns.2022, here::here("Data/backscatter/nasc_nearshore_2207RL.rds"))
  
  # Save data
  nasc_ns <- nasc.ns.2019 %>% bind_rows(nasc.ns.2021) %>% bind_rows(nasc.ns.2022) %>% 
    select(survey, transect:cps.nasc, starts_with("prop."), ends_with(".dens"))
  saveRDS(nasc_ns, here::here("Data/backscatter/nasc_nearshore.rds"))
} else {
  # Load processed data
  nasc_ns <- readRDS(here::here("Data/backscatter/nasc_nearshore.rds"))
}
