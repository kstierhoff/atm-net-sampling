# pacman::p_load(tidyverse, lubridate, here)
# pacman::p_load_github("kstierhoff/atmData")

# import.data <- TRUE

# Backscatter data
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
  saveRDS(nasc, here::here("Data/backscatter/nasc_core.rds"))
} else {
  # Load processed data
  nasc_core <- readRDS(here::here("Data/backscatter/nasc_core.rds"))
}
