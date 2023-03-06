# pacman::p_load(tidyverse, lubridate, here)
# pacman::p_load_github("kstierhoff/atmData")

# import.data <- TRUE

# Backscatter data
## Import backscatter data from {atmData}
nasc.2019 <- nasc_density_1907RL %>% 
  select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
         lat, long, cluster, cluster.distance, Region, 
         NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
  mutate(survey = "1907RL")
nasc.2021 <- nasc_density_2107RL %>% 
  select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
         lat, long, cluster, cluster.distance, Region, 
         NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
  mutate(survey = "2107RL")
nasc.2022 <- nasc_density_2207RL %>% 
  select(transect, transect.name, Interval, datetime, vessel.name, vessel.orig, 
         lat, long, cluster, cluster.distance, Region, 
         NASC, cps.nasc, starts_with("prop."), ends_with(".dens")) %>% 
  mutate(survey = "2207RL")

nasc <- nasc.2019 %>% bind_rows(nasc.2021) %>% bind_rows(nasc.2022) %>% 
  select(survey, transect:cps.nasc, starts_with("prop."), ends_with(".dens"))
