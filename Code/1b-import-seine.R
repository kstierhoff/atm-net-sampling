# Set data-------------------------------------------------
## Import set data
if (import.data) {
  # 2019
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/purse_seine_sets.Rdata")
  sets.lm.2019 <- lm.sets %>% 
    mutate(vessel.name = "LM",
           datetime = ymd_hms(paste(date, time)),
           survey.name = "1907RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.lbc.2019 <- lbc.sets %>% 
    mutate(vessel.name = "LBC",
           datetime = mdy_hm(datetime),
           survey.name = "1907RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.2019 <- bind_rows(sets.lm.2019, sets.lbc.2019)
  
  saveRDS(sets.2019, here::here("Data/seine/sets_1907RL.rds"))
  
  # 2021
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/purse_seine_sets.Rdata")
  sets.lm.2021 <- lm.sets %>% 
    mutate(vessel.name = "LM",
           datetime = ymd_hms(paste(date, time)),
           survey.name = "2107RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.lbc.2021 <- lbc.sets %>% 
    mutate(vessel.name = "LBC",
           datetime = mdy_hm(datetime),
           survey.name = "2107RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.2021 <- bind_rows(sets.lm.2021, sets.lbc.2021)
  
  saveRDS(sets.2021, here::here("Data/seine/sets_2107RL.rds"))
  
  # 2022
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/purse_seine_sets.Rdata")
  sets.lm.2022 <- lm.sets %>% 
    mutate(vessel.name = "LM",
           datetime = ymd_hms(paste(date, time)),
           survey.name = "2207RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.lbc.2022 <- lbc.sets %>% 
    mutate(vessel.name = "LBC",
           datetime = ymd_hms(paste(date, time)),
           survey.name = "2207RL") %>% 
    select(survey.name, key.set, datetime, lat, long, vessel.name)
  
  sets.2022 <- bind_rows(sets.lm.2022, sets.lbc.2022)
  
  saveRDS(sets.2022, here::here("Data/seine/sets_2207RL.rds"))
  
  # ggplot(sets.2022, aes(long, lat, colour = vessel.name)) + geom_point()
} else {
  sets.2019 <- readRDS(here::here("Data/seine/sets_1907RL.rds"))
  sets.2021 <- readRDS(here::here("Data/seine/sets_2107RL.rds"))
  sets.2022 <- readRDS(here::here("Data/seine/sets_2207RL.rds"))
}

sets.all <- sets.2019 %>% bind_rows(sets.2021) %>% bind_rows(sets.2022) 
# ggplot(sets.all, aes(long, lat, colour = vessel.name)) + geom_point() + facet_wrap(~survey.name)

# Catch data------------------------------------------------
## Import catch data
if (import.data) {
  # 2019
  load("C:/KLS/CODE/Github/estimATM/1907RL/Output/purse_seine_specimens.Rdata")
  specimens.lm.2019 <- lm.specimens %>% 
    mutate(survey.name = "1907RL") %>% 
    select(vessel.name = vessel_name, survey.name, key.set, scientificName, weightg, 
           standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.lbc.2019 <- lbc.specimens %>% 
    mutate(survey.name = "1907RL") %>% 
    select(vessel.name, survey.name, key.set, scientificName, weightg, standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.2019 <- bind_rows(specimens.lm.2019, specimens.lbc.2019)
  
  saveRDS(specimens.2019, here::here("Data/seine/specimens_1907RL.rds"))
  
  # 2021
  load("C:/KLS/CODE/Github/estimATM/2107RL/Output/purse_seine_specimens.Rdata")
  specimens.lm.2021 <- lm.specimens %>% 
    mutate(survey.name = "2107RL") %>% 
    select(vessel.name = vessel_name, survey.name, key.set, scientificName, weightg, 
           standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.lbc.2021 <- lbc.specimens %>% 
    mutate(survey.name = "2107RL") %>% 
    select(vessel.name, survey.name, key.set, scientificName, weightg, standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.2021 <- bind_rows(specimens.lm.2021, specimens.lbc.2021)
  
  saveRDS(specimens.2021, here::here("Data/seine/specimens_2107RL.rds"))
  
  # 2022
  load("C:/KLS/CODE/Github/estimATM/2207RL/Output/purse_seine_specimens.Rdata")
  specimens.lm.2022 <- lm.lengths %>% 
    mutate(survey.name = "2207RL") %>% 
    select(vessel.name = vessel_name, survey.name, key.set, scientificName, weightg, 
           standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.lbc.2022 <- lbc.lengths %>% 
    mutate(survey.name = "2207RL") %>% 
    select(vessel.name, survey.name, key.set, scientificName, weightg, standardLength_mm, forkLength_mm, totalLength_mm)
  
  specimens.2022 <- bind_rows(specimens.lm.2022, specimens.lbc.2022)
  
  saveRDS(specimens.2022, here::here("Data/seine/specimens_2207RL.rds"))
  
  # ggplot(specimens.2022, aes(totalLength_mm, weightg, colour = vessel.name)) +
  #   geom_point(alpha = 0.5) +
  #   facet_wrap(~scientificName, scales = "free")
  
} else {
  specimens.2019 <- readRDS(here::here("Data/seine/specimens_1907RL.rds"))
  specimens.2021 <- readRDS(here::here("Data/seine/specimens_2107RL.rds"))
  specimens.2022 <- readRDS(here::here("Data/seine/specimens_2207RL.rds"))
}

specimens.all <- specimens.2019 %>% bind_rows(specimens.2021) %>% bind_rows(specimens.2022) 

# ggplot(specimens.all, aes(totalLength_mm, weightg, colour = vessel.name)) +
#   geom_point(alpha = 0.5) +
#   facet_grid(survey.name~scientificName, scales = "free")

catch.all <- specimens.all %>% 
  group_by(survey.name, vessel.name, key.set, scientificName) %>% 
  summarise(totalWtKg = sum(weightg)/1000,
            totalCount = n()) %>% 
  left_join(select(sets.all, key.set, lat, long))

# ggplot(catch.all, aes(long, lat, size = totalWtKg, colour = vessel.name)) + 
#   geom_point(shape = 21) + 
#   coord_map() + 
#   facet_grid(survey.name~scientificName)


# Specimen data----------------------------------------------
## Import specimen data

  

# Import clf data -------------------------------------------
# "Output/clf_ts_proportions_seine.csv"
## Import set data
if (import.data) {
  # 2019
  clf.seine.2019 <- read_csv("C:/KLS/CODE/Github/estimATM/1907RL/Output/clf_ts_proportions_seine.csv") %>% 
    mutate(survey.name = "1907RL")
  # 2021
  clf.seine.2021 <- read_csv("C:/KLS/CODE/Github/estimATM/2107RL/Output/clf_ts_proportions_seine.csv") %>% 
    mutate(survey.name = "2107RL")
  # 2022
  clf.seine.2022 <- read_csv("C:/KLS/CODE/Github/estimATM/2207RL/Output/clf_ts_proportions_seine.csv") %>% 
    mutate(survey.name = "2207RL")
  
  saveRDS(clf.seine.2019, here::here("Data/seine/clf_seine_1907RL.rds"))
  saveRDS(clf.seine.2021, here::here("Data/seine/clf_seine_2107RL.rds"))
  saveRDS(clf.seine.2022, here::here("Data/seine/clf_seine_2207RL.rds"))
  
  # ggplot(sets.2022, aes(long, lat, colour = vessel.name)) + geom_point()
} else {
  clf.seine.2019 <- readRDS(here::here("Data/seine/clf_seine_1907RL.rds"))
  clf.seine.2021 <- readRDS(here::here("Data/seine/clf_seine_2107RL.rds"))
  clf.seine.2022 <- readRDS(here::here("Data/seine/clf_seine_2207RL.rds"))
}

clf.seine.all <- clf.seine.2019 %>% bind_rows(clf.seine.2021) %>% bind_rows(clf.seine.2022) 
