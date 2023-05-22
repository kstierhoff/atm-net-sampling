# Load trawl data
catch.trawl   <- readRDS(here::here("Data/trawl/catch_1907RL.rds"))
lengths.trawl <- readRDS(here::here("Data/trawl/lengths_1907RL.rds"))
clf.trawl     <- readRDS(here::here("Data/trawl/clf_1907RL.rds"))

# Load seine data
sets.seine    <- readRDS(here::here("Data/seine/sets_1907RL.rds"))
lengths.seine <- readRDS(here::here("Data/seine/specimens_1907RL.rds"))
clf.seine     <- readRDS(here::here("Data/seine/clf_seine_1907RL.rds"))

# Load backscatter
nasc.core <- readRDS(here::here("Data/backscatter/nasc_core_1907RL.rds"))
nasc.ns   <- readRDS(here::here("Data/backscatter/nasc_nearshore_1907RL.rds"))