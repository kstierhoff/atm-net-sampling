# Summarise catch by weight -----------------------------------------------
set.summ.wt <- catch.all %>% 
  # ungroup() %>% 
  # filter(vessel.name %in% seine.vessels) %>% 
  select(survey.name, vessel.name, key.set, scientificName, totalWtKg) %>% 
  tidyr::spread(scientificName, totalWtKg) 

# Add species with zero total weight
if (!has_name(set.summ.wt, "Engraulis mordax"))      {set.summ.wt$`Engraulis mordax`      <- 0}
if (!has_name(set.summ.wt, "Sardinops sagax"))       {set.summ.wt$`Sardinops sagax`       <- 0}
if (!has_name(set.summ.wt, "Scomber japonicus"))     {set.summ.wt$`Scomber japonicus`     <- 0}
if (!has_name(set.summ.wt, "Trachurus symmetricus")) {set.summ.wt$`Trachurus symmetricus` <- 0}
if (!has_name(set.summ.wt, "Clupea pallasii"))       {set.summ.wt$`Clupea pallasii`       <- 0}
# if (!has_name(set.summ.wt, "Atherinopsis californiensis")) {set.summ.wt$`Atherinopsis californiensis` <- 0}

# Calculate total weight of all CPS species
set.summ.wt <- set.summ.wt %>%
  ungroup() %>% 
  replace(is.na(.), 0) %>% 
  mutate(AllCPS = rowSums(select(., -(survey.name:key.set)))) %>%
  rename("PacHerring" = "Clupea pallasii",
         "Anchovy"    = "Engraulis mordax",
         "Sardine"    = "Sardinops sagax",
         "PacMack"    = "Scomber japonicus",
         "JackMack"   = "Trachurus symmetricus") %>% 
  left_join(select(sets.all, key.set, lat, long))

set.pie <- set.summ.wt %>% 
  project_df(to = 3310) 

set.pos <- filter(set.pie, AllCPS > 0) %>% 
  arrange(desc(X))

set.zero <- filter(set.pie, AllCPS == 0)

map.bounds.ns <- set.pie %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = 3310) %>%
  st_bbox()

pie.radius.ns <- as.numeric(abs(map.bounds.ns$ymin - map.bounds.ns$ymax)*pie.scale)

# Plot seine catch weight
ggplot() + 
  # Plot trawl pies
  geom_scatterpie(data = set.pos, 
                  aes(X, Y, group = key.set, r = pie.radius.ns),
                  cols = c("Anchovy","JackMack",
                           "PacHerring","PacMack","Sardine"),
                  color = 'black', alpha = 0.8) +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "J. Mackerel", 
                               "P. herring", "P. mackerel", "Sardine"),
                    values = c(anchovy.color, jack.mack.color, 
                               pac.herring.color, pac.mack.color, sardine.color)) +
  # Plot empty trawl locations
  geom_point(data = set.zero, aes(X, Y), 
             size = 2, shape = 21, fill = 'black', colour = 'white') +
  facet_wrap(~survey.name) + 
  coord_sf(crs = 3310)

# Plot seine acoustic proportions
ggplot() + 
  # Plot set pies
  geom_scatterpie(data = clf.seine.all, 
                  aes(X, Y, group = cluster, r = pie.radius.ns),
                  cols = c("prop.anch","prop.her","prop.jack",
                           "prop.mack", "prop.sar"),
                  color = 'black', alpha = 0.8) +
  # Configure scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "P. herring", "J. mackerel",
                               "P. mackerel", "Sardine"),
                    values = c(anchovy.color, pac.herring.color, jack.mack.color,  
                               pac.mack.color,  sardine.color)) +
  # # Plot empty trawl locations
  # geom_point(data = set.zero, aes(X, Y), 
  #            size = 2, shape = 21, fill = 'black', colour = 'white') +
  facet_wrap(~survey.name) + 
  coord_sf(crs = 3310)
