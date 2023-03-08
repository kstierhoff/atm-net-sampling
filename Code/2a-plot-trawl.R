# pacman::p_load(tidyverse, lubridate, here, scatterpie, patchwork, sf)
# source(here::here("Doc/settings.R"))

# Load map data
load(here("Data/map/basemap.Rdata"))

# Calculate pie radius based on latitude range
# Use clf data to resize map, if needed
map.bounds <- clf.all %>% 
  st_as_sf(coords = c("long", "lat"), crs = crs.geog) %>% 
  st_transform(crs = crs.proj) %>%
  st_bbox()

# Set pie scale
pie.radius <- as.numeric(abs(map.bounds$ymin - map.bounds$ymax)*pie.scale)

# Plot trawl data
## Pie charts-catch weight
# Summarise catch by cluster
cluster.pie <- haul.summ.wt %>% 
  select(-haul, -AllCPS) %>% 
  group_by(survey, cluster) %>% 
  summarise_all(list(sum)) %>% 
  ungroup() %>% 
  mutate(AllCPS = rowSums(select(., -survey, -cluster)),
         r = pie.radius) %>% 
  replace(is.na(.), 0) %>% 
  # Add lat/long
  left_join(select(clf.all, survey, cluster, lat, long, X, Y))

cluster.zero <- filter(cluster.pie, AllCPS == 0) 
cluster.pos <- filter(cluster.pie, AllCPS > 0) %>% 
  arrange(X)

## Pie charts-catch weight proportion
cluster.wt.grid <- base.map +
  # Plot trawl pies
  geom_scatterpie(data = cluster.pos, 
                  aes(X, Y, group = cluster, r = pie.radius),
                  cols = c("Anchovy", "JackMack", "PacHerring", 
                           "PacMack", "RndHerring", "Sardine"),
                  alpha = 0.8, sorted_by_radius = TRUE) +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "J. Mackerel", "P. herring", 
                               "P. mackerel", "R. herring", "Sardine"),
                    values = c(anchovy.color, jack.mack.color, pac.herring.color, 
                               pac.mack.color, rnd.herring.color, sardine.color)) +
  # Plot empty cluster locations
  geom_point(data = cluster.zero, aes(X, Y),
             size = 3, shape = 21, fill = 'black', colour = 'white') +
  ggtitle("Catch Proportions (Weight) by Trawl Cluster") +
  # Facet by survey
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"], map.bounds["xmax"])), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

ggsave(cluster.wt.grid, filename = here("Figs/cluster_grid.png"),
       height = 10, width = 15)

## Pie charts-acoustic proportion
clf.grid <- base.map +
  # Plot trawl pies
  geom_scatterpie(data = clf.pos, 
                  aes(X, Y, group = cluster, r = pie.radius),
                  cols = c("prop.anch.wg","prop.her.wg","prop.jack.wg",
                           "prop.mack.wg", "prop.rher.wg", "prop.sar.wg"),
                  alpha = 0.8) +
  # Plot empty trawl locations
  geom_point(data = clf.zero, aes(X, Y),
             size = 3, shape = 21, fill = 'black', colour = 'white') +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "P. herring", "J. mackerel",
                               "P. mackerel", "R. herring", "Sardine"),
                    values = c(anchovy.color, pac.herring.color, jack.mack.color,  
                               pac.mack.color, rnd.herring.color, sardine.color)) +
  ggtitle("Acoustic Proportions (Weight) by Trawl Cluster") +
  # Facet by survey
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, 
           xlim = unname(c(map.bounds["xmin"], map.bounds["xmax"])), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

ggsave(clf.grid, filename = here("Figs/clf_grid.png"),
       height = 10, width = 15)

# Format catch data for plotting
haul.summ.wt <- catch.trawl %>%
  select(survey, haul, cluster, scientificName, totalWeight) %>%
  tidyr::spread(scientificName, totalWeight) %>% 
  replace(is.na(.), 0) %>% 
  mutate(AllCPS = rowSums(select(., -survey, -haul, -cluster))) %>%
  rename("PacHerring" = "Clupea pallasii",
         "Anchovy"    = "Engraulis mordax",
         "Sardine"    = "Sardinops sagax",
         "PacMack"    = "Scomber japonicus",
         "JackMack"   = "Trachurus symmetricus",
         "RndHerring" = "Etrumeus acuminatus")
