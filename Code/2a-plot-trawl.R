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


## Pie charts-acoustic proportion
clf.grid <- base.map +
  # Plot trawl pies
  geom_scatterpie(data = clf.pos, 
                  aes(X, Y, group = cluster, r = pie.radius, colour = sample.type),
                  cols = c("prop.anch.wg","prop.her.wg","prop.jack.wg",
                           "prop.mack.wg", "prop.sar.wg"),
                  alpha = 0.8) +
  # Plot empty trawl locations
  geom_point(data = clf.zero, aes(X, Y),
             size = 3, shape = 21, fill = 'black', colour = 'white') +
  # Configure pie outline colors
  scale_colour_manual(name = "Sample type",
                      labels = c("Purse seine", "Trawl"),
                      values = c(seine.color, trawl.color),
                      guide = "none") +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "P. herring", "J. mackerel",
                               "P. mackerel", "Sardine"),
                    values = c(anchovy.color, pac.herring.color, jack.mack.color,  
                               pac.mack.color, sardine.color)) +
  # Facet by survey
  facet_grid(survey ~ sample.type) +
  coord_sf(crs = crs.proj, 
           xlim = unname(c(map.bounds["xmin"], map.bounds["xmax"])), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

ggsave(clf.grid, filename = here("Figs/clf_grid.png"),
       height = 15, width = 6)
