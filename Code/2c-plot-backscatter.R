# pacman::p_load(tidyverse, lubridate, here, scatterpie, patchwork, sf)
# pacman::p_load_gh("kstierhoff/atm")
# source(here::here("Doc/settings.R"))

# Load map data
if (!exists("base.map")) load(here("Data/map/basemap.Rdata"))

# Calculate pie radius based on latitude range
# Use clf data to resize map, if needed
map.bounds <- nasc_core %>% 
  st_as_sf(coords = c("long", "lat"), crs = crs.geog) %>% 
  st_transform(crs = crs.proj) %>%
  st_bbox()

# Summarize backscatter for plotting
# Average cps.nasc over defined interval
# Summarize by filename, not transect, so that renamed (i.e., strip.tx.chars == TRUE) transects get included.
nasc.plot <- nasc.core %>%
  select(survey, filename, transect.name, transect, int, datetime, lat, long, cps.nasc) %>% 
  group_by(survey, filename, transect.name, transect, int) %>% 
  summarise(
    lat  = lat[1],
    long = long[1],
    NASC = mean(cps.nasc)) %>% 
  # Create bins for defining point size in NASC plots
  mutate(bin       = cut(NASC, nasc.breaks, include.lowest = TRUE),
         bin.level =  as.numeric(bin)) %>% 
  # st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
  ungroup() %>% 
  project_df(to = crs.proj)

nasc.ns.plot <- nasc.ns %>%
  select(survey, filename, transect.name, transect, int, datetime, lat, long, cps.nasc) %>% 
  group_by(survey, filename, transect.name, transect, int) %>% 
  summarise(
    lat  = lat[1],
    long = long[1],
    NASC = mean(cps.nasc)) %>% 
  # Create bins for defining point size in NASC plots
  mutate(bin       = cut(NASC, nasc.breaks, include.lowest = TRUE),
         bin.level =  as.numeric(bin)) %>% 
  # st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
  ungroup() %>% 
  project_df(to = crs.proj)

# Select plot levels for backscatter data
nasc.levels.all <- unique(nasc.plot$bin.level)
nasc.labels.all <- nasc.labels[sort(nasc.levels.all)]
nasc.sizes.all  <- nasc.sizes[sort(nasc.levels.all)]
nasc.colors.all <- nasc.colors[sort(nasc.levels.all)]

# Map backscatter
## Core
nasc.map <- base.map +
  # Plot NASC data
  geom_point(data = nasc.plot, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75) +
  # Configure size and color scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all,labels = nasc.labels.all) +
  # Configure legend guides
  guides(fill = guide_legend(), size = guide_legend()) +
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"], map.bounds["xmax"])), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

ggsave(nasc.map, filename = here("Figs/nasc_grid.png"),
       height = 10, width = 15)

## Core
nasc.ns.map <- base.map +
  # Plot NASC data
  geom_point(data = nasc.ns.plot, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75) +
  # Configure size and color scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all,labels = nasc.labels.all) +
  # Configure legend guides
  guides(fill = guide_legend(), size = guide_legend()) +
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"], map.bounds["xmax"])), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

ggsave(nasc.ns.map, filename = here("Figs/nasc_ns_grid.png"),
       height = 10, width = 15)
