# Load processed data
## Core
nasc.core <- readRDS(here::here("Data/backscatter/nasc_core.rds"))
## Nearshore (not sure if needed)
nasc.ns <- readRDS(here::here("Data/backscatter/nasc_nearshore.rds"))

# nasc.lim <- filter(nasc.ns, survey == "1907RL", vessel.name == "LM") %>% 
#   select(lat) %>% 
#   range()

nasc.lim <- filter(set.pie, survey.name == "1907RL", vessel.name == "Lisa Marie") %>% 
  select(lat) %>% 
  range()

# ggplot() + 
#   geom_point(data = filter(nasc.core, survey == "1907RL"), aes(long, lat), colour = "red") + 
#   geom_point(data = filter(nasc.ns, survey == "1907RL", vessel.name == "LM"), 
#              aes(long, lat), fill = NA, colour = "blue", shape = 21) +
#   coord_map()

nasc.core <- nasc.core %>% 
  filter(survey == "1907RL", between(lat, min(nasc.lim) - 0.5, max(nasc.lim) + 1), transect < 107)

# Plot core and nearshore backscatter; purse seine set locations
ggplot() +
  geom_point(data = nasc.core.sub, aes(long, lat, colour = "red")) +
  geom_point(data = filter(nasc.ns, survey == "1907RL", vessel.name == "LM"),
             aes(long, lat), fill = NA, colour = "blue", shape = 21) +
  geom_point(data = filter(set.pie, survey.name == "1907RL", vessel.name == "Lisa Marie"),
             aes(long, lat), colour = "black", size = 4) +
  coord_map()

# Summarize backscatter for plotting
# Average cps.nasc over defined interval
# Summarize by filename, not transect, so that renamed (i.e., strip.tx.chars == TRUE) transects get included.
nasc.plot.sub <- nasc.core %>%
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

# Subset clf files for the ROI
clf.zero.sub <- filter(clf.all, survey == "1907RL", between(lat, min(nasc.lim) - 0.5, max(nasc.lim) + 1), CPS.num == 0)
clf.pos.sub  <- filter(clf.all, survey == "1907RL", between(lat, min(nasc.lim) - 0.5, max(nasc.lim) + 1), CPS.num > 0) %>% 
  arrange(X)

map.bounds.sub <- nasc.plot.sub %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = 3310) %>%
  st_bbox()

pie.radius.sub <- as.numeric(abs(map.bounds.sub$ymin - map.bounds.sub$ymax)*pie.scale)

# Load map data
if (!exists("base.map")) load(here("Data/map/basemap.Rdata"))

# Select plot levels for backscatter data
nasc.levels.all <- unique(nasc.plot.sub$bin.level)
nasc.labels.all <- nasc.labels[sort(nasc.levels.all)]
nasc.sizes.all  <- nasc.sizes[sort(nasc.levels.all)]
nasc.colors.all <- nasc.colors[sort(nasc.levels.all)]

# Map backscatter
## Core w/ trawls
nasc.core.trawl <- base.map +
  # Plot NASC data
  geom_point(data = nasc.plot.sub, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75) +
  # Configure size and color scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all,labels = nasc.labels.all) +
  # Configure legend guides
  guides(fill = guide_legend(), size = guide_legend()) +
  ggnewscale::new_scale_fill() +
  # Plot trawl pies
  geom_scatterpie(data = clf.pos.sub, 
                  aes(X, Y, group = cluster, r = pie.radius.sub),
                  cols = c("prop.anch","prop.jack",
                           "prop.her","prop.mack","prop.sar"),
                  color = 'black', alpha = 0.8) +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "J. Mackerel", 
                               "P. herring", "P. mackerel", "Sardine"),
                    values = c(anchovy.color, jack.mack.color, 
                               pac.herring.color, pac.mack.color, sardine.color)) +
  # Plot empty trawl locations
  geom_point(data = clf.zero.sub, aes(X, Y), 
             size = 2, shape = 21, fill = 'black', colour = 'white') +
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds.sub["xmin"]-150000, map.bounds.sub["xmax"]+100000)), 
           ylim = unname(c(map.bounds.sub["ymin"], map.bounds.sub["ymax"])))

## Core w/ purse seine sets
nasc.core.sets <- base.map +
  # Plot NASC data
  geom_point(data = nasc.plot.sub, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75) +
  # Configure size and color scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all,labels = nasc.labels.all) +
  # Configure legend guides
  guides(fill = guide_legend(), size = guide_legend()) +
  ggnewscale::new_scale_fill() +
  # Plot trawl pies
  geom_scatterpie(data = set.pos, 
                  aes(X, Y, group = key.set, r = pie.radius.sub),
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
  facet_wrap(~survey) +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "bold")) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds.sub["xmin"]-150000, map.bounds.sub["xmax"]+100000)), 
           ylim = unname(c(map.bounds.sub["ymin"], map.bounds.sub["ymax"])))

library(patchwork)
nasc.core.combo <- nasc.core.trawl + nasc.core.sets

nasc.core.combo

# nasc.core.sets

