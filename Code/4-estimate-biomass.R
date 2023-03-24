pacman::p_load(tidyverse, lubridate, here, scatterpie, patchwork, 
               sf, geosphere, mapview, shadowtext, mapview, swfscMisc)
pacman::p_load_gh("kstierhoff/atmData")
pacman::p_load_gh("kstierhoff/atm")

theme_set(theme_bw())

# Load settings
source(here::here("Doc/settings.R"))

# Load 2019 backscatter
nasc <- readRDS(here("Data/backscatter/nasc_core_1907RL.rds")) %>% 
  filter(vessel.name %in% c("RL"))

# ggplot(nasc, aes(long, lat, colour = factor(cluster))) +
#   geom_point(show.legend = FALSE) +
#   facet_wrap(~vessel.name) +
#   coord_sf()

# Load 2019 clf
clf <- readRDS(here("Data/trawl/clf_1907RL.rds"))

# Load 2019 clf.seine
clf.ns <- readRDS(here("Data/seine/clf_seine_1907RL.rds"))

# Quick plot
ggplot() +
  geom_point(data = nasc, aes(long, lat, size = log(cps.nasc+1)), 
             colour = "red", alpha = 0.2) +
  geom_point(data = clf, aes(long, lat), 
             colour = "blue", alpha = 0.7) +
  geom_point(data = clf.ns, aes(long, lat), 
             colour = "green", alpha = 0.7) +
  # geom_sf(data = clf.match, shape = 21, fill = NA, colour = "white") +
  # geom_sf(data = clf.ns.match, shape = 21, fill = NA, colour = "white") +
  coord_sf()

# Match nasc to clf
nasc.match <- nasc %>% 
  st_as_sf(coords = c("long","lat"), crs = crs.geog)

clf.match <- clf %>%
  filter(CPS.num > 0) %>% 
  st_as_sf(coords = c("long","lat"), crs = crs.geog)

clf.ns.match <- clf.ns %>%
  filter(CPS.num > 0) %>% 
  st_as_sf(coords = c("long","lat"), crs = crs.geog)

# Find nearest cluster ----------------------------
# Returns a vector of nearest clusters
nearest.cluster <- st_nearest_feature(nasc.match, clf.match)
nearest.cluster.ns <- st_nearest_feature(nasc.match, clf.ns.match)

# Expand clf to match nasc ------------------------
cluster.sp <- clf.match[nearest.cluster, ] %>% 
  select(geometry) %>% 
  as_Spatial()

cluster.ns.sp <- clf.ns.match[nearest.cluster.ns, ] %>% 
  select(geometry) %>% 
  as_Spatial()

# Make nasc sp
# Removing the other data (i.e., only retaining the geometry) decreases the size of nasc from 50 to 1 MB
nasc.sp <- nasc.match %>% 
  select(geometry) %>% 
  as_Spatial()

# Compute distances with {geosphere}
## Must be done in WGS84 projection
nasc.match <- nasc.match %>% 
  mutate(cluster = clf.match$cluster[nearest.cluster],
         cluster.distance = geosphere::distGeo(nasc.sp, cluster.sp)*0.000539957) 

nasc.match.ns <- nasc.match %>% 
  mutate(cluster = clf.ns.match$cluster[nearest.cluster.ns],
         cluster.distance = geosphere::distGeo(nasc.sp, cluster.ns.sp)*0.000539957) 

# Create cluster polygons
nasc.super.clusters <- nasc.match %>% 
  group_by(cluster) %>% 
  summarise() %>% 
  st_convex_hull()

nasc.super.clusters.ns <- nasc.match.ns %>% 
  group_by(cluster) %>% 
  summarise() %>% 
  st_convex_hull()

# Remove geometry
nasc.match    <- st_set_geometry(nasc.match, NULL) 
nasc.match.ns <- st_set_geometry(nasc.match.ns, NULL) 

# Add clusters and cluster distances to nasc from core and nearshore area
nasc <- nasc %>% 
  # select(-cluster, -cluster.distance) %>% 
  bind_cols(select(nasc.match, cluster.core = cluster, cluster.distance.core = cluster.distance)) %>% 
  bind_cols(select(nasc.match.ns, cluster.ns = cluster, cluster.distance.ns = cluster.distance)) 

nasc.lim <- filter(clf.ns) %>% 
  select(lat) %>% 
  range()

nasc <- nasc %>% 
  filter(between(lat, min(nasc.lim) - 0.5, max(nasc.lim) + 1), transect < 107)

# Plot nasc with corresponding core and nearshore clusters
# ggplot() +
#   geom_point(data = nasc, aes(long, lat, colour = factor(cluster.core)),
#              show.legend = FALSE) +
#   geom_sf(data = nasc.super.clusters, fill = NA) +
# #   geom_point(data = filter(nasc, cluster != cluster.core), aes(long, lat)) +
#   coord_sf()
# 
# ggplot() +
#   geom_point(data = nasc, aes(long, lat, colour = factor(cluster.ns)),
#              show.legend = FALSE) +
#   geom_sf(data = nasc.super.clusters.ns, fill = NA) +
#   coord_sf()

# Use retained backscatter data to resize map, if needed
map.bounds <- nasc %>% 
  st_as_sf(coords = c("long", "lat"), crs = crs.geog) %>% 
  st_transform(crs = crs.proj) %>%
  st_bbox()

nasc <- project_df(nasc, to = crs.proj)

if (!exists("base.map")) load(here("Data/map/basemap.Rdata"))

plot.core.clusters <- base.map +
  geom_point(data = nasc, aes(X, Y, colour = factor(cluster.core)),
             show.legend = FALSE) +
  geom_sf(data = nasc.super.clusters, fill = NA) +
  geom_shadowtext(data = filter(clf, CPS.num > 0), aes(X, Y, label = cluster, colour = factor(cluster)),
             bg.color = "white", fontface = "bold", show.legend = FALSE) +
  geom_shadowtext(data = filter(clf, CPS.num == 0), aes(X, Y, label = cluster),
                  colour = "gray50", bg.color = "white", show.legend = FALSE) +
  # geom_sf(data = filter(nasc.super.clusters, cluster %in% unique(nasc$cluster.core)), fill = NA) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

plot.ns.clusters <- base.map +
  geom_point(data = nasc, aes(X, Y, colour = factor(cluster.ns)),
             show.legend = FALSE) +
  geom_sf(data = nasc.super.clusters.ns, fill = NA) +
  geom_shadowtext(data = filter(clf.ns, CPS.num > 0), aes(X, Y, label = cluster, colour = factor(cluster)),
                  bg.color = "white", fontface = "bold", show.legend = FALSE) +
  geom_shadowtext(data = filter(clf.ns, CPS.num == 0), aes(X, Y, label = cluster),
                  colour = "gray50", bg.color = "white", show.legend = FALSE) +
  # geom_sf(data = filter(nasc.super.clusters, cluster %in% unique(nasc$cluster.core)), fill = NA) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)), 
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

cluster.plot.combo <- plot.core.clusters + plot.ns.clusters +
  plot_annotation(tag_levels = c('a'), tag_suffix = ')')

ggsave(cluster.plot.combo, filename = here("figs/fig_cluster_assignments_combo.png"),
       height = 7, width = 8)
# 
# hist.dist.core <- ggplot(nasc, aes(cluster.distance)) + geom_histogram()
# hist.dist.ns   <- ggplot(nasc, aes(cluster.distance.ns)) + geom_histogram()

# TODO-------------------------------------------------------------------------
# Create separate nasc dfs with core and nearshore proportions
nasc.core.final <- nasc %>% 
  select(transect, Interval, int, cluster.core, lat, long, X, Y, NASC, cps.nasc) %>% 
  left_join(select(clf, cluster.core = cluster, starts_with("sigmawg."), starts_with("prop."))) %>% 
  mutate( 
    anch.dens = cps.nasc*prop.anch / (4*pi*sigmawg.anch) / 1000,
    her.dens  = cps.nasc*prop.her  / (4*pi*sigmawg.her)  / 1000,
    jack.dens = cps.nasc*prop.jack / (4*pi*sigmawg.jack) / 1000,
    mack.dens = cps.nasc*prop.mack / (4*pi*sigmawg.mack) / 1000,
    sar.dens  = cps.nasc*prop.sar  / (4*pi*sigmawg.sar)  / 1000)

nasc.ns.final <- nasc %>% 
  select(transect, Interval, int, cluster.ns, lat, long, X, Y, NASC, cps.nasc) %>% 
  left_join(select(clf.ns, cluster.ns = cluster, starts_with("sigmawg."), starts_with("prop."))) %>% 
  mutate( 
    anch.dens = cps.nasc*prop.anch / (4*pi*sigmawg.anch) / 1000,
    her.dens  = cps.nasc*prop.her  / (4*pi*sigmawg.her)  / 1000,
    jack.dens = cps.nasc*prop.jack / (4*pi*sigmawg.jack) / 1000,
    mack.dens = cps.nasc*prop.mack / (4*pi*sigmawg.mack) / 1000,
    sar.dens  = cps.nasc*prop.sar  / (4*pi*sigmawg.sar)  / 1000)

# Recompute density from prop and sigma

dens.core.final <- ggplot(nasc.core.final, aes(X, Y, colour = factor(cluster.core), size = cps.nasc)) +
  geom_point(show.legend = FALSE) + 
  geom_sf(data = nasc.super.clusters, 
          fill = NA, inherit.aes = FALSE) +
  # geom_sf(data = filter(nasc.super.clusters, cluster %in% unique(nasc.core.final$cluster)), 
  #         fill = NA, inherit.aes = FALSE) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))
dens.ns.final <- ggplot(nasc.ns.final, aes(X, Y, colour = factor(cluster.ns), size = cps.nasc)) +
  geom_point(show.legend = FALSE) + 
  geom_sf(data = nasc.super.clusters.ns, fill = NA, inherit.aes = FALSE) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

dens.final.combo <- dens.core.final + dens.ns.final +
  plot_annotation(tag_levels = c('a'), tag_suffix = ')')

ggsave(dens.final.combo, filename = here("figs/fig_density_final_combo.png"),
       height = 7, width = 8)

## i.e., re-join with clf and clf.ns
# Re-summarise nasc.density to define strata
## Core
nasc.density.core <- nasc.core.final %>%
  select(transect, int, lat, long, 
         anch.dens, her.dens, jack.dens, mack.dens, sar.dens) %>% 
  # When transects have multiple sections (e.g., 041-1, 041-2), interval is not unique
  # So create a int.key for computing density along unique intervals for plotting
  mutate(int.key = paste(transect, int)) %>% 
  # Then summarize by transect and unique interval
  group_by(transect, int.key) %>%
  summarise(
    lat = lat[1],
    long = long[1],
    `Engraulis mordax`      = mean(anch.dens),
    `Clupea pallasii`       = mean(her.dens),
    `Trachurus symmetricus` = mean(jack.dens),
    `Scomber japonicus`     = mean(mack.dens),
    `Sardinops sagax`       = mean(sar.dens)) %>% 
  gather(scientificName, density, -transect, -int.key, -lat, -long) %>% 
  mutate(bin       = cut(density, dens.breaks, include.lowest = TRUE),
         bin.level = as.numeric(bin))

## Nearshore
nasc.density.ns <- nasc.ns.final %>%
  select(transect, int, lat, long, 
         anch.dens, her.dens, jack.dens, mack.dens, sar.dens) %>% 
  # When transects have multiple sections (e.g., 041-1, 041-2), interval is not unique
  # So create a int.key for computing density along unique intervals for plotting
  mutate(int.key = paste(transect, int)) %>% 
  # Then summarize by transect and unique interval
  group_by(transect, int.key) %>%
  summarise(
    lat = lat[1],
    long = long[1],
    `Engraulis mordax`      = mean(anch.dens),
    `Clupea pallasii`       = mean(her.dens),
    `Trachurus symmetricus` = mean(jack.dens),
    `Scomber japonicus`     = mean(mack.dens),
    `Sardinops sagax`       = mean(sar.dens)) %>% 
  gather(scientificName, density, -transect, -int.key, -lat, -long) %>% 
  mutate(bin       = cut(density, dens.breaks, include.lowest = TRUE),
         bin.level = as.numeric(bin))

# Summarise biomass density by transect and species
nasc.density.summ <- nasc.density.core %>% 
  group_by(scientificName, transect) %>% 
  summarise(lat = lat[1], density = mean(density)) %>% 
  mutate(pos = case_when(
    density > 0 ~ TRUE,
    TRUE ~ FALSE))

nasc.density.ns.summ <- nasc.density.ns %>% 
  group_by(scientificName, transect) %>% 
  summarise(lat = lat[1], density = mean(density)) %>% 
  mutate(pos = case_when(
    density > 0 ~ TRUE,
    TRUE ~ FALSE))

dens.plot <- ggplot(nasc.density.summ, aes(lat, log(density+1))) + 
  geom_line() + 
  geom_text(aes(label = transect, colour = pos), show.legend=FALSE) +
  facet_wrap(~scientificName, nrow = 1) + coord_flip() +
  ggtitle("Core")

dens.plot.ns <- ggplot(nasc.density.ns.summ, aes(lat, log(density+1))) + 
  geom_line() + 
  geom_text(aes(label = transect, colour = pos),show.legend = FALSE) +
  facet_wrap(~scientificName, nrow = 1) + coord_flip()+
  ggtitle("Nearshore")

# dens.plot + dens.plot.ns

# Create new strata with area estimates
strata.final <- bind_rows(
  data.frame(
    scientificName = "Clupea pallasii",
    stratum = 1,
    transect = 65:69),
  data.frame(
    scientificName = "Clupea pallasii",
    stratum = 2,
    transect = 78:106),
  data.frame(
    scientificName = "Engraulis mordax",
    stratum = 1,
    transect = 79:100),
  data.frame(
    scientificName = "Sardinops sagax",
    stratum = 1,
    transect = 65:95),
  data.frame(
    scientificName = "Scomber japonicus",
    stratum = 1,
    transect = 65:95),
  data.frame(
    scientificName = "Trachurus symmetricus",
    stratum = 1,
    transect = 65:98)
)

strata.final.ns <- bind_rows(
  data.frame(
    scientificName = "Clupea pallasii",
    stratum = 1,
    transect = 76:106),
  data.frame(
    scientificName = "Engraulis mordax",
    stratum = 1,
    transect = 91:101),
  data.frame(
    scientificName = "Sardinops sagax",
    stratum = 1,
    transect = 81:93),
  data.frame(
    scientificName = "Trachurus symmetricus",
    stratum = 1,
    transect = 65:82)
)

# Create stratum polygons
# Draw pseudo-transects --------------------------------------------------------
# Get transect ends, calculate bearing, and add transect spacing
tx.ends <- nasc %>% 
  group_by(transect.name, transect, vessel.name) %>% 
  summarise(
    lat.i  = lat[which.max(long)],
    long.i = max(long),
    lat.o  = lat[which.min(long)],
    long.o = min(long)) %>% 
  # left_join(select(tx.nn, transect.name, spacing)) %>% 
  mutate(
    brg = swfscMisc::bearing(lat.i, long.i,
                             lat.o, long.o)[1])

# Get original inshore transect ends -------------------------------------------
# Select original inshore waypoints
tx.i <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    grp = "original",
    loc = "inshore",
    order = 2)

# Get N and S inshore waypoints ------------------------------------------------
# Calculate inshore/north transects
spacing = 10

tx.i.n <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    lat  = swfscMisc::destination(lat, long, brg + 90, spacing/2, units = "nm")["lat"],
    long = swfscMisc::destination(lat, long, brg + 90, spacing/2, units = "nm")["lon"],
    grp = "north",
    loc = "inshore",
    order = 1)

# Calculate inshore/south transects
tx.i.s <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    lat  = swfscMisc::destination(lat, long, brg - 90, spacing/2, units = "nm")["lat"],
    long = swfscMisc::destination(lat, long, brg - 90, spacing/2, units = "nm")["lon"],
    grp = "south",
    loc = "inshore",
    order = 3)

# Combine all inshore transects
tx.i <- tx.i %>%
  bind_rows(tx.i.n) %>%
  bind_rows(tx.i.s) %>%
  arrange(transect, desc(order))

# Get original offshore transect ends ------------------------------------------
tx.o <- tx.ends %>% 
  select(-lat.i, -long.i) %>% 
  rename(lat = lat.o, long = long.o) %>% 
  mutate(
    grp = "original",
    loc = "offshore",
    order = 2)

# Get N and S offshore waypoints -----------------------------------------------
# Calculate offshore/north transects
tx.o.n <- tx.ends %>% 
  select(-lat.i, -long.i) %>% 
  rename(lat = lat.o, long = long.o) %>% 
  mutate(
    lat  = swfscMisc::destination(lat, long, brg + 90, spacing/2, units = "nm")["lat"],
    long = swfscMisc::destination(lat, long, brg + 90, spacing/2, units = "nm")["lon"],
    grp = "north",
    loc = "offshore",
    order = 3)

# Calculate offshore/south transects
tx.o.s <- tx.ends %>% 
  select(-lat.i, -long.i) %>% 
  rename(lat = lat.o, long = long.o) %>% 
  mutate(
    lat  = swfscMisc::destination(lat, long, brg - 90, spacing/2, units = "nm")["lat"],
    long = swfscMisc::destination(lat, long, brg - 90, spacing/2, units = "nm")["lon"],
    grp = "south",
    loc = "offshore",
    order = 1)

tx.o.final <- data.frame()

for (v in unique(tx.o$vessel.name)) {
  if (v == "SD1024") {
    tx.o.tmp <- filter(tx.o, vessel.name == v)
  } else {
    tx.o.tmp <- filter(tx.o, vessel.name == v) %>% 
      bind_rows(filter(tx.o.n, vessel.name == v)) %>% 
      bind_rows(filter(tx.o.s, vessel.name == v))
  }
  tx.o.final <- bind_rows(tx.o.final, tx.o.tmp)
}

tx.o <- tx.o.final %>% 
  arrange(desc(transect), desc(order))

# tx.o <- tx.o %>% 
#   bind_rows(tx.o.n) %>% 
#   bind_rows(tx.o.s) %>% 
#   arrange(desc(transect), desc(order))

# Assemble the final data frame with all waypoints -----------------------------
strata.points <- tx.i %>% 
  bind_rows(tx.o)   %>%
  mutate(key = paste(transect.name, grp)) 

# Convert to points
strata.points.sf <- st_as_sf(strata.points, coords = c("long","lat"), crs = crs.geog) 

# Read 5m bathymetry points shapefile
bathy_5m_points <- st_read(here("Data/GIS/isobath_5m_final.shp"))

# Read 5m bathymetry polygon
bathy_5m_poly <- bathy_5m_points %>% 
  summarise(do_union = FALSE) %>% 
  st_cast("POLYGON") %>% 
  st_make_valid()

# Read 20 m bathy polygon
bathy_20m_poly <- st_read(here("Data/GIS/bathy_20m_polygon.shp"))

# Create polygons
strata.super.polygons <- strata.points.sf %>%
  group_by(vessel.name) %>%
  summarise(do_union = FALSE) %>%
  st_cast("POLYGON") %>%
  st_make_valid() %>%
  st_difference(st_union(bathy_20m_poly)) %>%
  mutate(area = st_area(.))

# mapview(strata.super.polygons)

# Convert to lines
tx.lines.sf <- strata.points.sf %>% 
  group_by(transect, grp) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast("LINESTRING")

# Create polygons for clipping nearshore strata
# Get original inshore transect ends -------------------------------------------
# Select original inshore waypoints
tx.c <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    lat  = destination(lat, long, brg, -20, units = "nm")["lat"],
    long = destination(lat, long, brg, -20, units = "nm")["lon"],
    grp = "original",
    loc = "onshore",
    order = 2)

# Get N and S inshore waypoints ------------------------------------------------
# Calculate inshore/north transects
tx.c.n <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    lat1  = destination(lat, long, brg + 90, spacing/2, units = "nm")["lat"],
    long1 = destination(lat, long, brg + 90, spacing/2, units = "nm")["lon"],
    lat2  = destination(lat1, long, brg, -20, units = "nm")["lat"],
    long2 = destination(lat1, long, brg, -20, units = "nm")["lon"],
    grp = "north",
    loc = "onshore",
    order = 1) %>% 
  select(-lat, -long, -lat1, -long1) %>% 
  rename(lat = lat2, long = long2)

tx.c.s <- tx.ends %>% 
  select(-lat.o, -long.o) %>% 
  rename(lat = lat.i, long = long.i) %>% 
  mutate(
    lat1  = destination(lat, long, brg - 90, spacing/2, units = "nm")["lat"],
    long1 = destination(lat, long, brg - 90, spacing/2, units = "nm")["lon"],
    lat2  = destination(lat1, long, brg, -20, units = "nm")["lat"],
    long2 = destination(lat1, long, brg, -20, units = "nm")["lon"],
    grp = "south",
    loc = "onshore",
    order = 3) %>% 
  select(-lat, -long, -lat1, -long1) %>% 
  rename(lat = lat2, long = long2)

# Combine all inshore transects
tx.c <- tx.c %>% 
  bind_rows(tx.c.n) %>% 
  bind_rows(tx.c.s) %>%   
  arrange(transect, desc(order))

# Create final stratum polygons
# Read 20 m isobath points; used to create nearshore stratum polygons
nearshore.points.sf <- st_read(here("Data/GIS/bathy_20m_points.shp")) 

nearshore.points <- nearshore.points.sf %>% 
  mutate(
    long = as.data.frame(st_coordinates(.))$X,
    lat = as.data.frame(st_coordinates(.))$Y,
    order = seq(1, n())) %>% 
  select(-Id) %>% 
  st_set_geometry(NULL)

# Ensure strata objects do not exist
if (exists("strata.primary")) rm(strata.primary)

# Create final strata and calculate area
# Create df for transect-level stock info
nasc.stock <- data.frame()

## Do so for primary strata using Lasker catch data
for (i in unique(strata.final$scientificName)) {
  # Select each strata per species
  strata.sub <- filter(strata.final, scientificName == i) %>% 
    select(transect, stratum)  
  
  # Define strata to stock
  nasc.stock.temp <- strata.points %>% 
    filter(loc == "inshore", grp == "original") %>%
    mutate(stock = case_when(
      i == "Engraulis mordax" & lat >= stock.break.anch ~ "Northern",
      i == "Engraulis mordax" & lat <  stock.break.anch ~ "Central",
      i == "Sardinops sagax"  & lat >= stock.break.sar  ~ "Northern",
      i == "Sardinops sagax"  & lat <  stock.break.sar  ~ "Southern",
      i %in% c("Clupea pallasii","Scomber japonicus",
               "Trachurus symmetricus","Etrumeus acuminatus") ~ "All"),
      scientificName = i) %>% 
    select(transect, stock, scientificName) %>% 
    distinct()
  
  # Combine results
  nasc.stock <- bind_rows(nasc.stock, nasc.stock.temp)
  
  for (j in sort(unique(strata.sub$stratum))) {
    # Create offshore stratum polygons ----------------------------------------
    # Add stratum numbers and stock designation to strata.points
    primary.poly.temp <- strata.points %>% 
      left_join(strata.sub) %>%
      left_join(nasc.stock.temp) %>% 
      filter(stratum == j) %>%
      mutate(scientificName = i) %>% 
      ungroup()
    
    # Select the southern-most inshore point for j-th stratum
    primary.poly.j.s <- primary.poly.temp %>%
      filter(loc == "inshore") %>% 
      slice(1)
    
    # Select the northern-most inshore point for j-th stratum
    primary.poly.j.n <- primary.poly.temp %>%
      filter(loc == "inshore") %>% 
      slice(n())
    
    # Select only the original inshore waypoints for j-th stratum
    primary.poly.j.i <- primary.poly.temp %>% 
      filter(loc == "inshore", grp == "original") %>% 
      mutate(scientificName = i)
    
    # Create the final polygon
    primary.poly.j <- primary.poly.temp %>% 
      filter(loc == "offshore") %>% 
      bind_rows(primary.poly.j.s) %>%
      bind_rows(primary.poly.j.i) %>%
      bind_rows(primary.poly.j.n) %>% 
      st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
      group_by(stratum, stock) %>% 
      summarise(do_union = FALSE) %>% 
      st_cast("POLYGON") %>% 
      mutate(scientificName = i)
    
    # Create nearshore stratum polygons ----------------------------------------
    # Select the inshore portions of offshore strata
    nearshore.temp <- primary.poly.j.s %>%
      bind_rows(primary.poly.j.i) %>%
      bind_rows(primary.poly.j.n) %>% 
      select(long, lat, scientificName, stratum, loc, stock)
    
    # Create polygon for clipping nearshore strata -----------------------------
    # Get the points on land
    clip.temp <- tx.c %>%
      left_join(strata.sub) %>%
      left_join(nasc.stock.temp) %>%
      filter(stratum == j) %>%
      ungroup()
    
    # Combine with offshore points from nearshore polygons
    clip.poly.j <- nearshore.temp %>%
      arrange(desc(lat)) %>%
      mutate(scientificName = i,
             stratum = j,
             loc = "nearshore",
             stock = unique(primary.poly.temp$stock)) %>% 
      bind_rows(clip.temp) %>% 
      st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
      group_by(stratum, stock) %>% 
      summarise(do_union = FALSE) %>% 
      st_cast("POLYGON") %>% 
      mutate(scientificName = i) %>% 
      st_make_valid() %>% 
      st_difference(bathy_5m_poly)
    
    # Combine with other polygons ----------------------------------------
    if (exists("strata.primary")) {
      strata.primary <- rbind(strata.primary, primary.poly.j)
    } else {
      strata.primary <- primary.poly.j
    }
  }
}

# ggplot(strata.primary) + geom_sf(aes(fill = stratum)) + facet_wrap(~scientificName)

# Save strata polygons
save(strata.primary, 
     file = here("Output/strata_primary_raw.Rdata"))

# Clip primary polygons using the 5 m isobathy polygon -------------------------
strata.primary <- strata.primary %>% 
  st_make_valid() %>% 
  st_difference(st_union(bathy_20m_poly)) %>% 
  st_difference(st_union(bathy_5m_poly)) %>% 
  ungroup() %>% 
  mutate(area = st_area(.)) #%>% left_join(select(strata.summ, scientificName, stratum, vessel.name))

# Ensure that strata are ordered by number
strata.primary <- arrange(strata.primary, scientificName, stratum)

# Save clipped primary polygons
save(strata.primary, 
     file = here("Output/strata_primary_final.Rdata"))

# Write final strata to shapefile
strata.primary.sub <- strata.primary %>% 
  filter(!sf::st_geometry_type(.) == "GEOMETRYCOLLECTION") %>% 
  mutate(area = round(as.numeric(area))) %>% 
  ungroup() 

# Write polygons to shapefile
strata.primary.sub %>% 
  select(-area) %>% 
  sf::st_write(here::here("Output/strata_primary.shp"),
               delete_layer = TRUE)

# Convert polygons to points and add coordinates -------------------------------
strata.primary.points  <- strata.primary %>% 
  st_cast("MULTIPOINT") %>%
  st_cast("POINT") %>%
  mutate(
    long = as.data.frame(st_coordinates(.))$X,
    lat = as.data.frame(st_coordinates(.))$Y,
    grp = paste(scientificName, stock, stratum)) %>% 
  st_set_geometry(NULL)


# Save final strata points
save(strata.primary.points,  
     file = here("Output/strata_points_primary.Rdata"))

# Write primary stata points to CSV
write.csv(strata.primary.points,  
          file = here("Output/strata_points_primary.csv"),
          quote = FALSE, row.names = FALSE)

# Summarize nasc.strata by stock
strata.summ.primary <- strata.primary %>% 
  select(scientificName, stratum, stock, area) %>%
  mutate(area = as.numeric(area)) %>% 
  st_set_geometry(NULL)


# Write stata summaries to CSV
write.csv(strata.summ.primary,  
          file = here("Output/strata_summary_primary.csv"),
          quote = FALSE, row.names = FALSE)

## Do so for primary strata using Lisa Marie purse-seine catch data
# Ensure strata objects do not exist
if (exists("strata.primary.ns")) rm(strata.primary.ns)

# Create final strata and calculate area
# Create df for transect-level stock info
nasc.stock.ns <- data.frame()


for (i in unique(strata.final.ns$scientificName)) {
  # Select each strata per species
  strata.sub <- filter(strata.final.ns, scientificName == i) %>% 
    select(transect, stratum)  
  
  # Define strata to stock
  nasc.stock.temp <- strata.points %>% 
    filter(loc == "inshore", grp == "original") %>%
    mutate(stock = case_when(
      i == "Engraulis mordax" & lat >= stock.break.anch ~ "Northern",
      i == "Engraulis mordax" & lat <  stock.break.anch ~ "Central",
      i == "Sardinops sagax"  & lat >= stock.break.sar  ~ "Northern",
      i == "Sardinops sagax"  & lat <  stock.break.sar  ~ "Southern",
      i %in% c("Clupea pallasii","Scomber japonicus",
               "Trachurus symmetricus","Etrumeus acuminatus") ~ "All"),
      scientificName = i) %>% 
    select(transect, stock, scientificName) %>% 
    distinct()
  
  # Combine results
  nasc.stock.ns <- bind_rows(nasc.stock.ns, nasc.stock.temp)
  
  for (j in sort(unique(strata.sub$stratum))) {
    # Create offshore stratum polygons ----------------------------------------
    # Add stratum numbers and stock designation to strata.points
    primary.poly.temp <- strata.points %>% 
      left_join(strata.sub) %>%
      left_join(nasc.stock.temp) %>% 
      filter(stratum == j) %>%
      mutate(scientificName = i) %>% 
      ungroup()
    
    # Select the southern-most inshore point for j-th stratum
    primary.poly.j.s <- primary.poly.temp %>%
      filter(loc == "inshore") %>% 
      slice(1)
    
    # Select the northern-most inshore point for j-th stratum
    primary.poly.j.n <- primary.poly.temp %>%
      filter(loc == "inshore") %>% 
      slice(n())
    
    # Select only the original inshore waypoints for j-th stratum
    primary.poly.j.i <- primary.poly.temp %>% 
      filter(loc == "inshore", grp == "original") %>% 
      mutate(scientificName = i)
    
    # Create the final polygon
    primary.poly.j <- primary.poly.temp %>% 
      filter(loc == "offshore") %>% 
      bind_rows(primary.poly.j.s) %>%
      bind_rows(primary.poly.j.i) %>%
      bind_rows(primary.poly.j.n) %>% 
      st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
      group_by(stratum, stock) %>% 
      summarise(do_union = FALSE) %>% 
      st_cast("POLYGON") %>% 
      mutate(scientificName = i)
    
    # Create nearshore stratum polygons ----------------------------------------
    # Select the inshore portions of offshore strata
    nearshore.temp <- primary.poly.j.s %>%
      bind_rows(primary.poly.j.i) %>%
      bind_rows(primary.poly.j.n) %>% 
      select(long, lat, scientificName, stratum, loc, stock)
    
    # Create polygon for clipping nearshore strata -----------------------------
    # Get the points on land
    clip.temp <- tx.c %>%
      left_join(strata.sub) %>%
      left_join(nasc.stock.temp) %>%
      filter(stratum == j) %>%
      ungroup()
    
    # Combine with offshore points from nearshore polygons
    clip.poly.j <- nearshore.temp %>%
      arrange(desc(lat)) %>%
      mutate(scientificName = i,
             stratum = j,
             loc = "nearshore",
             stock = unique(primary.poly.temp$stock)) %>% 
      bind_rows(clip.temp) %>% 
      st_as_sf(coords = c("long","lat"), crs = crs.geog) %>% 
      group_by(stratum, stock) %>% 
      summarise(do_union = FALSE) %>% 
      st_cast("POLYGON") %>% 
      mutate(scientificName = i) %>% 
      st_make_valid() %>% 
      st_difference(bathy_5m_poly)
    
    # Combine with other polygons ----------------------------------------
    if (exists("strata.primary.ns")) {
      strata.primary.ns <- rbind(strata.primary.ns, primary.poly.j)
    } else {
      strata.primary.ns <- primary.poly.j
    }
  }
}

# ggplot(strata.primary.ns) + geom_sf(aes(fill = factor(stratum))) + facet_wrap(~scientificName)

# Save strata polygons
save(strata.primary.ns, 
     file = here("Output/strata_primary_ns_raw.Rdata"))

# Clip primary polygons using the 5 m isobathy polygon -------------------------
strata.primary.ns <- strata.primary.ns %>% 
  st_make_valid() %>% 
  st_difference(st_union(bathy_20m_poly)) %>% 
  st_difference(st_union(bathy_5m_poly)) %>% 
  ungroup() %>% 
  mutate(area = st_area(.)) #%>% left_join(select(strata.summ, scientificName, stratum, vessel.name))

# Ensure that strata are ordered by number
strata.primary.ns <- arrange(strata.primary.ns, scientificName, stratum)

# Save clipped primary polygons
save(strata.primary.ns, 
     file = here("Output/strata_primary_ns_final.Rdata"))

# Write final strata to shapefile
strata.primary.ns.sub <- strata.primary.ns %>% 
  filter(!sf::st_geometry_type(.) == "GEOMETRYCOLLECTION") %>% 
  mutate(area = round(as.numeric(area))) %>% 
  ungroup() 

# Write polygons to shapefile
strata.primary.ns.sub %>% 
  select(-area) %>% 
  sf::st_write(here::here("Output/strata_primary_ns.shp"),
               delete_layer = TRUE)

# Convert polygons to points and add coordinates -------------------------------
strata.primary.ns.points  <- strata.primary.ns %>% 
  st_cast("MULTIPOINT") %>%
  st_cast("POINT") %>%
  mutate(
    long = as.data.frame(st_coordinates(.))$X,
    lat = as.data.frame(st_coordinates(.))$Y,
    grp = paste(scientificName, stock, stratum)) %>% 
  st_set_geometry(NULL)


# Save final strata points
save(strata.primary.ns.points,  
     file = here("Output/strata_points_primary_ns.Rdata"))

# Write primary stata points to CSV
write.csv(strata.primary.ns.points,  
          file = here("Output/strata_points_primary_ns.csv"),
          quote = FALSE, row.names = FALSE)

# Summarize nasc.strata by stock
strata.summ.primary.ns <- strata.primary.ns %>% 
  select(scientificName, stratum, stock, area) %>%
  mutate(area = as.numeric(area)) %>% 
  st_set_geometry(NULL)


# Write stata summaries to CSV
write.csv(strata.summ.primary.ns,  
          file = here("Output/strata_summary_primary_ns.csv"),
          quote = FALSE, row.names = FALSE)

# Compute point estimates
# Save final nasc data frame used for point and bootstrap estimates
save(nasc, file = here("Output/nasc_final.Rdata"))
write_csv(nasc, file = here("Output/nasc_final.csv"))

if (exists("nasc.summ.strata")) rm(nasc.summ.strata)
if (exists("point.estimates")) rm(point.estimates)

# Calculate point estimates for each species
for (i in unique(strata.final$scientificName)) {
  # Subset strata for species i
  strata.temp <- filter(strata.final, scientificName == i) %>% 
    select(transect, stratum)
  
  # Add stratum numbers to nasc
  nasc.temp <- nasc.core.final %>%
    left_join(strata.temp) %>% 
    filter(!is.na(stratum))
  
  # ggplot() +
  #   geom_sf(data = filter(strata.primary, scientificName == i)) +
  #   geom_point(data = nasc.temp, aes(long, lat, colour = stratum))
  
  # Summarise nasc by stratum
  nasc.temp.summ <- nasc.temp %>% 
    group_by(stratum) %>% 
    summarise(
      n_samples = n(),
      mean_nasc = mean(cps.nasc)) %>% 
    mutate(scientificName = i) %>% 
    select(scientificName, everything())
  
  # Combine nasc summaries
  if (exists("nasc.summ.strata")) {
    nasc.summ.strata <- bind_rows(nasc.summ.strata, 
                                  nasc.temp.summ)
  } else {
    nasc.summ.strata <- nasc.temp.summ
  }
  
  # Create data frame with stratum and area (m^2)  
  stratum.info <- strata.primary %>%
    filter(scientificName == i) %>%
    select(stratum, area) %>%
    mutate(area = as.numeric(area)) %>%
    st_set_geometry(NULL) 
  
  # Compute point estimates
  # Currently has na.rm = TRUE for calculating biomass
  if (exists("point.estimates")) {
    point.estimates <- bind_rows(point.estimates,
                                 data.frame(scientificName = i,
                                            estimate_point(nasc.temp, stratum.info, species = i)))
  } else {
    point.estimates <- data.frame(scientificName = i,
                                  estimate_point(nasc.temp, stratum.info, species = i))
  }
}

# Save results
save(point.estimates, 
     file = here("Output/biomass_point_estimates.Rdata")) 

# Save strata nasc summaries to CSV
write_csv(nasc.summ.strata, here("Output/nasc_strata_summary.csv"))

# Add stock designations to point estimates
point.estimates     <- left_join(point.estimates, strata.summ.primary)

# Summarize point estimates (by stocks)
pe <- point.estimates %>%
  group_by(scientificName, stock) %>%
  summarise(
    area          = sum(area),
    biomass.total = sum(biomass.total)) %>%
  bind_rows(point.estimates) %>%
  mutate(area = area * 2.915533e-07) %>%
  arrange(scientificName, stock, stratum) %>%
  mutate(stratum = case_when(
    is.na(stratum) ~ "All",
    TRUE ~  as.character(stratum))) %>% 
  select(scientificName, stock, stratum, area, biomass.total) %>%
  rename(
    Species            = scientificName,
    Stock              = stock,
    Stratum            = stratum,
    Area               = area,
    biomass.mean.point = biomass.total)

# Save point estimates
save(pe, file = here("Output/biomass_point_estimates_final.Rdata"))
write_csv(pe, here("Output/biomass_point_estimates_final.csv"))

if (exists("nasc.summ.strata.ns")) rm(nasc.summ.strata.ns)
if (exists("point.estimates.ns")) rm(point.estimates.ns)

# Calculate point estimates for each species
for (i in unique(strata.final.ns$scientificName)) {
  # Subset strata for species i
  strata.temp <- filter(strata.final.ns, scientificName == i) %>% 
    select(transect, stratum)
  
  # Add stratum numbers to nasc
  nasc.temp <- nasc.ns.final %>%
    left_join(strata.temp) %>% 
    filter(!is.na(stratum))
  
  # ggplot() +
  #   geom_sf(data = filter(strata.primary.ns, scientificName == i)) +
  #   geom_point(data = nasc.temp, aes(long, lat, colour = stratum))
  
  # Summarise nasc by stratum
  nasc.temp.summ <- nasc.temp %>% 
    group_by(stratum) %>% 
    summarise(
      n_samples = n(),
      mean_nasc = mean(cps.nasc)) %>% 
    mutate(scientificName = i) %>% 
    select(scientificName, everything())
  
  # Combine nasc summaries
  if (exists("nasc.summ.strata.ns")) {
    nasc.summ.strata.ns <- bind_rows(nasc.summ.strata.ns, 
                                  nasc.temp.summ)
  } else {
    nasc.summ.strata.ns <- nasc.temp.summ
  }
  
  # Create data frame with stratum and area (m^2)  
  stratum.info <- strata.primary.ns %>%
    filter(scientificName == i) %>%
    select(stratum, area) %>%
    mutate(area = as.numeric(area)) %>%
    st_set_geometry(NULL) 
  
  # Compute point estimates
  # Currently has na.rm = TRUE for calculating biomass
  if (exists("point.estimates.ns")) {
    point.estimates.ns <- bind_rows(point.estimates.ns,
                                 data.frame(scientificName = i,
                                            estimate_point(nasc.temp, stratum.info, species = i)))
  } else {
    point.estimates.ns <- data.frame(scientificName = i,
                                  estimate_point(nasc.temp, stratum.info, species = i))
  }
}

# Save results
save(point.estimates.ns, 
     file = here("Output/biomass_point_estimates_ns.Rdata")) 

# Save strata nasc summaries to CSV
write_csv(nasc.summ.strata.ns, here("Output/nasc_strata_summary_ns.csv"))

# Add stock designations to point estimates
point.estimates.ns     <- left_join(point.estimates.ns, strata.summ.primary.ns)

# Summarize point estimates (by stocks)
pe.ns <- point.estimates.ns %>%
  group_by(scientificName, stock) %>%
  summarise(
    area          = sum(area),
    biomass.total = sum(biomass.total)) %>%
  bind_rows(point.estimates.ns) %>%
  mutate(area = area * 2.915533e-07) %>%
  arrange(scientificName, stock, stratum) %>%
  mutate(stratum = case_when(
    is.na(stratum) ~ "All",
    TRUE ~  as.character(stratum))) %>% 
  select(scientificName, stock, stratum, area, biomass.total) %>%
  rename(
    Species            = scientificName,
    Stock              = stock,
    Stratum            = stratum,
    Area               = area,
    biomass.mean.point = biomass.total)

# Save point estimates
save(pe.ns, file = here("Output/biomass_point_estimates_final_ns.Rdata"))
write_csv(pe.ns, here("Output/biomass_point_estimates_final_ns.csv"))

# Combine results from both methods
pe.combo <- pe %>% 
  left_join(select(pe.ns, Species, Stock, Stratum, Area.ns = Area, biomass.mean.point.ns = biomass.mean.point))
pe.combo.summ <- filter(pe.combo, Stratum == "All") %>% 
  mutate(pct.diff = (biomass.mean.point - biomass.mean.point.ns)/biomass.mean.point*100)
write_csv(pe.combo, file = here("Output/point_estimates_combo.csv"))
write_csv(pe.combo.summ, file = here("Output/point_estimates_combo_summ.csv"))

# Plot biomass density, stratum boundaries, and catch proportions by stock and net sample type
nasc.density.core <- project_df(ungroup(nasc.density.core), to = crs.proj)
nasc.density.ns <- project_df(ungroup(nasc.density.ns), to = crs.proj)

# Select legend objects 
dens.levels.all <- unique(sort(c(unique(nasc.density.core$bin.level),unique(nasc.density.ns$bin.level))))
dens.labels.all <- dens.labels[dens.levels.all]
dens.sizes.all  <- dens.sizes[dens.levels.all]
dens.colors.all <- dens.colors[dens.levels.all]

# Plot trawl sampling data
dens.plot.final <- base.map +
  geom_sf(data = st_transform(strata.primary, 3310),
          aes(colour = factor(stratum)), fill = NA, size = 1) +
  scale_colour_discrete('Stratum') + 
  # Plot zero nasc data
  geom_point(data = filter(nasc, cps.nasc == 0), aes(X, Y),
             colour = 'gray50',size = 0.15, alpha = 0.5) +
  # Plot NASC data
  geom_point(data = filter(nasc.density.core, density > 0), aes(X, Y, size = bin, fill = bin),
             shape = 21, alpha = 0.75) +
  # Configure size and colour scales
  scale_size_manual(name = bquote(atop(Biomass~density, ~'(t'~'nmi'^-2*')')),
                    values = dens.sizes.all, labels = dens.labels.all) +
  scale_fill_manual(name = bquote(atop(Biomass~density, ~'(t'~'nmi'^-2*')')),
                    values = dens.colors.all, labels = dens.labels.all) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"]))) + 
  facet_wrap(~scientificName, nrow = 1) + 
  ggtitle("Biomass density - Trawl catch") +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "italic"))

# Plot seine sampling data
dens.plot.final.ns <- base.map +
  geom_sf(data = st_transform(strata.primary.ns, 3310),
          aes(colour = factor(stratum)), fill = NA, size = 1) +
  scale_colour_discrete('Stratum') + 
  # Plot zero nasc data
  geom_point(data = filter(nasc, cps.nasc == 0), aes(X, Y),
             colour = 'gray50',size = 0.15, alpha = 0.5) +
  # Plot NASC data
  geom_point(data = filter(nasc.density.ns, density > 0), aes(X, Y, size = bin, fill = bin),
             shape = 21, alpha = 0.75) +
  # Configure size and colour scales
  scale_size_manual(name = bquote(atop(Biomass~density, ~'(t'~'nmi'^-2*')')),
                    values = dens.sizes.all, labels = dens.labels.all) +
  scale_fill_manual(name = bquote(atop(Biomass~density, ~'(t'~'nmi'^-2*')')),
                    values = dens.colors.all, labels = dens.labels.all) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"]))) + 
  facet_wrap(~scientificName, nrow = 1) + 
  ggtitle("Biomass density - Purse-seine catch") +
  theme(strip.background.x = element_blank(),
        strip.text.x       = element_text(face = "italic"))

# Combine plots
dens.plot.final.combo <- dens.plot.final / dens.plot.final.ns

# Save plot
ggsave(dens.plot.final.combo, file = here("Figs/fig_biomass_dens_combo.png"),
       width = 14, height = 12)

# Summarize NASC for plotting
nasc.plot.core <- nasc.core.final %>%
  select(transect, int, lat, long, cps.nasc) %>% 
  group_by(transect, int) %>% 
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

nasc.plot.ns <- nasc.ns.final %>%
  select(transect, int, lat, long, cps.nasc) %>% 
  group_by(transect, int) %>% 
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
nasc.levels.all <- sort(unique(nasc.plot.core$bin.level))
nasc.labels.all <- nasc.labels[sort(nasc.levels.all)]
nasc.sizes.all  <- nasc.sizes[sort(nasc.levels.all)]
nasc.colors.all <- nasc.colors[sort(nasc.levels.all)]
  
# Plot averaged cps NASC with nearest neighbor polygons
plot.core.clusters.nasc <- base.map +
  geom_point(data = nasc.plot.core, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75) +
  # Configure size and colour scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all, labels = nasc.labels.all) +
  geom_sf(data = nasc.super.clusters, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_scatterpie(data = clf, 
                  aes(X, Y, group = cluster, r = 15000),
                  cols = c("prop.anch","prop.her","prop.jack",
                           "prop.mack", "prop.sar"),
                  alpha = 0.5) +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "P. herring", "J. mackerel",
                               "P. mackerel", "Sardine"),
                    values = c(anchovy.color, pac.herring.color, jack.mack.color,  
                               pac.mack.color, sardine.color)) +
  # # Configure legend guides
  # guides(fill   = guide_legend(order = 1), 
  #        size   = guide_legend(order = 1)) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))

plot.ns.clusters.nasc <- base.map +
  geom_point(data = nasc.plot.ns, aes(X, Y, size = bin, fill = bin), 
             shape = 21, alpha = 0.75, show.legend = FALSE) +
  # Configure size and colour scales
  scale_size_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.sizes.all,labels = nasc.labels.all) +
  scale_fill_manual(name = bquote(atop(italic(s)[A], ~'(m'^2 ~'nmi'^-2*')')),
                    values = nasc.colors.all, labels = nasc.labels.all) +
  geom_sf(data = nasc.super.clusters.ns, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_scatterpie(data = clf.ns, 
                  aes(X, Y, group = cluster, r = 15000),
                  cols = c("prop.anch","prop.her","prop.jack",
                           "prop.sar"),
                  alpha = 0.5, show.legend = FALSE) +
  # Configure trawl scale
  scale_fill_manual(name = 'Species',
                    labels = c("Anchovy", "P. herring", "J. mackerel",
                               "Sardine"),
                    values = c(anchovy.color, pac.herring.color, jack.mack.color,  
                               sardine.color)) +
  coord_sf(crs = crs.proj, # CA Albers Equal Area Projection
           xlim = unname(c(map.bounds["xmin"]-100000, map.bounds["xmax"]+100000)),  
           ylim = unname(c(map.bounds["ymin"], map.bounds["ymax"])))


nasc.plot.combo <- plot.core.clusters.nasc + plot.ns.clusters.nasc +
  plot_annotation(tag_levels = c('a'), tag_suffix = ')')

ggsave(nasc.plot.combo, filename = here("figs/fig_nasc_cluster_combo.png"),
       height = 7, width = 8)

