# Processing controls ----------------------------------------------------
## Settings in this section control various behaviors and tasks used in the main data processing scripts
import.data    <- T # Process purse seine data, if present
save.figs      <- T # Process near backscatter data; typically TRUE

# Survey plan info --------------------------------------------------------
wpt.types            <- c(Adaptive = "Adaptive", Compulsory = "Compulsory", 
                          Nearshore = "Nearshore", #Offshore = "Offshore",
                          Saildrone = "Saildrone")
wpt.colors           <- c(Adaptive = "#FF0000", Compulsory = "#000000",  
                          Nearshore = "#FF33F5", #Offshore = "#FFA500",
                          Saildrone = "#00FFFF") 
wpt.linetypes        <- c(Adaptive = "dashed", Compulsory = "solid",
                          Nearshore = "solid", #Offshore = "dashed", 
                          Saildrone = "solid")

# Growth model parameters ------------------------------------------------------
model.season  <- "summer" # spring or summer; for selecting growth model parameters
model.type    <- "glm"    # lm, nlm, or glm; for selecting growth model

# Mapping preferences -----------------------------------------------------
# Turn off S2 processing in sf
sf::sf_use_s2(FALSE)

# Coordinate reference systems for geographic and projected data
crs.geog <- 4326 # WGS84
crs.proj <- 3310 # California Albers Equal Area

# Map landmarks
label.list <- c("Monterey Bay","San Francisco","Cape Flattery","Crescent City",
                "Newport","Point Conception","Cape Mendocino","Columbia River",
                "Cape Blanco","Bodega Bay","Westport","Fort Bragg",
                "Morro Bay","Long Beach","Cape Scott","San Diego",
                "Ensenada","Punta Eugenia","El Rosario","Cabo San Lucas",
                "Punta Abreojos","San Carlos")

# Figure preferences ------------------------------------------------------
# Define species to be analysed
cps.spp            <- c("Clupea pallasii","Engraulis mordax","Sardinops sagax",
                        "Scomber japonicus","Trachurus symmetricus", 
                        "Etrumeus acuminatus")

# Set species colors
sardine.color      <- '#FF0000'
anchovy.color      <- '#00CD66'
jack.mack.color    <- '#0000FF'
jacksmelt.color    <- '#A020F0'
pac.mack.color     <- '#00FFFF'
pac.herring.color  <- '#F5DEB3'
rnd.herring.color  <- '#F0B81D'
other.color        <- 'gray'

# Set gear type colors
seine.color <- "white"
trawl.color <- "black"
  
# Trawl proportion plots
scale.pies <- FALSE   # Scale pie charts (TRUE/FALSE)
pie.scale  <- 0.0125  # 0.01-0.02 works well for coast-wide survey (i.e., summer), larger values (~0.03) for spring

# Trawl
# For legend objects
trawl.breaks       <- c(0, 1, 10, 25, 50, 500, 1000, 10000) 
trawl.labels       <- c("<1", "1-10", "10-25", "25-50", "50-500", "500-1000", ">1000") 
trawl.sizes        <- c(1, 2, 3, 4, 5, 6, 7) 

# NASC
# For legend objects
nasc.breaks        <- c(0, 1, 200, 500, 2000, 5000, 20000, 50000, 20000000)
nasc.labels        <- c("0","1-200", "200-500", "500-2000", "2000-5000", 
                        "5000-20,000", "20,000-50,000", ">50,000")
nasc.scale         <- 0.55 # Scale percentage (smaller for larger scale)
nasc.sizes         <- c(0.1, 0.25, 2, 3, 4, 5, 6, 7)*nasc.scale
nasc.colors        <- c("#000000", "#C2E6F2", "#1E90FF", "#FFFF00", "#FF8C00", 
                        "#FF0000", "#FFC0CB", "#FFFFFF")

# Acoustic biomass density map
dens.breaks        <- c(0, 1, 10, 100, 500, 1000, 10000, 50000, 1000000)
dens.labels        <- c("0-1", "1-10", "10-100", "100-500", "500-1000",
                        "1000-10,000", "10,000-50,000", ">50,000")
dens.colors        <- c("#000000", "#1E90FF", "#FFFF00", "#FF8C00", 
                        "#FF0000", "#FFC0CB", "#FFFFFF", "#00FF00") # for legend colors
dens.sizes         <- c(0.25, 1, 2.25, 3, 4.25, 5.5, 6.5, 7.5) # for legend sizes

# Catch map
# For legend objects
catch.breaks       <- c(0, 10, 100, 500, 1000)
catch.labels       <- c("0-10", "10-100", "100-500", "500-1000")
catch.pie.sizes    <- c(1, 2, 3, 4, 5, 6)

# Cluster relative length frequency
# Set number of columns in facet plot
lf.ncols <- 5

# Data sources ------------------------------------------------------------
# Backscatter data info
# Survey vessels that collected acoustic data (a character vector of vessel abbreviations)
nasc.vessels           <- c("RL","LM","LBC") #c("RL","LBC","LM","SD") 
nasc.vessels.nearshore <- c("LBC", "LM")

# Define columns to use for a fixed integration depth (if cps.nasc is not present)
# Options include 0-100 (by 5), 100, 150, 250, and 350 m.
# Defined by the atm::extract_csv() function.
nasc.depth.cps   <- "NASC.250"

# Purse seine data info
# Survey vessels that collected purse seine data
seine.vessels          <- c("LBC","LM")
# Use seine data to apportion backscatter
use.seine.data         <- TRUE

# Interval length (m); from Echoview
nasc.interval          <-  100    

# Number of intervals over which to summarize NASC
nasc.summ.interval     <- 2000/nasc.interval 

# Stock boundaries --------------------------------------------------------
stock.break.anch <- c("Cape Mendocino" = 40.50)  # Latitude of Cape Mendocino
stock.break.sar  <- c("Big Sur" = 36.153) # Latitude of Morro Bay (based on revised sardine model)
# stock.break.sar  <- 37.674 # Latitude of San Francisco, based on differences in length dist.
# stock.break.sar  <- c("Pt. Conception" = 34.46) # Latitude of Pt. Conception (or change based on SST)