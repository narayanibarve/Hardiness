
# PLEASE CHANGE THIS AS NEEDED. NO NEED TO KEEP MY NAME NOR ANY OF THE COMMENTS

# Project: Hardiness zones
# Date: 23-08-2022 
# Marlon E. Cobos

# Setting up -------------------------------------------------------------------
# packages
library(raster)  
library(prism)   
library(spThin)  
library(rgdal)

#install.packages("remotes")
#remotes::install_github("fmachados/grinnell")
library(grinnell)
library(ellipsenm)

# set directory to save results
set("YOUR/DIRECTORY")  # change this and place inside the data
# ------------------------------------------------------------------------------



# Obtaining environmental data and pre-processing records ----------------------
# environmental data (obtaining normals for 30 years)
# set a directory fro prism 
prism_set_dl_dir("prism_data")

# getting the data
get_prism_normals(type = "tmax", resolution = "4km", annual = TRUE)
get_prism_normals(type = "tmin", resolution = "4km", annual = TRUE)
get_prism_normals(type = "tmean", resolution = "4km", annual = TRUE)
get_prism_normals(type = "ppt", resolution = "4km", annual = TRUE)
get_prism_normals(type = "vpdmax", resolution = "4km", annual = TRUE)
get_prism_normals(type = "vpdmin", resolution = "4km", annual = TRUE)

# locating raster data in groups by month
prism_files <- pd_to_file(list.files("prism_data", pattern = "_bil$"))
prism <- stack(prism_files)

# rewrite layers for simulation
dir.create("simulation_var")

nam <- names(prism)
nam <- gsub("_30.*", "", nam)
nam <- paste0("simulation_var/", nam, ".tif")

re <- lapply(1:nlayers(prism), function(x) {
  writeRaster(prism[[x]], filename = nam[x], format = "GTiff")
})
# ------------------------------------------------------------------------------


# Occurrence data thinning -----------------------------------------------------
# occurrence data
spnames <- c("Abies_balsamea", "Carnegiea_gigantea", "Juglans_nigra",
             "Magnolia_grandiflora")

# thinning in loop 
occ_thin <- lapply(spnames, function(x) {
  ## read spatial file
  spps <- paste0(x, "_reduced")
  occ <- readOGR("data", layer = spps)
  
  ## get coordinates
  occr <- coordinates(occ)
  occd <- data.frame(species = x, longitude = occr[, 1], latitude = occr[, 2])
  
  ## thinned result file names
  octhin <- paste0(x, "_thinned")
  
  ## spatial thinning
  occt <- thin(loc.data = occd, lat.col = "latitude", long.col = "longitude",   
               spec.col = "species", thin.par = 16, reps = 5, 
               locs.thinned.list.return = FALSE, write.files = TRUE,
               max.files = 1, out.dir = ".", out.base = octhin)  
})
# ------------------------------------------------------------------------------



# M simulations ----------------------------------------------------------------
# points
spps <- list(A_b = read.csv("Abies_balsamea_thinned.csv"),
             C_g = read.csv("Carnegiea_gigantea_thinned.csv"),
             J_n = read.csv("Juglans_nigra_thinned.csv"),
             M_g = read.csv("Magnolia_grandiflora_thinned.csv"))

# variables
lvar <- list.files("simulation_var", pattern = ".tif$", full.names = TRUE)

all_vars <- stack(lvar)

# species
sspnames <- names(spps)

# standard deviations to be tested
sds <- 3:5

# simulations in loop
all_ms <- lapply(1:length(spps), function(w) {
  spname <- sspnames[w]
  
  ## thinning points
  datred <- thin_data(data = spps[[spname]], longitude = "Longitude", 
                      latitude = "Latitude", thin_distance = 20)
  
  ## SD 1-5, kernel = normal, 125 dispersal events
  res_folder125 <- paste0(spname, "_125de_M_normal_SD", sds)
  
  M_species_125de <- lapply(1:length(sds), function(x) {
    M_simulationR(data = datred,
                  current_variables = all_vars, 
                  dispersal_kernel = "normal", 
                  kernel_spread = sds[x], max_dispersers = 4, 
                  dispersal_events = 125, access_threshold = 5, 
                  out_format = "GTiff", 
                  output_directory = res_folder125[x])
  })
  
  ## SD 1-5, kernel = normal, 250 dispersal events
  res_folder250 <- paste0(spname, "_250de_M_normal_SD", sds)
  
  M_species_125de <- lapply(1:length(sds), function(x) {
    M_simulationR(data = datred, 
                  current_variables = all_vars, 
                  dispersal_kernel = "normal", 
                  kernel_spread = sds[x], max_dispersers = 4, 
                  dispersal_events = 250, access_threshold = 5, 
                  out_format = "GTiff", 
                  output_directory = res_folder250[x])
  })
})


# get Ms
ab <- rast("A_b_250de_M_normal_SD5/accessible_area_M.tif")
cg <- rast("C_g_250de_M_normal_SD5/accessible_area_M.tif")
jn <- rast("J_n_250de_M_normal_SD5/accessible_area_M.tif")
mg <- rast("M_g_250de_M_normal_SD5/accessible_area_M.tif")

# make them polygons
ab <- as.polygons(ab)
cg <- as.polygons(cg)
jn <- as.polygons(jn)
mg <- as.polygons(mg)

# write Ms back to directory
dir.create("Final_M")

writeVector(ab, filename = "Final_M/A_b_M.shp")
writeVector(cg, filename = "Final_M/C_g_M.shp")
writeVector(jn, filename = "Final_M/J_n_M.shp")
writeVector(mg, filename = "Final_M/M_g_M.shp")
# ------------------------------------------------------------------------------
