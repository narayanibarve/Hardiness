#####
# Initial cleaning of occurrence data
## Excluding: records with no coordinates, duplicates, records with 
## corrdinates (0, 0)
#####

# data needed 
## (1) Species occurrences
## This file must contain the following columns (in that order): 
## ID, Species_name, Longitud, Latitud.

# defining working directory
setwd("C:/Users/c361n270/Desktop/Hardiness_zones") 
# change this to your working directory

#Reading species occurrences
D <- list.files(path = "Species_occurrences/GBIF", pattern = ".csv$", 
                full.names = T)
nam <- list.files(path = "Species_occurrences/GBIF", pattern = ".csv$", 
                            full.names = F)
Finalnam <- paste0("Species_occurrences/Final/", nam)

for (i in 1:length(D)){
  occurrences <- read.csv(D[i]) # occurrences 
  
  occurrences <- occurrences[, c("acceptedScientificName", "decimalLongitude", "decimalLatitude")]
  colnames(occurrences) <- c("Species", "Longitude", "Latitude")
  
  # Excluding duplicates
  occurrences$code <-  paste(occurrences$Species, occurrences$Longitude, # concatenating columns of interest
                             occurrences$Latitude, sep = "_")
  
  occurrences <- occurrences[!duplicated(occurrences$code), 1:4] # erasing duplicates
  occurrences <- na.omit(occurrences[, 1:3])
  
  # Excluding records with (0, 0) coordinates
  occurrences <- occurrences[occurrences$Longitude != 0 & occurrences$Latitude != 0, ]
  
  # Excluding recors with low level of precision (<= 2 decimals)
  ## samll function to detect precision 
  ## (from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r)
  # decimalplaces <- function(x) {
  #   if (abs(x - round(x)) > .Machine$double.eps^0.5) {
  #     nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  #   } else {
  #     return(0)
  #   }
  # }
  # 
  # occurrences <- occurrences[sapply(occurrences$Longitude, decimalplaces) >= 2 & # keep only the ones with more than 1 decimals
  #                              sapply(occurrences$Latitude, decimalplaces) >= 2, ]
  
  # saving the new set of occurrences inside continents and area of interest
  write.csv(occurrences, Finalnam[i], row.names = FALSE)
}

##--------------------------------------------------------------------------
#Keeping only the occurrences from the USA

library(rgdal)

Finalnam1 <- paste0("Species_occurrences/Final/USA_only/", nam)

usa <- readOGR(dsn = ".", layer = "STATES")
usa <- usa[!usa$STATE_NAME %in% c("Hawaii", "Alaska"), ]

D <- list.files(path = "Species_occurrences/Final", pattern = ".csv$", 
                full.names = T)
nam <- list.files(path = "Species_occurrences/Final", pattern = ".csv$", 
                  full.names = F)

for (i in 1:length(nam)){
  occ <- read.csv(D[i])
  occ1 <- occ[, c("Species", "Longitude", "Latitude")]
  
  occ1 <- SpatialPointsDataFrame(coords = occ1[, 2:3], data = occ,
                                 proj4string = usa@proj4string, match.ID = F)
  
  occ1 <- occ1[usa, ]
  
  if (nrow(occ1@data) > 0){
    write.csv(occ1@data, Finalnam1[i], row.names = F)
  }
}
