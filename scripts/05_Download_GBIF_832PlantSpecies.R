####
# PROJECT: Hardiness Zones 
# Performs: Downloading GBIF data of US species trees
# Date: 09-06-2022

# package
library(rgbif)

# defining working directory
## project folder
setwd("C:/Users/c361n270/Desktop/Hardiness_zones") # Your folder

spvector <- as.character(read.csv("Tree_species.csv", header = T)[, 1]) # binomial names

## folder for occurrences 
dir.create("Species_occurrences")
setwd("C:/Users/c361n270/Desktop/Hardiness_zones/Species_occurrences/GBIF") 

## Getting info species by species 
occ_count <- list() # object to save info on number of georeferenced records per species 

for (i in 1:length(spvector)) {
  sps <- try(name_lookup(query = spvector[i], rank = "species", 
                         return = "data", limit = 100), silent = TRUE) # information about the species
  
  sps_class <- class(sps)
  
  # avoiding errors from GBIF (e.g., species name not in GBIF)
  if(sps_class[1] == "try-error") {
    occ_count[[i]] <- c(Species = spvector[i], keys = 0, counts = 0) # species not in GBIF
    cat("species", spvector[i], "is not in the GBIF database\n")
    
  }else {
    keys <- sps$data$key # all keys returned
    counts <- vector() # object to save info on number of records per key
    
    for (j in 1:length(keys)) { # testing if keys return records
      counts[j] <- occ_count(taxonKey = keys[j], georeferenced = TRUE) 
    }
    
    if (sum(counts) == 0) { # if no info, tell the species
      occ_count[[i]] <- c(Species = spvector[i], keys = "all", counts = 0) # species not in GBIF
      cat("species", spvector[i], "has no goereferenced data\n")
      
    }else { # if it has info, use the key with more records, which is the most useful
      if (length(keys) == 1) { # if it is only one key
        key <- keys # detecting species key 
        occ_count[[i]] <- cbind(spvector[i], counts) # count how many records
        
      }else { # if its more than one key
        keysco <- cbind(keys, counts)
        keysco <- keysco[order(keysco[, 2]), ]
        key <- keysco[dim(keysco)[1], 1] # detecting species key that return information
        occ_count[[i]] <- c(Species = spvector[i], keysco[dim(keysco)[1], ])# count how many records
      }
      
      occ <- try(occ_search(taxonKey = key, return = "data", limit = 10000), silent = TRUE) # getting the data from GBIF
      occ_class <- class(occ)
      
      # avoiding errors from GBIF
      while (occ_class[1] == "try-error") {
        occ <- try(occ_search(taxonKey = key, return = "data", limit = 10000), silent = TRUE) # getting the data from GBIF
        occ_class <- class(occ)
        
        if(occ_class[1] != "try-error") {
          break()
        }
      }
      
      # following steps
      occ_g <- occ
      #occ_g <- occ_g[, c(1, 2, 4, 3, 5:dim(occ_g)[2])] # reordering longitude and latitude
      
      # keeping only unique georeferenced records. IF NO FILTERING IS NEEDED, PUT A # IN FRONT OF THE NEXT 3 LINES
      occ_g <- occ_g$data[!is.na(occ_g$data$decimalLatitude) & !is.na(occ_g$data$decimalLongitude), ] # excluding no georeferences
      #occ_g <- occ_g$data[!duplicated(paste(occ_g$data$scientificName, occ_g$data$decimalLatitude, # excluding duplicates
      #                                 occ_g$data$decimalLongitude, sep = "_")), ]
      
      # writing file
      file_name <- paste(gsub(" ", "_", spvector[i]), "csv", sep = ".") # csv file name per each species
      write.csv(occ_g, file_name, row.names = FALSE) # writing inside each genus folder
      
      cat(i, "of", length(spvector), "species\n") # counting species per genus 
    }
  }
}

