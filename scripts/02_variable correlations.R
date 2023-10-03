#####
# Evaluation of variables correlation
#####

# Description
## The following script helps to measure correlation among distinct raster layers
## to be used as predictors in ecological niche modeling. This process needs to
## be performed in the area of model calibration.

## Note: All variables must have the same projection, extent, and resolution.
## We suggest to work with Geographic projections WGS84, with no planar projection.

## The main processes are performed the package ellisenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/ellipsenm.

# loading needed packages (packages will be automatically installed if required)
suppressWarnings({
  if(!require(raster)){
    install.packages("raster")
    library(raster)
  }
}) 
# Installing and loading packages
if(!require(devtools)){
  install.packages("devtools")
}
if(!require(ellipsenm)){
  devtools::install_github("marlonecobos/ellipsenm")
}
library(ellipsenm)

# assuming that you installed ellipsenm, load it, if not installed see 
# https://github.com/marlonecobos/ellipsenm for instructions
library(ellipsenm)

setwd("")
#######################################################################################
# Preparing directory and data
##################

# defining working directory
setwd("C:/Narayani/Projects/Hardiness/threshold_maps/") # change this to your working directory

# IF YOU HAVE THE DATA IN YOUR DIRECTORY AS DESCRIBED BELOW, USE THIS
# variables need to be saved in a subdirectory named "bio", variables must be in 
# ascii format (.asc)

## reading data
varaibles_list <- list.files(path = getwd(), pattern = ".asc", # vector of variables
                             full.names = TRUE)


mnlist = varaibles_list[grep("MN",varaibles_list)]
mxlist = varaibles_list[grep("MX",varaibles_list)]

## making VPD lists. 
vnlist = varaibles_list[grep("VN",varaibles_list)]
vxlist = varaibles_list[grep("VX",varaibles_list)]


# v1 = basename(varaibles_list)
# v2 = strsplit(v1,"_")
# v3 = sapply(v2, "[", 1)
# v4 = cbind(varaibles_list, v3)
# v5 = v4[order(as.numeric(v4[,2])), ]
Sorted_MnList <- sortvars(mnlist)
Sorted_MxList <- sortvars(mxlist)
Sorted_VnList <- sortvars(vnlist)
Sorted_VxList <- sortvars(vxlist)



## AllVar_list <- sortvars(varaibles_list)
AllVar_list <- c(Sorted_MnList, Sorted_MxList)
VPD_list <- sortvars(c(Sorted_VnList,Sorted_VxList))

sortvars <- function(Varfilenames)
{
  v1 = basename(Varfilenames)
  v2 = strsplit(v1,"_")
  v3 = sapply(v2, "[", 1)
  v4 = cbind(Varfilenames, v3)
  v5 = v4[order(as.numeric(v4[,2])), 1]
  return(v5)
  
}

### making the x and y label readable, by adding only the threshold. 
makelabels <- function(varnames)
{
  s1 <- strsplit(varnames, "_")
  s2 <- sapply(s1,"[", 3)
  s3 <- gsub("\\.","-",s2 )
  #names(mn_variables) = s3
  return(s3)
}
#######################################################################################
# Varriable correlation analysis
##################

# functions help
help(variable_correlation)
# to save matrix correlation results see the arguments "save" and "name"

## For minumum temperature
mn_variables <- stack(Sorted_MnList) # stack of variables
# correlation matrix
cors_mn <- variable_correlation(mn_variables)
# checking the table
View(cors_mn) 
## Changing the labels only to threshold.
mn_labels = makelabels(names(mn_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_mn <- variable_correlation(mn_variables, correlation_limit = 0.95,
                              corrplot = TRUE, magnify_to = 1.5, labeltext = mn_labels,
                              main = "\n Minimum Temperature")



## For maximum temperature
mx_variables <- stack(Sorted_MxList) # stack of variables
# correlation matrix
cors_mx <- variable_correlation(mx_variables)
# checking the table
View(cors_mx) 
## Changing the labels only to threshold.
mx_labels = makelabels(names(mx_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_mx <- variable_correlation(mx_variables, correlation_limit = 0.9,
                                 corrplot = TRUE, magnify_to = 1.5, labeltext = mx_labels,
                                 main = "\n Maximum Temperature")



## For minimum and maximum temperature
all_variables <- stack(AllVar_list) # stack of variables
# correlation matrix
cors_all <- variable_correlation(all_variables)
# checking the table
View(cors_all) 
## Changing the labels only to threshold.
all_labels = makelabels(names(all_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_all <- variable_correlation(all_variables, correlation_limit = 0.9,
                                 corrplot = TRUE, magnify_to = 1.5, labeltext = all_labels,
                                 main = "\n Minimum and Maximum Temperature")



## For minimum and maximum temperature
vpd_variables <- stack(VPD_list) # stack of variables
# correlation matrix
cors_vpd <- variable_correlation(vpd_variables)
# checking the table
View(cors_vpd) 
## Changing the labels only to threshold.
vpd_labels = makelabels(names(vpd_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_vpd <- variable_correlation(vpd_variables, correlation_limit = 0.95,
                                  corrplot = TRUE, magnify_to = 1.5, labeltext = vpd_labels,
                                  main = "\n Minimum and Maximum VPD")







############################

## For minumum VPD temperature
vn_variables <- stack(Sorted_VnList) # stack of variables
# correlation matrix
cors_vn <- variable_correlation(vn_variables)
# checking the table
View(cors_vn) 
## Changing the labels only to threshold.
vn_labels = makelabels(names(vn_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_vn <- variable_correlation(vn_variables, correlation_limit = 0.95,
                                 corrplot = TRUE, magnify_to = 1.5, labeltext = vn_labels,
                                 main = "\n Minimum VPD")


## For maximum VPD temperature
vx_variables <- stack(Sorted_VxList) # stack of variables
# correlation matrix
cors_vx <- variable_correlation(vx_variables)
# checking the table
View(cors_vx) 
## Changing the labels only to threshold.
vx_labels = makelabels(names(vx_variables))

# analysis and ploting (all values above 0.75 will be magnified)
cors1_vx <- variable_correlation(vx_variables, correlation_limit = 0.95,
                                 corrplot = TRUE, magnify_to = 1.5, labeltext = vx_labels,
                                 main = "\n Maximum VPD")




# Preparing sets of variables
##################

# simple raster PCA
## functions help
library(kuenm)
help(kuenm_rpca)
setwd("C:/Narayani/Projects/Hardiness/threshold_maps/")

## For all temperature files 
## preparing function's arguments
#var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
#var_folder <- getwd() # name of folder with variables to be combined in distinct sets
#var_folder <- list.files("./", pattern = ".asc$") # name of folder with variables to be combined in distinct sets
var_folder <- AllVar_list
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_temp" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## This is for all threshold files in the threshold maps folder
## i.e all minimum and maximum temperature
## runing PCA
res1 = kuenm_rpca_m(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)

comp1 = raster("./PCA_results_temp/Initial/pc_1.asc")
comp2 = raster("./PCA_results_temp/Initial/pc_2.asc")
comp3 = raster("./PCA_results_temp/Initial/pc_3.asc")

plot(comp1, main = "First PCA component - Temperature")
plot(comp2, main = "Second PCA component - Temperature")
plot(comp3, main = "Third PCA component - Temperature")


################### Maximum Temperature STARTS ###################
## preparing function's arguments
#var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
## var_folder <- "./MaxTemp/" # name of folder with variables to be combined in distinct sets
var_folder <- mxlist # name of folder with variables to be combined in distinct sets 
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_tempmax" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## This is for all threshold files in the threshold maps folder
## i.e all minimum and maximum temperature
## runing PCA
kuenm_rpca_m(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)

comp1_max = raster("./PCA_results_tempmax/Initial/pc_1.asc")
comp2_max = raster("./PCA_results_tempmax/Initial/pc_2.asc")
comp3_max = raster("./PCA_results_tempmax/Initial/pc_3.asc")

plot(comp1_max, main = "First PCA component - Maximum temperature")
plot(comp2_max, main = "Second PCA component - Maximum temperature")
plot(comp3_max, main = "Third PCA component - Maximum temperature")

################### Maximum Temperature Ends ############





################### Minimum Temperature STARTS ###################
## preparing function's arguments
#var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
var_folder <- mnlist # name of folder with variables to be combined in distinct sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_tempmin" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## This is for all threshold files in the threshold maps folder
## i.e all minimum and maximum temperature
## runing PCA
kuenm_rpca_m(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)

comp1_min = raster("./PCA_results_tempmin/Initial/pc_1.asc")
comp2_min = raster("./PCA_results_tempmin/Initial/pc_2.asc")
comp3_min = raster("./PCA_results_tempmin/Initial/pc_3.asc")

plot(comp1_min, main = "First PCA component - Minimum temperature")
plot(comp2_min, main = "Second PCA component - Minimum temperature")
plot(comp3_min, main = "Third PCA component - Minimum temperature")

################### Minimum Temperature Ends ############




################### VPD minimum STARTS ###################
## preparing function's arguments
#var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
## var_folder <- "./MaxTemp/" # name of folder with variables to be combined in distinct sets
var_folder <- vnlist # name of folder with variables to be combined in distinct sets 
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_vpdmin" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## This is for all threshold files in the threshold maps folder
## i.e all minimum and maximum temperature
## runing PCA
kuenm_rpca_m(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)

comp1_max = raster("./PCA_results_vpdmin/Initial/pc_1.asc")
comp2_max = raster("./PCA_results_vpdmin/Initial/pc_2.asc")
comp3_max = raster("./PCA_results_vpdmin/Initial/pc_3.asc")

plot(comp1_max, main = "First PCA component - VPD minimum")
plot(comp2_max, main = "Second PCA component - VPD minimum")
plot(comp3_max, main = "Third PCA component - VPD minimum")

################### VPD minimum Ends ############



################### VPD Maximum STARTS ###################
## preparing function's arguments
#var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
## var_folder <- "./MaxTemp/" # name of folder with variables to be combined in distinct sets
var_folder <- vxlist # name of folder with variables to be combined in distinct sets 
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_vpdmax" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## This is for all threshold files in the threshold maps folder
## i.e all minimum and maximum temperature
## runing PCA
kuenm_rpca_m(variables = var_folder, in.format = in_format, var.scale = scalev, 
             write.result = writer, out.format = out_format, out.dir = out_folder,
             n.pcs = n_pcs)

comp1_max = raster("./PCA_results_vpdmax/Initial/pc_1.asc")
comp2_max = raster("./PCA_results_vpdmax/Initial/pc_2.asc")
comp3_max = raster("./PCA_results_vpdmax/Initial/pc_3.asc")

plot(comp1_max, main = "First PCA component - VPD maximum")
plot(comp2_max, main = "Second PCA component - VPD maximum")
plot(comp3_max, main = "Third PCA component - VPD maximum")

################### VPD Maximum Ends ############




## Doing the bivariate plot between min and max temperature

comp1_min <- raster("./PCA_results_tempmin/Initial/pc_1.asc")
comp1_max <- raster("./PCA_results_tempmax/Initial/pc_1.asc")
comp1_vpmn <- raster("./PCA_results_vpdmin/Initial/pc_1.asc")
comp1_vpmx <- raster("./PCA_results_vpdmax/Initial/pc_1.asc")

comp1stack <- stack(comp1_min, comp1_max, comp1_vpmn, comp1_vpmx) 

comp_pts <- rasterToPoints(comp1stack)
comp_pts <- data.frame(comp_pts)
names(comp_pts) <- c("x", "y", "tempmin", "tempmax", "vpdmn", "vpdmx")
head(comp_pts)

minplotx <- min(comp_pts[,c(3,5)])
maxplotx = max(comp_pts[,c(3,5)])

minploty <- min(comp_pts[,c(4,6)])
maxploty = max(comp_pts[,c(4,6)])


## Plotting minimum and maximum in bivariate plot
plot(comp_pts[,3], comp_pts[,4], col="red", pch = 20, xlim = c(minplotx, maxplotx), 
    ylim = c(minploty, maxploty), xlab = "Minimum (Temperature / VPD)", 
    ylab = "Maximum (Temperature / VPD)")

points(comp_pts[,5], comp_pts[,6], col="blue", pch = 20)

## Keeping X axis as minimum temperate and Y axis with other 3 variables. 
xlim_mn = min(comp_pts[,3])
xlim_mx = max(comp_pts[,3])

ylim_mn = min(comp_pts[,c(4,5,6)])
ylim_mx = max(comp_pts[,c(4,5,6)])

plot(comp_pts[,3], comp_pts[,4], col="red", pch = 20, xlim = c(minplotx, maxplotx), 
     ylim = c(minploty, maxploty), xlab = "Minimum Temperature", 
     ylab = "Maximum Temperature / VPD min / VPD max)")

## VPD min
points(comp_pts[,3], comp_pts[,5], col="blue", pch = 20)
## VPD max
points(comp_pts[,3], comp_pts[,6], col="green", pch = 20)


## With ggplot bivariate
ggplot(comp_pts, aes(tempmin,tempmax )) +
  geom_point(colour = 'red') +
  geom_point(data = comp_pts, aes(vpdmn,vpdmx), colour = 'blue', size = 1)
             


## Temperature
plot(comp_pts[,3], comp_pts[,4], col="red", pch = 20, xlab = "Minimum Temperature", ylab = "Maximum Temperature" )
## VPD
plot(comp_pts[,5], comp_pts[,6], col="blue", pch = 20, xlab = "VPD Minimum", ylab = "VPD Maximum" )



#Bivariate plots of all combinations. 

#dt1 = comp_pts[,3:6]
dt1 = comp_pts
Bivariate <- function(dt1)
{
  
  for (i in 1:ncol(dt1))
  {
    for (j in 1:i)
    {
      print(paste("i ", i, "j ", j))
      if (i !=j)
      {
        #cat(paste("Generating plot for ", names(dt1)[j], " and ", names(dt1)[i], "\n", sep = ""))
        #print(paste(i,"_",j))
        #jpeg(filename=paste("bi_",dt1[j],"_",dt1[i],".jpg", sep = ""),width=1200,height=800) 
        # plot(bpt1[,j+2],bpt1[,i+2], pch= 15, col = "red",xlab=compname[j],ylab=compname[i])
        plot(dt1[,j+2],dt1[,i+2], pch= 15, col = "red",xlab=names(dt1)[j],ylab=names(dt1)[i])
        #points(tm[,j],tm[,i], pch= 19, col = "blue")
        #dev.off()
      }
    
    }
  
  } ## for i 
} ## Bivariate

Bivariate(dt1)





# raster PCA with projections
## preparing function's arguments
var_folder <- "bio_10m_esri" # name of folder with variables to be combined in distinct sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_proj1" # name of folder that will contain the sets 
project <- TRUE
proj_folder <- "Project"
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## runing PCA
kuenm_rpca(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           project = project, proj.vars = proj_folder, n.pcs = n_pcs)





















#' Evaluates correlation among variables
#'
#' @description variable_correlation helps in evaluating correlation among
#' distinct variables.
#'
#' @param variables RasterStack, RasterBrick, or matrix. If matrix, columns
#' represent distinct variables for analysis, otherwise, group of raster layers.
#' @param sample_size (numeric) sample size to be taken from all variables;
#' default = 10000.
#' @param correlation_limit (numeric) absolute value of correlation limit;
#' default = 0.8.
#' @param save (logical) whether or not to save the results; default = FALSE.
#' @param name (character) name of the csv files to be writen;
#' default = "correlation".
#' @param corrplot (logical) whether or not to plot the results; default = FALSE.
#' @param magnify_to (numeric) optional value to be used to magnify all values
#' with absolute correlations above \code{correlation_limit}. Default = NULL.
#' @param ... other arguments to be passed to \code{\link[corrplot]{corrplot}}.
#' Arguments "type", "tl.col", and "tl.srt" are fixed.
#'
#' @return
#' A correlation matrix. If argument \code{corrplot} = TRUE correlation values
#' are shown in a plot.
#'
#' @usage
#' variable_correlation(variables, sample_size = 10000, correlation_limit = 0.8,
#'                      save = FALSE, name = "correlation", corrplot = FALSE,
#'                      magnify_to = NULL, ...)
#'
#' @details
#' If \code{magnify_to} is defined and \code{save} = TRUE, an additional csv
#' file named as "\code{name}_magnified.csv" will be written.
#'
#' @export
#'
#' @examples
#' # raster layers of environmental data
#' vars <- raster::stack(list.files(system.file("extdata", package = "ellipsenm"),
#'                                  pattern = "bio", full.names = TRUE))
#'
#' # simple correlation matrix
#' cors <- variable_correlation(variables, sample_size = 5000)
#'
#' # correlation matrix and plot (values correlated above |0.8| are magnified)
#' cors <- variable_correlation(variables, sample_size = 5000, corrplot = TRUE,
#'                              magnified = 2)
#'
#' # to save results check arguments "save" and "name"

variable_correlation <- function(variables, sample_size = 10000,
                                 correlation_limit = 0.8, save = FALSE,
                                 name = "correlation", corrplot = FALSE,
                                 magnify_to = NULL, labeltext = NULL, ...) {
  # -----------
  # detecting potential errors
  if (missing(variables)) {
    stop("Argument 'variables' is necessary to perform the analysis")
  }
  var_class <- class(variables)[1]
  if (!var_class %in% c("matrix", "RasterStack", "RasterBrick")) {
    stop("'variables' must be either 'matrix', 'RasterStack', or 'RasterBrick'")
  }
  
  # -----------
  # preparing data
  if (var_class != "matrix") {
    ## getting data from the variables
    variables <- na.omit(raster::values(variables))
    
    ## sample of 10000 values if more pixels exist (optional)
    if (nrow(variables) > sample_size) {
      variables <- variables[sample(1:nrow(variables), sample_size), ]
    }
  }
  
  # -----------
  # analyses
  ## correlation matrix calculation
  correlation_matrix <- cor(variables)
  
  ## detecting correlated varaibles more easily
  correlation_matrix1 <- correlation_matrix # making other table with results
  
  max_cor <- correlation_limit # maximum value of correlation allowed
  
  if (!is.null(magnify_to)) {
    mv <- magnify_to
    for (i in 1:dim(correlation_matrix1)[2]) { #correlated values will turn into 2 for easier detection
      for (j in 1:dim(correlation_matrix1)[1]) {
        correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] < -max_cor,
                                            -mv, correlation_matrix1[j, i])
        correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] > max_cor,
                                            mv, correlation_matrix1[j, i])
      }
    }
  }
  
  # -----------
  # write and plot
  # saving correlation matrix
  if (save == TRUE) {
    write.csv(correlation_matrix, paste0(name, ".csv"), row.names = TRUE)
    
    if (!is.null(magnify_to)) {
      name <- paste0(name, "_magnified.csv")
      write.csv(correlation_matrix1, name, row.names = TRUE)
    }
  }
  
  # plotting
  if (corrplot == TRUE) {
    if (!is.null(magnify_to)) {
      cor_mat <- correlation_matrix1/mv
    } else {
      cor_mat <- correlation_matrix
    }
    ## Changing the x and y label text 
    if (!is.null(labeltext)) {
      colnames(cor_mat) = labeltext
      rownames(cor_mat) = labeltext
    }
    
    corrplot::corrplot(cor_mat, type = "upper", tl.col = "black",
                       tl.srt = 45, ...)
  }
  
  # -----------
  # return results
  return(correlation_matrix)
}



#' Principal componens for raster layers and projections
#'
#' @description kuenm_rpca performs a principal component analysis with a set of variables and
#' produces raster layers of them. If needed the pricipal components are projected to other
#' scenarios.
#'
#' @param variables (character or RasterStack) if character, name of the folder where raster layers are located.
#' If RasterStack, stack of raster layers to be used in principal component analyses.
#' @param in.format (character) valid only if \code{variables} is character. Format of variables in the directory.
#' Options are "ascii", "GTiff", and "EHdr" = bil.
#' @param var.scale (logical) wheter or not to scale variables before performing principal component
#' analyses. Default = TRUE.
#' @param write.result (logical) whether or not to write PCA results and raster layers (PCs) in \code{out.dir}.
#' @param out.dir (character) valid if \code{write.result} = TRUE. Name of the folder to be created to save the
#' results of the analyses. Default = "PCA_results".
#' @param out.format (character) if \code{write.result} = TRUE, format of variables to be written in distinct
#' sets inside \code{out.dir}. Options are "ascii", "GTiff", and "EHdr" = bil. Default = "GTiff".
#' @param project (logical) whether or not to project the species niche to other scenario(s).
#' If TRUE, argument \code{proj.variables} needs to be defined. Default = FALSE.
#' @param proj.vars (character or RasterStack) if character, name of the folder where subfolders with environmental
#' variables of scenarios for projections are (useful if multiple projections are needed). If RasterStack, object
#' containing stacked variables of only one projection scenario. Variables must correspond with variables in \code{vars.folder}
#' (i.e., their names must correspond but they should represent conditions in other scenario).
#' @param n.pcs (numeric) number of principal components to be returned as rasters. By default all principal
#' components are returned as RasterLayers.
#'
#' @return
#' A list containing PCA loadings and PCA summary as matrices, as well as one or multiple (if projected) RasterStacks
#' of principal components.
#'
#' If \code{write.result} = TRUE, all results are written in \code{out.dir}.
#'
#' @details
#' If \code{var.scale} = TRUE, variables are centered to cero and scaled using \code{\link[base]{scale}}.
#'
#' @usage
#' kuenm_rpca(variables, in.format, var.scale = TRUE, write.result = TRUE,
#'            out.format = "GTiff", out.dir = "PCA_results", project = FALSE,
#'            proj.vars, n.pcs)
#'
#' @export
#'
#' @examples
#' # Data
#' variab <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "Mbio_", full.names = TRUE))
#' names(variab) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' proj_var <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                      pattern = "Gbio_", full.names = TRUE))
#' names(proj_var) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' # Example with no projection
#' npcs <- 3
#'
#' rpca <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, n.pcs = npcs)
#'
#' # Example with projection
#' project <- TRUE
#'
#' rpca1 <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, project = project,
#'                     proj.vars = proj_var, n.pcs = npcs)


## runing PCA
#kuenm_rpca(variables = var_folder, in.format = in_format, var.scale = scalev, 
#           write.result = writer, out.format = out_format, out.dir = out_folder,
#           project = project, proj.vars = proj_folder, n.pcs = n_pcs)



kuenm_rpca_m <- function(variables, in.format, var.scale = TRUE, write.result = TRUE, out.format = "GTiff",
                       out.dir = "PCA_results", project = FALSE, proj.vars, n.pcs) {
  
  # testing potential errors
  if (missing(variables)) {
    stop("Argument variables must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(proj.vars)) {
      stop("If projections are needed, argument proj.vars must be defined. See functions help.")
    }
  }
  
  # formatting
  if (class(variables)[1] == "character") {
    if (missing(in.format)) {
      stop("Argument variables is a character, in.format needs to be defined.")
    }
    if (in.format == "ascii") {
      patt <- ".asc$"
    }
    if (in.format == "GTiff") {
      patt <- ".tif$"
    }
    if (in.format == "EHdr") {
      patt <- ".bil$"
    }
  }
  
  if (!missing(write.result)) {
    if (out.format == "ascii") {
      patt1 <- ".asc"
    }
    if (out.format == "GTiff") {
      patt1 <- ".tif"
    }
    if (out.format == "EHdr") {
      patt1 <- ".bil"
    }
  }
  
  # reading variables
  if (class(variables)[1] == "character") {
    # var <- list.files(variables, pattern = patt, full.names = TRUE)
    # variables <- raster::stack(var)
    variables <- raster::stack(variables)
  }
  
  var_points <- na.omit(raster::values(variables))
  
  # pca analyses
  if (var.scale == TRUE) {
    pca <- prcomp(var_points, center = TRUE, scale = TRUE)
  } else {
    pca <- prcomp(var_points, center = TRUE, scale = FALSE)
  }
  
  scores <- pca$x
  
  if (missing(n.pcs)) {
    n.pcs <- length(var)
  }
  
  pcras <- list()
  
  if (write.result == TRUE) {
    cat("\nWriting raster PCs in Output folder, please wait...\n")
    dir.create(out.dir)
    
    pca_fol <- paste(out.dir, "Initial", sep = "/")
    dir.create(pca_fol)
  }
  
  for (i in 1:n.pcs) {
    pcra <- variables[[1]]
    pcra[!is.na(raster::values(pcra))] <- scores[, i]
    
    if (write.result == TRUE) {
      filenam <- paste(pca_fol, "/pc_", i, patt1, sep = "")
      raster::writeRaster(pcra, filenam, format = out.format)
    }
    
    pcras[[i]] <- pcra
  }
  
  pcras <- do.call(raster::stack, pcras)
  names(pcras) <- paste0("pc_", 1:dim(pcras)[3])
  
  StdDev <- pca$sdev
  VarExp <- pca$sdev^2/sum(pca$sdev^2)
  CumVar <- cumsum(VarExp)
  SumPCAMat <- rbind(StdDev, VarExp, CumVar)
  colnames(SumPCAMat) <- paste("PC", seq(1, length(StdDev)), sep = "")
  row.names(SumPCAMat) <- c("Standard deviation", "Proportion of Variance",
                            "Cumulative Proportion")
  
  if (write.result == TRUE) {
    sink(paste(paste(pca_fol, "pca_results.txt", sep = "/")))
    cat("Principal component analysis results\n")
    cat("\nPCA loadings\n")
    print(pca$rotation, nrow = nrow(pca$rotation))
    #cat(pca$rotation)
    
    cat("\n\nPCA summary\n")
    print(SumPCAMat)
    #cat(SumPCAMat)
    sink()
  }
  
  # pca results to be returned
  loadings <- pca$rotation
  respca <- SumPCAMat
  
  # projecting PCs
  if (project == TRUE) {
    ppcrass <- list()
    
    if (write.result == TRUE) {
      cat("\nProjecting and writing projected raster PCs in Output folder, please wait...\n")
    } else {
      cat("\nProjecting raster PCs\n")
    }
    
    if (class(proj.vars)[1] == "character") {
      proj_dirs <- list.dirs(proj.vars, recursive = FALSE)
      proj_names <- list.dirs(proj.vars, recursive = FALSE, full.names = FALSE)
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }
    if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }
    
    
    for (h in 1:length(proj_dirs)) {
      if (class(proj.vars)[1] == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
        p_stack <- raster::stack(pvar)
      }
      if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
        p_stack <- proj.vars
      }
      if (write.result == TRUE) {
        dir.create(fol_names[h])
      }
      
      ppcras <- list()
      
      p_stackp <- na.omit(raster::values(p_stack))
      colnames(p_stackp) <- names(pca[[4]])
      p_pcs <- predict(pca, newdata = p_stackp)
      
      for (i in 1:n.pcs) {
        pcra <- p_stack[[1]]
        pcra[!is.na(raster::values(pcra))] <- p_pcs[, i]
        
        if (write.result == TRUE) {
          filenam <- paste(fol_names[h], "/pc_", i, patt1, sep = "")
          raster::writeRaster(pcra, filenam, format = out.format)
        }
        
        ppcras[[i]] <- pcra
      }
      
      ppcrass[[h]] <- do.call(raster::stack, ppcras)
      names(ppcrass[[h]]) <- paste0("pc_", 1:dim(ppcrass[[h]])[3])
    }
    
    names(ppcrass) <- paste("PCRasters", proj_names, sep = "_")
  }
  
  if (project == TRUE) {
    results <- c(list(loadings, respca, pcras), ppcrass)
    names(results)[1:3] <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }else {
    results <- list(loadings, respca, pcras)
    names(results) <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }
  
  if (write.result == TRUE) {
    cat("\nRaster PCA finished. Check your output directory", paste(getwd(), out.dir, sep = "/"), "\n")
  }
  
  return(results)
}









#############################################################
################ This is from the pacakge ###################




#' Principal componens for raster layers and projections
#'
#' @description kuenm_rpca performs a principal component analysis with a set of variables and
#' produces raster layers of them. If needed the pricipal components are projected to other
#' scenarios.
#'
#' @param variables (character or RasterStack) if character, name of the folder where raster layers are located.
#' If RasterStack, stack of raster layers to be used in principal component analyses.
#' @param in.format (character) valid only if \code{variables} is character. Format of variables in the directory.
#' Options are "ascii", "GTiff", and "EHdr" = bil.
#' @param var.scale (logical) wheter or not to scale variables before performing principal component
#' analyses. Default = TRUE.
#' @param write.result (logical) whether or not to write PCA results and raster layers (PCs) in \code{out.dir}.
#' @param out.dir (character) valid if \code{write.result} = TRUE. Name of the folder to be created to save the
#' results of the analyses. Default = "PCA_results".
#' @param out.format (character) if \code{write.result} = TRUE, format of variables to be written in distinct
#' sets inside \code{out.dir}. Options are "ascii", "GTiff", and "EHdr" = bil. Default = "GTiff".
#' @param project (logical) whether or not to project the species niche to other scenario(s).
#' If TRUE, argument \code{proj.variables} needs to be defined. Default = FALSE.
#' @param proj.vars (character or RasterStack) if character, name of the folder where subfolders with environmental
#' variables of scenarios for projections are (useful if multiple projections are needed). If RasterStack, object
#' containing stacked variables of only one projection scenario. Variables must correspond with variables in \code{vars.folder}
#' (i.e., their names must correspond but they should represent conditions in other scenario).
#' @param n.pcs (numeric) number of principal components to be returned as rasters. By default all principal
#' components are returned as RasterLayers.
#'
#' @return
#' A list containing PCA loadings and PCA summary as matrices, as well as one or multiple (if projected) RasterStacks
#' of principal components.
#'
#' If \code{write.result} = TRUE, all results are written in \code{out.dir}.
#'
#' @details
#' If \code{var.scale} = TRUE, variables are centered to cero and scaled using \code{\link[base]{scale}}.
#'
#' @usage
#' kuenm_rpca(variables, in.format, var.scale = TRUE, write.result = TRUE,
#'            out.format = "GTiff", out.dir = "PCA_results", project = FALSE,
#'            proj.vars, n.pcs)
#'
#' @export
#'
#' @examples
#' # Data
#' variab <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "Mbio_", full.names = TRUE))
#' names(variab) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' proj_var <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                      pattern = "Gbio_", full.names = TRUE))
#' names(proj_var) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' # Example with no projection
#' npcs <- 3
#'
#' rpca <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, n.pcs = npcs)
#'
#' # Example with projection
#' project <- TRUE
#'
#' rpca1 <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, project = project,
#'                     proj.vars = proj_var, n.pcs = npcs)


kuenm_rpca <- function(variables, in.format, var.scale = TRUE, write.result = TRUE, out.format = "GTiff",
                       out.dir = "PCA_results", project = FALSE, proj.vars, n.pcs) {
  
  # testing potential errors
  if (missing(variables)) {
    stop("Argument variables must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(proj.vars)) {
      stop("If projections are needed, argument proj.vars must be defined. See functions help.")
    }
  }
  
  # formatting
  if (class(variables)[1] == "character") {
    if (missing(in.format)) {
      stop("Argument variables is a character, in.format needs to be defined.")
    }
    if (in.format == "ascii") {
      patt <- ".asc$"
    }
    if (in.format == "GTiff") {
      patt <- ".tif$"
    }
    if (in.format == "EHdr") {
      patt <- ".bil$"
    }
  }
  
  if (!missing(write.result)) {
    if (out.format == "ascii") {
      patt1 <- ".asc"
    }
    if (out.format == "GTiff") {
      patt1 <- ".tif"
    }
    if (out.format == "EHdr") {
      patt1 <- ".bil"
    }
  }
  
  # reading variables
  if (class(variables)[1] == "character") {
    var <- list.files(variables, pattern = patt, full.names = TRUE)
    variables <- raster::stack(var)
  }
  
  var_points <- na.omit(raster::values(variables))
  
  # pca analyses
  if (var.scale == TRUE) {
    pca <- prcomp(var_points, center = TRUE, scale = TRUE)
  } else {
    pca <- prcomp(var_points, center = TRUE, scale = FALSE)
  }
  
  scores <- pca$x
  
  if (missing(n.pcs)) {
    n.pcs <- length(var)
  }
  
  pcras <- list()
  
  if (write.result == TRUE) {
    cat("\nWriting raster PCs in Output folder, please wait...\n")
    dir.create(out.dir)
    
    pca_fol <- paste(out.dir, "Initial", sep = "/")
    dir.create(pca_fol)
  }
  
  for (i in 1:n.pcs) {
    pcra <- variables[[1]]
    pcra[!is.na(raster::values(pcra))] <- scores[, i]
    
    if (write.result == TRUE) {
      filenam <- paste(pca_fol, "/pc_", i, patt1, sep = "")
      raster::writeRaster(pcra, filenam, format = out.format)
    }
    
    pcras[[i]] <- pcra
  }
  
  pcras <- do.call(raster::stack, pcras)
  names(pcras) <- paste0("pc_", 1:dim(pcras)[3])
  
  StdDev <- pca$sdev
  VarExp <- pca$sdev^2/sum(pca$sdev^2)
  CumVar <- cumsum(VarExp)
  SumPCAMat <- rbind(StdDev, VarExp, CumVar)
  colnames(SumPCAMat) <- paste("PC", seq(1, length(StdDev)), sep = "")
  row.names(SumPCAMat) <- c("Standard deviation", "Proportion of Variance",
                            "Cumulative Proportion")
  
  if (write.result == TRUE) {
    sink(paste(paste(pca_fol, "pca_results.txt", sep = "/")))
    cat("Principal component analysis results\n")
    cat("\nPCA loadings\n")
    print(pca$rotation)
    
    cat("\n\nPCA summary\n")
    print(SumPCAMat)
    sink()
  }
  
  # pca results to be returned
  loadings <- pca$rotation
  respca <- SumPCAMat
  
  # projecting PCs
  if (project == TRUE) {
    ppcrass <- list()
    
    if (write.result == TRUE) {
      cat("\nProjecting and writing projected raster PCs in Output folder, please wait...\n")
    } else {
      cat("\nProjecting raster PCs\n")
    }
    
    if (class(proj.vars)[1] == "character") {
      proj_dirs <- list.dirs(proj.vars, recursive = FALSE)
      proj_names <- list.dirs(proj.vars, recursive = FALSE, full.names = FALSE)
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }
    if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }
    
    
    for (h in 1:length(proj_dirs)) {
      if (class(proj.vars)[1] == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
        p_stack <- raster::stack(pvar)
      }
      if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
        p_stack <- proj.vars
      }
      if (write.result == TRUE) {
        dir.create(fol_names[h])
      }
      
      ppcras <- list()
      
      p_stackp <- na.omit(raster::values(p_stack))
      colnames(p_stackp) <- names(pca[[4]])
      p_pcs <- predict(pca, newdata = p_stackp)
      
      for (i in 1:n.pcs) {
        pcra <- p_stack[[1]]
        pcra[!is.na(raster::values(pcra))] <- p_pcs[, i]
        
        if (write.result == TRUE) {
          filenam <- paste(fol_names[h], "/pc_", i, patt1, sep = "")
          raster::writeRaster(pcra, filenam, format = out.format)
        }
        
        ppcras[[i]] <- pcra
      }
      
      ppcrass[[h]] <- do.call(raster::stack, ppcras)
      names(ppcrass[[h]]) <- paste0("pc_", 1:dim(ppcrass[[h]])[3])
    }
    
    names(ppcrass) <- paste("PCRasters", proj_names, sep = "_")
  }
  
  if (project == TRUE) {
    results <- c(list(loadings, respca, pcras), ppcrass)
    names(results)[1:3] <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }else {
    results <- list(loadings, respca, pcras)
    names(results) <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }
  
  if (write.result == TRUE) {
    cat("\nRaster PCA finished. Check your output directory", paste(getwd(), out.dir, sep = "/"), "\n")
  }
  
  return(results)
}


