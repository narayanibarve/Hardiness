
library(raster)

# defining working directory
setwd("C:/Narayani/Projects/Hardiness/threshold_maps/") # 


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

## This is with 3 classes in the component file. (c(0,0.33,0.66,1)) 
  comp_pts$tempmin_q = DoQuantiles(comp_pts$tempmin)
  comp_pts$tempmax_q = DoQuantiles(comp_pts$tempmax)
  comp_pts$vpdmn_q = DoQuantiles(comp_pts$vpdmn)
  comp_pts$vpdmx_q = DoQuantiles(comp_pts$vpdmx)
  comp_pts$FinalClass = as.numeric(paste(comp_pts$tempmin_q, comp_pts$tempmax_q, comp_pts$vpdmn_q, comp_pts$vpdmx_q, sep = ""))
  
  r1 = raster("./PCA_results_tempmax/Initial/pc_1.asc")
  r2 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalClass)
  writeRaster(r2, "./threshold_maps/map/Classified_comp.asc")

  write.csv(comp_pts, "./comp_classified.csv", row.names = FALSE)

## Save individual env variables as classified maps
  
  #comp_pts$FinalTmp = as.numeric(paste(comp_pts$tempmin_q, comp_pts$tempmax_q, sep = ""))
  i1 = rasterize(comp_pts[, c(1,2)], r1, field = as.numeric(comp_pts$tempmin_q))
  writeRaster(i1, "./map/TN_cls.asc")
  
  i2 = rasterize(comp_pts[, c(1,2)], r1, field = as.numeric(comp_pts$tempmax_q))
  writeRaster(i2, "./map/TX_cls.asc")
  
  i4 = rasterize(comp_pts[, c(1,2)], r1, field = as.numeric(comp_pts$vpdmn_q))
  writeRaster(i4, "./map/VN_cls.asc")
  
  i5 = rasterize(comp_pts[, c(1,2)], r1, field = as.numeric(comp_pts$vpdmx_q))
  writeRaster(i5, "./map/VX_cls.asc")
  
    
## Only temperature related variables. 
comp_pts$FinalTmp = as.numeric(paste(comp_pts$tempmin_q, comp_pts$tempmax_q, sep = ""))
r3 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalTmp)
writeRaster(r3, "./map/Classified_comp_Tmp.asc")


## Only temperature related variables.max,min
comp_pts$FinalTmxTmn = as.numeric(paste(comp_pts$tempmax_q, comp_pts$tempmin_q, sep = ""))
r3 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalTmxTmn)
writeRaster(r3, "./map/Classified_comp_TmxTmn.asc")




## Only VPD related variables. 
comp_pts$FinalVpd = as.numeric(paste(comp_pts$vpdmx_q, comp_pts$vpdmn_q, sep = ""))
r4 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalVpd)
writeRaster(r4, "./map/Classified_comp_Vpd.asc")


## Only TempMin, VPDMax  variables. 
comp_pts$FinalTmnVmx = as.numeric(paste(comp_pts$tempmin_q, comp_pts$vpdmx_q, sep = ""))
r5 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalTmnVmx)
writeRaster(r5, "./map/Classified_comp_TmnVmx.asc", overwrite=TRUE)



## Tmn, Tmx, Vmx
comp_pts$FinalTmnTmxVmx = as.numeric(paste(comp_pts$tempmin_q, comp_pts$tempmax_q, comp_pts$vpdmx_q, sep = ""))
r6 = rasterize(comp_pts[, c(1,2)], r1, field = comp_pts$FinalTmnTmxVmx)
writeRaster(r6, "./map/Classified_comp_TmnTmxVmx.asc")


##########################################
## This is for quartiles. (4 classes)

comp_4cls <- rasterToPoints(comp1stack)
comp_4cls <- data.frame(comp_4cls)
names(comp_4cls) <- c("x", "y", "tempmin", "tempmax", "vpdmn", "vpdmx")
head(comp_4cls)



## This is with 3 classes in the component file. (c(0,0.33,0.66,1)) 
comp_4cls$tempmin_q = DoQuantiles(comp_4cls$tempmin)
comp_4cls$tempmax_q = DoQuantiles(comp_4cls$tempmax)
comp_4cls$vpdmn_q = DoQuantiles(comp_4cls$vpdmn)
comp_4cls$vpdmx_q = DoQuantiles(comp_4cls$vpdmx)
comp_4cls$FourFour = as.numeric(paste(comp_4cls$tempmin_q, comp_4cls$tempmax_q, comp_4cls$vpdmn_q, comp_4cls$vpdmx_q, sep = ""))

r1 = raster("./PCA_results_tempmax/Initial/pc_1.asc")
r2 = rasterize(comp_4cls[, c(1,2)], r1, field = comp_4cls$FourFour)
writeRaster(r2, "./map/FourFour.asc", overwrite=TRUE)

write.csv(comp_pts, "./comp_classified.csv", row.names = FALSE)








DoQuantiles <- function(data1)
{
  #maxval = max(data1, na.rm = TRUE)
  qvals = quantile(data1, na.rm = T, probs = c(0,0.33,0.66,1))
  class1 = which(data1 >= qvals[[1]] & data1 < qvals[[2]])
  class2 = which(data1 >= qvals[[2]] & data1 < qvals[[3]])
  class3 = which(data1 >= qvals[[3]] & data1 <= qvals[[4]])
  #class4 = which(data1 >= qvals[[4]] & data1 <= qvals[[5]])
  data1[class1] = "1"
  data1[class2] = "2"  
  data1[class3] = "3"  
  #data1[class4] = "4"  
  return(data1)
}

