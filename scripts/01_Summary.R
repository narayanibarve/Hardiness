
## In office desktop
#setwd("Z:/mybook/narayani/mosquito_ufbi/prism/")

## On Wallace
setwd("/datapool/data/narayani/mosquito_ufbi/prism_daily_download/")

# YearList <- 2016:2020

library(raster)



calculate <- function(VarName, ThresholdList)
{
  
  for (Threshold in ThresholdList)
  {
    
    GetSuitableDays(Threshold, VarName)
      
  }
  
}

## calculate ("tmax", seq(25,27))
## calculate ("tmin", seq(-15, 15))








GetSuitableDays <- function(Threshold, VarName)
{
  YearList <- 2016:2020
  #YearList <- c(2020)
  
  # For Office desktop
  # FileList <- list.files(paste("./", VarName, sep=""), pattern = ".asc", full.names =  TRUE)
  
  FileList <- list.files(paste("./", VarName, sep=""), full.names =  TRUE)
  
  
  ## Get the required file names from the total file name. 
  CurFileList <- c()
  for (i in 1:length(YearList))
  {
    #print(paste("Year ", YearList[i]))
    # subfilelist <- FileList[grep(paste("_", YearList[i], sep = ""), FileList)]
    # print(paste("length of days ", length(subfilelist)))
    # CurFileList <- c(CurFileList, subfilelist)
    
    subdirlist <- FileList[grep(paste("_", YearList[i], sep = ""), FileList)]
    #print(paste("length of days ", length(subfilelist)))
    #CurFileList <- c(CurFileList, subfilelist)
    
    subfilelist = paste(subdirlist, "/", basename(subdirlist), ".bil" , sep = "")
    #print(paste("length of days ", length(subfilelist)))
    CurFileList <- c(CurFileList, subfilelist)
    
    
  }
  
  ## r_l1 = sapply(CurFileList[1:366], raster)
  ## s1  = stack(r_l1)
  
  
  ## Now we will read each file and reclass it for threshold
  ## if the value is below the threshold then it is 1 otherwise 0 in case of min temp
  ## if the value is above the threshold then it is 1 otherwise 0 in case of max temp. 
  ## Suitable = 0 and unsuitable = 1 
  
  s1 <- stack()
  ##for (fileindex in 1:length(CurFileList[1:5]))
  for (fileindex in 1:length(CurFileList))
  {
   
    FileName <- CurFileList[fileindex] 
    ## Read raster
    r1 <- raster(FileName)
    #print(FileName)
    print(fileindex)
    ## For max temperature. 
    if (VarName == "tmax")
    {
      lbltext = "Maximum Temperature"
      ## Assigning 0 for unsuitable and 1 for suitable days for growth. 
      r1cl <- reclassify(r1, c(-100,Threshold,1, Threshold,100,0) )
    }
    
    if (VarName == "tmin")
    {
      ## Assigning 0 for unsuitable and 1 for suitable days for growth. 
      lbltext = "Minimum Temperature"
      r1cl <- reclassify(r1, c(-100,Threshold,0, Threshold,100,1) )
    }
    
    
    par(mfrow=c(1,2))
    #print(paste("range ", range(values(r1), na.rm = TRUE))
    #print(fileindex)
    #plot(r1)
    #plot(r1cl, zlim=c(0,1))
    
    
    s1 <- stack(s1, r1cl)
  }
  
  TotalDays <- sum(s1)
  
  ## OpFileName <- paste("./results/", VarName, "_", min(YearList), "_", max(YearList), ".asc", sep = "") 
  ## writeRaster(TotalDays, OpFileName, overwrite=TRUE)
  
  OpFileName <- paste("/datapool/data/narayani/hardiness/threshold_maps/", 
                      VarName, "_", Threshold , "_", min(YearList), "_", max(YearList), ".asc", sep = "") 
  writeRaster(TotalDays, OpFileName, overwrite=TRUE)
  
  ## Save as jpg file. 
  jpgfilename = paste("/datapool/data/narayani/hardiness/threshold_maps/", VarName, "_", Threshold, "_2016_2020.jpg", sep = "")
  jpeg(filename = jpgfilename, width=8, height=5, units = "in", res = 150)
  plot(TotalDays, main = paste(lbltext, " ", Threshold, "C \n", "2016-2020 (1826 days)", sep = ""), legend.args = list(text = 'Suitable days', side=2, font=2, line=.5, cex=0.8))
  #plot(r1, main = "Maximum Temperature 25C \n 2016-2020 (1826 days)") 
  dev.off()
  
  

}  



## calculate ("tmax", seq(25,27))
## calculate ("tmin", seq(-15, 15))


