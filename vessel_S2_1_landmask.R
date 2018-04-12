#' @title Land mask for vessel detection on Sentinel-2 images (batch mode)
#' @description First step of the automatic vessel detection/classification (batch) procedure, where land areas on Sentinel-2 (Level 1C) images are removed 
#' @author Ursa Kanjir
#' @date "'r format(Sys.Date())'"
#' @return a raster image with land areas (larger cloud areas also) removed, only sea areas left for further processing 

library(raster)
library(rgdal)
library(gdalUtils)
library(pastecs)

#-----Function for finding local minima in ndwi----------
locmin <- function(band) 
{
  dens_band <- density(band)
  ts_y <- ts(dens_band$y)  #makes it a time series
  tp <- turnpoints(ts_y)  # analysing turning points (peaks or pits)
  
  points(dens_band$x[tp$tppos], dens_band$y[tp$tppos], col="red")
  tab_pnts <- cbind(dens_band$x[tp$tppos], dens_band$y[tp$tppos])
  tab_pnts_x <- tab_pnts[,1]
  thresh_x <- tab_pnts_x < 0.7 & tab_pnts_x > 0
  
  threshold <- min(tab_pnts[thresh_x,1]) 
  return(threshold)
}
#----------------------

mainDir <- "define your directory where your Sentinel-2 data is stored"
setwd(mainDir)
outDir <- "define your directory where you want to save your output (land masks)"
ptm <- proc.time()

s2.folders <- list.files(mainDir)

for (i in 1:length(s2.folders)){
  #load bands that will be processed in the next steps - bands 2,3,4,8 and 11
  fileList <- list.files(s2.folders[i], pattern = "*_B02.jp2$|*_B03.jp2$|*_B04.jp2$|_B08.jp2$|_B11.jp2$", recursive = TRUE, full.names = TRUE)

  #-------------If Sentinel-2 are not yet in TIF format conver them using this procedure--------------------------
  #-------------If they already are in TIF format you can skip next 12 lines-------------------------------------- 
   gdal_setInstallation()  #if gdal is not found this will be executed to find a working GDAL that has the right drivers
   valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
   if(require(raster) && require(rgdal) && valid_install)
   {
     for(j in 1:length(fileList)){
       src_dataset <- fileList[j]
       otpt_dataset <- paste0(substr(fileList[j],1, nchar(fileList[j])-3), "tif")
       # Original gdal_translate call:
       gdal_translate(src_dataset, otpt_dataset, of="GTiff", output_Raster = T) 
     }
   }
  #--------------------------------------------------------------------------------------------------------
  
  #different naming according to older or newer Sentinel-2 images
  if (grepl("OPER", fileList[i])) {  #naming for older S-2 images than 6.12.2016
    names1 <- substr(s2.folders[i],1, nchar(s2.folders[i])-50)
  } else {    #naming for newer S-2 images than 6.12.2016
    names1 <- substr(s2.folders[i],1, nchar(s2.folders[i])-46)
  }
  
  tifList <- list.files(s2.folders[i], pattern = "*.tif$", recursive = TRUE, full.names = TRUE)
  
  #stack of all five bands is not functioning as SWIR band has spatial resolution of 20 m (others are 10 m)
  #downsample of Green band on 20 m, same resolution as SWIR
  
  downs.green <- aggregate(raster(tifList[2]), fact=2, FUN=mean)
  #plot(downs.green)
  #writeRaster(downs.green, filename = paste0(ourDir, "bgreen_downsampled_mean.tif"), format="GTiff", overwrite=T)
  swir <- raster(tifList[5])
  
  #Calculate MNDWI from the (downsampled) Green and SWIR band
  mndwi <- (downs.green - swir)/(downs.green + swir) 
  #plot(mndwi)
  writeRaster(mndwi, filename = paste0(outDir, names1, "_mndwi"), format = "GTiff", overwrite=TRUE)

  hist_mndwi <- hist(mndwi)
  
  #calculate threshold by finding local minima
  minima <- locmin(mndwi)
 
  #apply threshold on the mndwi mask
  maskmndwi <- mndwi > minima  
  #plot(maskmndwi)

  #apply smoothing on the mask - to avoid too much false alarms in later stages (IMPORTANT for later phases!!!)
  #f = matrix(1, nrow=5, ncol=5) 
  
  smooth_maskmndwi_1 <- focal(maskmndwi, w=matrix(1,5,5), fun=max, pad = TRUE) 
  smooth_maskmndwi <- focal(smooth_maskmndwi_1, w = matrix(1,7,7), fun = min, pad = TRUE)
  
  #plot(smooth_maskmndwi)
  
  writeRaster(smooth_maskmndwi, filename= paste0(outDir, names1, "_mask_mndwi.tif"), format="GTiff", overwrite=TRUE)
  
  #Binary mask upsample back on 10 m spatial resolution
  resampleFactor <- .5  # reduce the cell size by 50% and double the number of rows and columns.      
  inCols <- ncol(smooth_maskmndwi)
  inRows <- nrow(smooth_maskmndwi)
  resampledRaster <- raster(ncol=(inCols / resampleFactor), nrow=(inRows / resampleFactor))
  extent(resampledRaster) <- extent(smooth_maskmndwi)
  maskmndwi10m <- resample(smooth_maskmndwi, resampledRaster, method='bilinear',filename=paste0(outDir, "testOutResamp2.tif"),overwrite=TRUE) #in bilinearna in nearest neighbor isto naredita!!! no value tam kjer je land, zakaj???
  crs(maskmndwi10m) <- "+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  #plot(maskmndwi10m)
  #writeRaster(maskmndwi10m, filename = paste0(outDir, names1, "/maskmaskmndwi10m.tif"), format="GTiff", overwrite=T)
  
  #stack layers of the same 10 m resolution into multispectral band
  img1 <- stack(raster(tifList[1]), raster(tifList[2]), raster(tifList[3]), raster(tifList[4]))
  writeRaster(img1, filename = paste0(outDir, names1, "_stack.tif"), format = "GTiff", overwrite=TRUE)
  img2 <- img1
  
  #Apply mask on the all bands
  for (k in 1:nlayers(img2)) 
  {
    img2[[k]][maskmndwi10m[[1]][,]==0] <- -9999
  }
  
  plotRGB(img2, r=3, g=2, b=1, stretch="hist")
  writeRaster(img2, filename = paste0(outDir, names1, "_seamask.tif"), format = "GTiff", overwrite=TRUE)

}

proc.time() - ptm
