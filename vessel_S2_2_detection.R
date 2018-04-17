#' @title Vessel detection module from Sentinel-2 landmasked images (batch mode)
#' @description Algorithm for automatic vessel detection which takes land mask image as an input and returns shapefile polygons of 
#potential vessels with a set of attributes (in a batch mode). 
#' @author Ursa Kanjir
#' @date "'r format(Sys.Date())'"
#' @return a shapefile polygons with a set of calculated geometrical and spectral attributes 

require(rgeos)
require(raster)
require(rgdal)
require(sp)

#for this code to function you need to call gdal_polygonizeR.R function (added in the VesselDet_S2 list)

mainDir <- "The directory where your landmasked images are stored"  
outDir <- "The directory where you want to store your outputs" 

#img input data is land mask created in the previous step
lm_List <- list.files(mainDir, pattern = "*_seamask.tif$", recursive = TRUE, full.names = TRUE)

ptm <- proc.time()

for (i in 1:length(lm_List)){

  name_full <- sub(".*/", "", lm_List[i])
  name <- substr(name_full, 1, nchar(name_full) - 4)
  name
  #load an image
  img <- stack(lm_List[i])
  #calculate vessel index on the image
  #vessel index is a subtraction of NIR and R band
  vess_indx <- (img[[4]] - img[[3]])
  #hist_vess_indx <- hist(vess_indx)
  #plot(hist_vess_indx)
  #writeRaster(vess_indx, filename = paste0(outDir, name, "vess_indx.tif"), format="GTiff", overwrite=TRUE)
  
  #to avoid loosing rubber and very bright boats from the image, which in the case of S-2 are not detected using only vessel index
  #we apply an additional condition where values of all 4 bands are higher than certain threshold value of mean/standard deviation
  rasterRescale<-function(r){
    ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
  }
  
  r2<- rasterRescale(img) #this step can take a lot of memory when processing a large image

  #sum of all the bands in a certain pixel
  img_tog <- sum(r2)
  #plot(img_tog)
  #writeRaster(img_tog, filename = paste0(outDir, name, "_testS2_okt.tif"), format="GTiff", overwrite=TRUE)
  
  #cellStat - statistics of the cells of a raster object
  img.mean <- cellStats(img_tog, 'mean', na.rm=TRUE) #calculation of mean of all the pixels (memory safe)
  img.sd <- cellStats(img_tog, 'sd', na.rm=TRUE) #calculation of standard deviation of all the pixels (memory safe)
  
  img.meja <- img.mean + 0.55*img.sd #0.55 is pragmatically selected threshold for S-2 images

  #the condition is different for every img separately
  vessels <- (vess_indx > 0) | (img_tog > img.meja)
  #plot(vessels)
  #writeRaster(vessels, filename = paste0(outDir, name, "_vessels.tif"), format="GTiff", overwrite=TRUE)
  
  rm(r2, vess_indx) #these two variables are large, therefore they shoudl be removed
  
  #Enhancing vessel mask to avoid errors in the further steps
  vessels_smooth <- focal(vessels, w = matrix(1,3,3), fun = max, pad = TRUE)
  vessels_smooth_fin <- focal(vessels_smooth, w = matrix(1,3,3), fun = min, pad = TRUE)
  #plot(vessels_smooth_fin)
  
  writeRaster(vessels_smooth_fin, filename = paste0(outDir, name, "_poss_vessel.tif"), format="GTiff", overwrite=TRUE)

  ################Obtain vector polygons from raster of possible vectors
   
  #Fast polygonisation calculation - BUT needs gdal_polygonize.py
  #Function  downloaded from https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
  system.time(vect_vessels <- gdal_polygonizeR(vessels_smooth_fin, outshape = paste0(outDir, name, "_vector_temp")))
  #vector and raster layer has same spatial extent (projection)
  
  #using a zero width buffer cleans up many topology problems
  vect_vessels <- gBuffer(vect_vessels, byid=TRUE, width=0)
  
  #calculate area of shape polygons 
  system.time(vess_area <- gArea(vect_vessels, byid = TRUE))
  #assign area values to the polygon you are interested in
  vect_vessels@data$area <- vess_area
  
  #define minimal and maximal area value for different sensors
  #The condition for ideal max and min vessel size are mentioned in the article
  res <- 10 
  min.val <- (res^2 * 3) #m2 (after Bannister and Neyland, 2015)
  max.val <- 30000  #to potrebno še stestirat!
  
  print(paste0("minimal vessel size is ", min.val, "m2"))
  print(paste0("maximal vessel size is ", max.val, "m2")) 
  #very big ships are detected per partes, cause their surface is very heterogene
  
  #remove all those polygons with DN value = 0 (we focus only on those with the value = 1 - possible ship)
  #remove all the polygons that are smaller than min and bigger than max values
  vect_vessels <- vect_vessels[(vect_vessels@data$area > min.val &
                                  vect_vessels@data$area < max.val & 
                                  vect_vessels@data$DN == 1),]
  #writeOGR(vect_vessels, paste0(outDir, name, "_vector_vmesni"), paste0(name, "_vector_vmesni"), driver="ESRI Shapefile", overwrite_layer = TRUE)
  
  print(paste0("Image ", i, " includes ", length(vect_vessels), " possible vessel objects."))
  
  #------------------------------
  #CALCULATE geometrical attributes of all obtained polygons
  
  #Calculate position, heading, width, length
  #set angle sequence 
  set_angle <- 10 # set angle sequence - for how much it should move between 0-180
  angle <- seq (0, 179, set_angle)
  #set empty matrix where calculated angles will be later stored
  fin_res <- list()
  mat_res <- matrix(nrow=length(angle), ncol=2) # matrix with the results of distances of every angle
  colnames(mat_res) <- c("heading", "length")
  mat_res2 <- matrix(nrow=1, ncol=5)
  colnames(mat_res2) <- c("x", "y", "heading", "length", "width")
  
  #add ID column to the vector table
  vect_vessels@data$idloop <- seq(1,length(vect_vessels),1)
  #head(vect_vessels@data)
  
  for (j in 1:length(vect_vessels)){
    
    vect_vess.p <- vect_vessels[(vect_vessels@data$idloop%in%j),]  #so it selects only one polygon out of many
    plot(vect_vess.p)
    
    #Calculate vessel polygon width and length
    #add centroid to the layer
    vec_centr <- gCentroid(vect_vess.p, byid = TRUE) #vec_ves_s
    points(vec_centr)  #plot points

    for (k in 1:length(angle)){
      angle2 <- angle[k]
      
      if (angle2 <= 45){
        xval = 500 * tan(angle2*pi/180)
        yval = 500
      } else if (45 < angle2 & angle2 < 135){
        xval = 500
        yval = 500/tan(angle2*pi/180)
      }else if (135 <= angle2){
        xval = tan(angle2*pi/180)*(-500)
        yval = -500
      } 

      #make vertical line
      x_coord <- c(vec_centr$x+xval, vec_centr$x-xval)
      y_coord <- c(vec_centr$y+yval, vec_centr$y-yval)

      l <- cbind(x_coord, y_coord)
      sl <- Line(l)

      s <- Lines(list(sl), ID = "a")
      myLine <- SpatialLines(list(s), proj4string = vec_centr@proj4string) #with proj4string you define whic projection new line should have
      #plot(myLine, add = TRUE)

      #cut vertical line
      transect <- gIntersection(vect_vess.p, myLine, byid = TRUE)
      plot(transect, col = "blue", add=TRUE)
    
      #calculate the length of the intersected line
      length_tr1 <- gLength(transect)
      
      mat_res[k,] <- c(angle2, length_tr1)
      #print(paste(k, "out of ", length(angle)))

    }
  
    #find max value and assign its length
    mat_res1 <- c(mat_res[which.max(mat_res[,2]),1], max(mat_res[,2]))
    #plot only the longest line in the segment
    
    #make perpendicular line on the above one
    mat_res.p <- mat_res1[1] + 90
      if (mat_res.p <= 45){
       xval = 500 * tan(mat_res.p*pi/180)
        yval = 500
      } else if (45 < mat_res.p  & mat_res.p  < 135){
        xval = 500
        yval = 500/tan(mat_res.p *pi/180)
      }else if (135 <= mat_res.p ){
       xval = tan(mat_res.p * pi/180)*(-500)
       yval = -500
      } 

    x_coord2 <- c(vec_centr$x+xval, vec_centr$x-xval)
    y_coord2 <- c(vec_centr$y+yval, vec_centr$y-yval)
    
    l2 <- cbind(x_coord2, y_coord2)
    sl2 <- Line(l2)
    s2 <- Lines(list(sl2), ID = "a")
    myLine2 <- SpatialLines(list(s2), proj4string = vec_centr@proj4string)
    #plot(myLine2, add=TRUE)
    
    #cut horizontal line
    transect2 <- gIntersection(vect_vess.p, myLine2)
    plot(transect2, col = "red", add = TRUE )
   
    #calculate the length of the horizontal intersected line
    length_tr2 <- gLength(transect2)
    
    mat_res2 <- c(vec_centr@coords[1,1], vec_centr@coords[1,2], mat_res1, length_tr2)
    print(j)

    fin_res[[j]] <- mat_res2

  }
  
  #assign polygon position (x, y), heading, length and width into an attribute table 
  fin_res.m <- matrix(unlist(fin_res), ncol = 5, byrow = T)
  vect_vessels@data$x <- fin_res.m[,1]
  vect_vessels@data$y <- fin_res.m[,2]
  vect_vessels@data$heading <- fin_res.m[,3]
  vect_vessels@data$length <- fin_res.m[,4]
  vect_vessels@data$width <- fin_res.m[,5]

  #calculate ratio of segments (for vessels between 4-7)
  vect_vessels@data$ratio <- vect_vessels@data$length/vect_vessels@data$width 
  #calculate elipsoidity of segments (for vessels between 0.89 - 0.96)
  vect_vessels@data$elipse <- (vect_vessels@data$length^2-vect_vessels@data$width^2)/(vect_vessels@data$length^2+vect_vessels@data$width^2) 
  
  head(vect_vessels)
  
  #calculate spectral values (mean) for all polygons
  ####time consuming
  system.time(val <- extract(img, vect_vessels, fun=mean)) 
  #If not working: pastecs package has a function of the same name, might call an error
  #to solve this type .rs.unloadPackage("pastecs") in the command line, to uninstall pastecs package
  
  head(val)
  #assign those values as separate columns to the original vessel shapefile
  vect_vessels@data$ch1 <- val[,1]
  vect_vessels@data$ch2 <- val[,2]
  vect_vessels@data$ch3 <- val[,3]
  vect_vessels@data$ch4 <- val[,4]
 
  #head(vect_vessels)
  
  #write vector with all new attributes to a file
  writeOGR(vect_vessels, paste0(outDir, name, "_vector_att"), paste0(name, "_vector_att2"), driver="ESRI Shapefile", overwrite_layer = TRUE)
  
  }

proc.time() - ptm

