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
  
  r2<- rasterRescale(img) #this step can take a lot of memory in the case of a large image

  #sum of all the bands in a certain pixel
  img_tog <- sum(r2)
  #plot(img_tog)
  #writeRaster(img_tog, filename = paste0(outDir, name, "_testS2_okt.tif"), format="GTiff", overwrite=TRUE)
  
  #cellStat - statistics of the cells of a raster object!
  img.mean <- cellStats(img_tog, 'mean', na.rm=TRUE) #izračunamo mean vsem pikslom, je memory safe
  #drugi način (ki ni kao memory safe) je mean(values(img_tog), na.rm = T)
  img.sd <- cellStats(img_tog, 'sd', na.rm=TRUE)
  
  img.meja <- img.mean + 0.55*img.sd #0.55 is pragmatically selected threshold for S-2 images

  #Sedaj stvar deluje za vse senzorje enako, saj išče mejo samodejno in je drugačna za vsak img posebej
  vessels <- (vess_indx > 0) | (img_tog > img.meja) #2.4 #3.3 #2.78
  #stestiraj še za ostale senorje, ampak po moje bi moralo delati!!!!
  #plot(vessels)
  #writeRaster(vessels, filename = paste0(outDir, name, "_vessels.tif"), format="GTiff", overwrite=TRUE)
  
  rm(r2, vess_indx) #ta dva sta velika in jo zbrišem, da ne kradeta spomin
  
  #Izboljšamo masko plovil (popraviti LUKNJE v segmentih), da ne pride do napak pri kasnejših fazah!!!
  vessels_smooth <- focal(vessels, w = matrix(1,3,3), fun = max, pad = TRUE)
  vessels_smooth_fin <- focal(vessels_smooth, w = matrix(1,3,3), fun = min, pad = TRUE)
  #plot(vessels_smooth_fin)
  
  writeRaster(vessels_smooth_fin, filename = paste0(outDir, name, "_poss_vessel2.tif"), format="GTiff", overwrite=TRUE)

  ################Obtain vector polygons from raster of possible vectors
   
  #FIRST OPTION _ time consuming calculation from raster to vector (NE DELA za prevelike imgje - S2 npr)
  #current.time <- Sys.time()
  #vect_vess <- rasterToPolygons(vessels, dissolve = TRUE)
  #Sys.time() - current.time
  
  #SECOND OPTION - short calculation - BUT needs gdal_polygonize.py on your computer
  #funkcija potegnjena iz https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
  system.time(vect_vessels <- gdal_polygonizeR(vessels_smooth_fin, outshape = paste0(outDir, name, "_vector_temp")))
  # za sentinel2 velik file poligonizira v 15ih sekundah! 
  #vektorski in rasterski sloj imata sedaj enako projekcijo in torej extent()
  #POMEMBNO! če datoteka z istim imenom že obstaja bo koda prekinjena! briši sproti!
  
  #using a zero width buffer cleans up many topology problems in R (da se ne pojavljajo napake!)
  vect_vessels <- gBuffer(vect_vessels, byid=TRUE, width=0)
  
  #calculate area of shape polygons 
  system.time(vess_area <- gArea(vect_vessels, byid = TRUE))
  #assign area values to the polygon you are interested in
  vect_vessels@data$area <- vess_area
  
  #define minimal and maximal value for different sensors
  #glej tisti del o najmanjših možnih površinah za zaznavo enega konkretnega objekta v preglednem članku
  
  # LOČJIVOST JE RAZLIČNA ZA VSAK SENZOR POSEBEJ! DOLOČI ŠE ZA OSTALE SENZORJE!!!!
  #različno poimenovanje glede na senozor
  if (grepl("S2", lm_List[i])) {  #za vse  S2 posnetke
    res <- 10 #za S- 2
  } else if (grepl("ge_|wv2_", lm_List[i])){    #za vse ostale VHR senzorje
    res <- 0.5
  } else if (grepl("ik_", lm_List[i])){
    res <- 1
  } else if (grepl("qcbd_", lm_List[i])){
    res <- 0.75
  } else print("Used sensor resolution is not known!")
  
  print(paste0("Ločljivost obravnavanega senzorja je ", res, "m"))
  
  if (grepl("S2", lm_List[i])) {  #za vse  S2 posnetke
    min.val <- (res^2 * 3) #m2 (po Bannister and Neyland, 2015)
    max.val <- 30000  #to potrebno še stestirat!
  } else if (grepl("ge_|wv2_", lm_List[i])){    #za vse ostale VHR senzorje
    min.val <- (res^2 * 5) #to sem določila empirično
    max.val <- 400 #to tudi
  } else if (grepl("ik_", lm_List[i])){
    min.val <- (res^2 * 4) #to sem določila empirično
    max.val <- 1000
  } else if (grepl("qcbd_", lm_List[i])){
    min.val <- (res^2 * 4) #to sem določila empirično
    max.val <- 750
  } else print("Used sensor is not known!")
  print(paste0("minimalna površina plovil je ", min.val, "m2"))
  print(paste0("maximalna površina plovil je ", max.val, "m2"))
  #max.val mora biti manjše, kot je površina največje ladje če ne je preveč napak 
  #tudi če zazna velike ladje, jih zazna delno po navadi, ker so njihove površine zelo heterogene
  #pri prevelikih objektih pa je itak toliko napak, da ni vredno
  
  #remove all those polygons with DN value = 0 (potrebujemo samo tiste z value = 1 - possible ship)
  #remove all the polygons that are smaller tham ...m2 and bigger than ...m2 - ODVISNO OD SENZORJA IN NJEGOVE LOČLJIVOSTI
  vect_vessels <- vect_vessels[(vect_vessels@data$area > min.val &
                                  vect_vessels@data$area < max.val & 
                                  vect_vessels@data$DN == 1),]
  #writeOGR(vect_vessels, paste0(outDir, name, "_vector_vmesni"), paste0(name, "_vector_vmesni"), driver="ESRI Shapefile", overwrite_layer = TRUE)
  
  print(paste0("Image ", i, " includes ", length(vect_vessels), " possible vessel objects."))
  
  #set angle sequence
  set_angle <- 10 # set angle sequence - for how much it should move between 0-180
  angle <- seq (0, 179, set_angle)
  #set empty matrix where calculated angles will be later stored
  fin_res <- list()
  mat_res <- matrix(nrow=length(angle), ncol=2) # matrix with the results of distances of every angle
  colnames(mat_res) <- c("heading", "length")
  mat_res2 <- matrix(nrow=1, ncol=5)
  colnames(mat_res2) <- c("x", "y", "heading", "length", "width")
  
  #create a buffer geometry - popravi poligone (v prvi različici so se mi podvajali)
  #polygone1 <- gBuffer(vect_vessels, byid=TRUE, width=0)
  #polygone2 <- gBuffer(vect_vessels, byid=TRUE, width=0)
  #vec_ves <- gIntersection(Polygone1, Polygone2, byid=TRUE)

  #add ID column to the vector table (kako preprosto!)
  vect_vessels@data$idloop <- seq(1,length(vect_vessels),1)
  #head(vect_vessels@data)
  #rm(i)
  
  for (j in 1:length(vect_vessels)){
  #for (j in 17:26){ # ZA TEST!   
    #j=168
    vect_vess.p <- vect_vessels[(vect_vessels@data$idloop%in%j),]  #so it selects only one polygon out of many
    plot(vect_vess.p)
    
    #Calculate vessel polygon width and beam
    #add centroid to the layer
    vec_centr <- gCentroid(vect_vess.p, byid = TRUE) #vec_ves_s
    points(vec_centr)  #plot points
    
    #k = 2
    #angle = 25
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
      #dolga pot
      #x_length <- transect@bbox[1,1]-transect@bbox[1,2]
      #y_length <- transect@bbox[2,1]-transect@bbox[2,2]
      #length_tr1 <- abs(sqrt(x_length^2 + y_length^2))
      
      #kratka pot
      length_tr1 <- gLength(transect)
      
      mat_res[k,] <- c(angle2, length_tr1)
      #print(paste(k, "out of ", length(angle)))

    }
  
    #find max value and assign its length
    mat_res1 <- c(mat_res[which.max(mat_res[,2]),1], max(mat_res[,2]))
    #plot only the longest line in the segment
    #iz Centroida (ki ima znane koordinate) bi s pomočjo razdalje in kota izračunala robne koordinate
    
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
    #dolga pot
    #x_length <- transect2@bbox[1,1]-transect2@bbox[1,2]
    #y_length <- transect2@bbox[2,1]-transect2@bbox[2,2]
    #length_tr1 <- abs(sqrt(x_length^2 + y_length^2))
    #kratka pot, ne vedno pravilna, zakaj???
    length_tr2 <- gLength(transect2)
    
    mat_res2 <- c(vec_centr@coords[1,1], vec_centr@coords[1,2], mat_res1, length_tr2)
    print(j)

    fin_res[[j]] <- mat_res2

  }
  
  #pripiši heading, length in width v atributno tabelo 
  #heading
  #vect_vessels@data$heading <- lapply(fin_res, "[", "heading")
  #to je super, če bi ohranila list, ampak list noče pripisati v att tabelo kasneje v fazi writeOGR
 
  #x, y, heading, length and width
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
  ####to dela RES dolgo!  Kako bi se to dalo pohitriti? glej tiling mogoče? paralelizacija!
  system.time(val <- extract(img, vect_vessels, fun=mean)) 
  #ČE NE DELA: pastecs package ima isto ime ene funkcije, zato lahko pride tukaj do napake!
  # v tem primeru samo odtipkaj .rs.unloadPackage("pastecs"), da odinštaliraš tisti paket
  
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

