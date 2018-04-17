#' @title Vessel classification module from Sentinel-2 detection results (batch mode)
#' @description Algorithm for vessel classification which takes vessel detection polygons as an input and returns classified vessels 
# based on their attributes using a decision tree algorithm (in a batch mode). 
#' @author Ursa Kanjir
#' @date "'r format(Sys.Date())'"
#' @return a shapefile classified polygons that represent vessels and other false alarms 

require(rgeos)
require(sp)
library(rpart)
library(rpart.plot)

ptm <- proc.time()

mainDir <- "Select your directory with shp detection results from prior step" 
setwd(mainDir)
trainDir <- "Select your directory where your training samples are stored"
outDir <- "Select directory where you wanna save your results"


class_List <- list.files(mainDir, pattern = "*_vector_att.shp", recursive = TRUE, full.names = TRUE)
class_List
# With this algorithm we obtain classification trees from test data (classification tree)
#clas_tree results are calculated for each image with test data separately

for (i in 1:length(class_List)) {
 
  name_full <- sub(".*/", "", class_List[i])
  name <- substr(name_full, 1, nchar(name_full) - 28)
  name
  
  object.shp <- readOGR(class_List[i])
  head(object.shp)

  # Polygin has to obtain column called "class" in which there are around 40% of train data approximately equally distributed classes
  train.obj <- object.shp[!is.na(object.shp$class),]
  train.obj.df <- data.frame(train.obj)
  
  # those are shapefiles with no class values (test data) 
  test.obj <- object.shp[is.na(object.shp$class),] 
  
  table(train.obj$class) #How many examples there is in each class
  class_tree <- rpart(class~ area+length+width+ratio+ch1+ch2+ch3+ch4, data = train.obj, method = "class")
  
  #visualise and save classification tree as image
  #plot(test_model)
  #text(test_model, cex=0.8)
  png(filename = paste0(outDir, name, "Classification_tree.png"), width=450, height = 450)
  rpart.plot(class_tree, type=3, extra = 101)    #rpart.plot package for better visualisation of decision trees
  dev.off()
  
  #prp(test_model, type = 1, extra = 1, branch = 1) #yet another version of visualisation
  pred.train <- predict(class_tree, train.obj, type="class")
  
  #Assign unclassified segments class accordingt to classification tree rules
  pred.test <- predict(class_tree, test.obj, type="class")
  
  #confusion matrix - how good is a classification based on visually selected classes
  table(train.obj.df[,15], pred.train)
  
  object.shp[is.na(object.shp$class), "class"]<- pred.test
  writeOGR(object.shp, dsn = paste0(outDir, name, "vess_class"),layer = paste0(name, "vess_class"), driver = "ESRI Shapefile", overwrite_layer = T)  
  } 

plotcp(class_tree)

