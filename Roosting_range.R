### Use of exclusive roosting ranges in disc-winged bats, Thyroptera tricolor
#Authors: Silvia Chaves-Ramírez, Maria Sagot, Hellen Solís Hernández, Gloriana Chaverri

#Install and load the following libraries

library(tidyverse)
library(sf)
library(dplyr)
library(ggplot2)
library(magrittr)
library(sp)
library(rgdal)
library(leaflet)
library(ggthemes)
library(raster)
library(ks)
library(rgeos)
library(ggmap)
library(stringr)
library(maps)
library(mapdata)
library(rpostgis)
library(adehabitatHR)
library(maptools)
library(spatstat)
library(spatialEco)
library(dismo)
library(units)
library(tmap)
library(GISTools)
library(viridis)
library(patchwork)

#DATA

# Read the csv file (should be in your working directory)
points_groups <- read.csv("30UbiG_date.csv", 
                       stringsAsFactors = FALSE) #Contains spatial location data (GPS) of 11 social groups of Thyroptera tricolor

leaves <- read.csv("leaves.csv", 
                   stringsAsFactors = FALSE) #Contains data on spatial locations (GPS) of roosts available 20 samples in the study area

#Groups identified by an ID
colnames(points_groups)
points_groups

colnames(leaves)
leaves

#Convert points_groups object to SpatialPointsDataFrame
# SpatialPointsDataFrame objects don't like missing values
# Remove rows with NA's
points_groups<- points_groups[!is.na(points_groups$x) & !is.na(points_groups$y),]



# create a copy as a 3 column Data frame
points_groups.sp <- points_groups[, c("ID", "x_", "y_")]
points_groups.sp
plot(points_groups.sp$x_, points_groups.sp$y_, asp = 1)

leaves.sp <- leaves[, c("Muestreo", "x", "y")]
leaves.sp
plot(leaves$x, leaves$y, asp = 1)


#Convert ID to Factor
points_groups.sp$ID <- factor(points_groups.sp$ID)

##Convert hoja to Factor
leaves$Hoja <- factor(leaves$Hoja)


#Establish the projection of spatial data
prj <- '+init=epsg:5367'
CRS("+init=epsg:5367") 


summary(points_groups.sp)
view(points_groups.sp)

summary(leaves.sp)
view(leaves.sp)



#Create SpatialPointsDataFrame

spdf_groups <- SpatialPointsDataFrame(coordinates(cbind(points_groups.sp$x_, points_groups.sp$y_)), data = points_groups.sp, proj4string =  CRS(prj))

spdf_leaves <- SpatialPointsDataFrame(coordinates(cbind(leaves$x, leaves$y)), data = leaves, proj4string =  CRS(prj))

groups_utm <- spTransform(spdf_groups, CRS(prj))
leaves_utm <- spTransform(spdf_leaves, CRS(prj))

spdf_leaves_utm <- sp::spTransform(spdf_leaves, sp::CRS("+init=epsg:5367"))


####################################################################################


#Split a shapefile based on attributes (in this case group ID) to make a series of spatial layers 

group_sep <- split(x = spdf_groups, f = spdf_groups$ID, drop = FALSE)

bw <- lapply(group_sep, FUN = function(x){ks::Hlscv(x@coords)})
#create the kernel density estimate. Then in the same function call we’ll convert the kernel density estimate into a raster layer that is spatially explicit

group_kde <-mapply(group_sep,
                  SIMPLIFY = FALSE,
                  FUN = function(x,y){
                    raster(kde(x@coords,h=y))})



##Let’s find the 50% contour. To do this we’ll write a custom function. The function will take the kde and probability as inputs. The below function call getContour will make the 95% KDE for us.

# This code makes a custom function called getContour. 
# Inputs:
#    kde = kernel density estimate
#    prob = probability - default is 0.95

getContour <- function(kde, prob = 0.50){
  # set all values 0 to NA
  kde[kde == 0]<-NA
  # create a vector of raster values
  kde_values <- raster::getValues(kde)
  # sort values 
  sortedValues <- sort(kde_values[!is.na(kde_values)],decreasing = TRUE)
  # find cumulative sum up to ith location
  sums <- cumsum(as.numeric(sortedValues))
  # binary response is value in the probability zone or not
  p <- sum(sums <= prob * sums[length(sums)])
  # Set values in raster to 1 or 0
  kdeprob <- raster::setValues(kde, kde_values >= sortedValues[p])
  # return new kde
  return(kdeprob)
}



group_50kde <- lapply(group_kde,
                     FUN = getContour,prob = 0.50)



group_50poly <- lapply(group_50kde, 
                      FUN = function(x){
                        x[x==0]<-NA
                        y <- rasterToPolygons(x, dissolve = TRUE)
                        return(y)
                      })


group_50poly <- mapply(group_50poly, names(group_50poly), 
                      SIMPLIFY = FALSE,
                      FUN = function(x,y){x@polygons[[1]]@ID <- y
                      return(x)})


group_50poly <- do.call(rbind, group_50poly)


group_50poly$ID <- getSpPPolygonsIDSlots(group_50poly)


groupkde_utm <- sp::spTransform(group_50poly, sp::CRS("+init=epsg:5367"))
view(groupkde_utm)
 plot(groupkde_utm)

#Estimation of the size of roosting range (50%)

rgeos::gArea(groupkde_utm, byid = TRUE)/10000

########

##To export SpatialPolygonsDataFrame layer to gis
class(groupkde_utm)

#writeOGR(obj=groupkde_utm , "gispoly" , dsn="/grupo_1",  layer="grupos", driver="ESRI Shapefile")
########################################################################################

##Create Point-in-Polygon Overlay (para todos los polígonos juntos)

#Create the subsets first (by sampling)
leaves_subset <- subset(spdf_leaves, Muestreo == "M1")


#Generate plots with overlapping points for each sample with all groups
#png(file="1M1.png",
    #width=600, height=350)
#plot(group_50poly, col= group_50poly$ID)
#points(leaves_subset)
#dev.off()


### Graphics with ggplot

## Figure 2. Roosting ranges (50% EDK) for 11 social groups of Thyroptera tricolor (colored polygons). 



#Fortify data frame manually

leaves <- as(spdf_leaves, "data.frame") 
leaves_subset <- as(leaves_subset, "data.frame")  
pg_reales<- as( spdf_groups,"data.frame")
names (pg_reales)[1] = "Group"


#Roosts available
#Plot the leaves available in 20 samples 
hj <- ggplot(leaves) + coord_equal()+
   geom_point(data=as.data.frame(leaves), aes(x,y), 
              color="black", size=0.6, shape=4, alpha = 0.2)+ theme_bw() + theme(panel.background = element_rect( fill = "transparent", size = 1)) 

hj
# Add Plot the roosting ranges of the 11 groups of Thyroptera tricolor 

HJ <-hj + geom_polygon(data=group_50poly, aes(x = long, y = lat, group = id, fill= id), alpha=0.3)+
  labs( x = "Longitude", y = "Latitude")+ coord_equal()+  theme_bw() 
HJ

#Add spatial points of each group of the same color as the corresponding roosting range
Q  <- HJ + geom_point(data=as.data.frame(pg_reales),  aes(x =x_, y = y_, group = Group, fill= Group, color = Group), alpha=0.4)+ labs( x = "Longitude", y = "Latitude")+ theme_bw() #+
Q
Q1 <- Q + scale_colour_manual(
  values = cols,
  aesthetics = c("colour", "fill"))
Q1
 

