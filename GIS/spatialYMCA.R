######################### %%%%%%%%%%%%% ###########################

#### HANDLING SPATIAL DATA IN R AND BASIC SPATIAL CALCULATIONS ####
#
#

######################## %%%%%%%%%%%%%%% ##########################

# In spatial analysis, special softwares have been developped to ease the use of spatial data as GRASS, QGis or ArcGis
# Those softwares are really good and allow very good computations on spatial files
# sometimes more quickly and easily than in R
# However, the fact that inside those software we can use functions already designed, make us wondering about what exact computation is behind
# and how to reproduce our analyses.

# For this, if you know how to code in Python language or use the console you can track and/or code exactly what you want.
# However, in Ecology and other Conservation Biology fields, the use of R is more than common and allow for reproducible examples, help, creation of functions and tracking what we have done.
# Moreover, it can allow to change easily some stuff such as the input data and/or some computation details.

# This is why the use of R for handling and analysing spatial data is more and more common nowdays
# and also why it should be important to know how to use it.

# As such, we are going through this script to do some example of handling, computing and mapping some spatial data in R
# using different functions and packages and compare it with QGis (as it is free software).

## remark: You can move inside the script to any parts by selecting the part you want on the bootom left of this script






                            #### %%%%%%%%%%%%% ####


                      #### A - Loading spatial data ####


                           #### %%%%%%%%%%%%% ####



# -----------------------------------------------------------------------
    # 1° Loading spatial data VECTOR LAYERS in R ####


# easy example and comparison with **sf package** ----

library(rgdal)

library(sf) # Linking to GEOS 3.6.1, GDAL 2.2.3, proj.4 4.9.3


# to collect French department I used a governmental website
# that actually provide OpenStreet Map data (but from the national cadastre)
#https://www.data.gouv.fr/fr/datasets/contours-des-departements-francais-issus-d-openstreetmap/

setwd("mypath")
# replace here your path for the file
#setwd("C:/Users/cazam/Desktop/YMCA sp")

# open with rgdal 
departement <- readOGR("./departements-20180101.shp",
                       layer="departements-20180101")

departement
  # we have several information here
class(departement)
  # what type of file is departement


# open with sf (always begging with st_)
depart_sf <- st_read("./departements-20180101.shp")
# replace here your path for the file

  # what we have in depart_sf is all the information
  # with **sf package** you open a spatial file as a data frame
  # the geometry is actually another column
  # thus, if you work with sf you must go on during all the analysis
  # as R does not recognize the equivalence between a sf file and a gdal file
depart_sf
  # as before you have several information directly


# ---------------------------------------------------------------------------
      # 2° Loading spatial data RASTER LAYERS in R ####

# you will need the raster package
# install.packages("raster")

# download CLC12 from https://land.copernicus.eu/pan-european/corine-land-cover/clc-2012?tab=download in the raster type
library(raster)
clc <- raster("./g250_clc12_V18_5.tif")
clc
  # several information






                                #### %%%%%%%%%%%%% ####


                      #### B - Handling projection system ####


                                #### %%%%%%%%%%%%% ####



# ----------------------------------------------------------------------------
    # 1° How to reach the projection system ####

proj4string(departement)
depart_sf
#proj4string(depart_sf)
# Error in (function (classes, fdef, mtable)  : 
# unable to find an inherited method for function ‘proj4string’ for signature ‘"sf"’
projection(clc)


# ----------------------------------------------------------------------------
    # 2° How to change the projection system ####

#plot(clc)
#plot(departement,add=T)
  # not working because they are not in the same projection
proj4string(clc) == proj4string(departement)


# it is easier to change the projection of a vector layer rather than a raster layer
departement_laea <- spTransform(departement,CRS(proj4string(clc)))
#plot(departement_laea,add=T) # it takes a little while but it works
#departement_laea <- spTransform(departement,CRS("+init=epsg:3857"))
  # possibility of using a reference number for each CRS

# with sf
dpsf_laea <- depart_sf %>% 
  st_transform(crs = st_crs(clc)) 
dpsf_laea


# if you really want to change the projection of a raster (function QGis::Warp)
?projectRaster
#clc_wgs <- projectRaster(clc,crs=proj4string(departement))
#proj4string(clc_wgs) == proj4string(departement)
#plot(clc_wgs)




                        #### %%%%%%%%%%%%% ####


                      #### C - DATA REQUEST  ####


                        #### %%%%%%%%%%%%% ####



# ----------------------------------------------------------------------------
    # 1° VECTOR LAYERS ####

      # a - How to reach the attribute table ----


# in order to ease the next steps
# we are going to use the LAEA projection (Europe)

# with gdal
departement_laea@data
head(departement_laea@data)
head(departement_laea@polygons) # you can access also the geometry


# with sf package you can examine the structure
  # it works with the tidyverse package
  library(tidyverse)
# and tidyverse can work also with gdal
glimpse(depart_sf)
names(departement_laea@data)
nrow(departement_laea)
glimpse(departement_laea)



########
# STOP #
########

# if you only need the data for analysis or calculations you can extract the data frame and take it from there
dat <- departement_laea@data
class(dat)
# and the coordinates if needed (for example in some models you have to use the locations)
head(coordinates(departement_laea)) # centroids of polygons
nrow(coordinates(departement_laea))
nrow(departement_laea)


# TO CHECK FOR CENTROIDS ####
#plot(departement_laea)
#library(rgeos)
#cp <- gCentroid(departement_laea)
#cp
#coordinates(cp)
#plot(cp,add=T,col="red") # we created the centroid of our vector file
# if we want the centroid for each polygon
#cp.each <- gCentroid(departement_laea,byid = T)
#plot(cp.each,add=T,col="green",pch=20)     
#coordinates(cp.each)
#nrow(coordinates(cp.each))
#cbind(coordinates(cp.each),coordinates(departement_laea))
  # some small differences





      # b - Select variables ----


# with gdal it works as a data frame
departement.study <- departement_laea[,c("nom")]
class(departement.study)
head(departement.study@data)

# with sf and tidyverse
depart.site <- dpsf_laea %>%
  select('Departement' = nom)
depart.site

# let's select ONLY one department:
departement.study.H <- departement.study[which(departement.study@data$nom == "Creuse"),]
departement.study.H
#plot(departement.study.H)

# with sf and tidyverse
depart.site.H <- depart.site %>%
  filter(Departement == "Creuse")
depart.site.H
#plot(depart.site.H)




      # c - First easy and quick (and ugly) map ----

# WITHS SF ####
# set the theme for ggplot (the library has been opened with tidyverse)
theme_set(theme_minimal())

# with sf and ggplot2
depart.site %>%
  ggplot() +
  geom_sf() +
  labs(title = 'Nos départements français',
       subtitle = 'Source: gouv, OpenStreet Map')
  # can take time to display

# WITH GDAL (SP) and ggplot2 ####
plot(departement.study)

# with gdal and ggplot2 ## Romain Lorellière for fortify actions ##
fptot <- fortify(departement.study) # can take a little time
data.tot <- departement.study@data
data.tot$id <- rownames(data.tot)
fptot <- merge(fptot,data.tot,by="id") # can take a little time

head(fptot)

plotgd <- ggplot() + 
  geom_polygon(data=fptot,aes(x=long,y=lat,group=group,fill=nom)) + 
  coord_equal() + 
  coord_quickmap() +
  theme_bw() +
  theme(legend.position='none')
plotgd

# of course once you have this basis you can add and remove what you want with the common ggplot2 commands


# ZOOM EXAMPLE ####
# just if we want to zoom in on the metropole

plotgd.metro <- ggplot() + 
  geom_polygon(data=fptot,aes(x=long,y=lat,group=group,fill=nom)) + 
  coord_equal() +
  coord_quickmap(xlim = c(3e+06,4.5e+06),ylim = c(2e+06,3.5e+06)) +
  theme_bw() +
  theme(legend.position='none')
plotgd.metro


# see the gif on this page :
# https://statnmap.com/fr/2018-07-14-initiation-a-la-cartographie-avec-sf-et-compagnie/
# to see the difference between geographical and planar coordinates



      # d - Some quick and easy calculations on vector layers ----


# Select Loire country side ####

# 44, 85, 49, 72, 53 departments ID number

head(departement_laea@data)
p.Loire <- departement_laea[which(departement_laea@data$code_insee == 44 |
                                    departement_laea@data$code_insee == 85 |
                                    departement_laea@data$code_insee == 49 |
                                    departement_laea@data$code_insee == 72 |
                                    departement_laea@data$code_insee == 53),]
plot(p.Loire)


# with sf
psfL <- dpsf_laea %>%
  filter(code_insee == 44 | code_insee == 85 |code_insee == 49 |code_insee == 72 |code_insee == 53) %>%
  select('Departement' = nom) 
psfL

psfL %>%
  ggplot() +
  geom_sf() +
  labs(title = 'Les pays de la Loire',
       subtitle = 'Source: gouv')




# Let's add some data in those departments ####
# Number of castles by departments (random data) 

p.Loire$nbChato <- c(12,2,6,4,7)
head(p.Loire@data)

# with gdal and ggplot2 ## Romain Lorellière for fortify actions ##
fptot <- fortify(p.Loire) # can take a little time
data.tot <- p.Loire@data
data.tot$id <- rownames(data.tot)
fptot <- merge(fptot,data.tot,by="id") # can take a little time

head(fptot)

plotchato <- ggplot() + 
  geom_polygon(data=fptot,aes(x=long,y=lat,group=group,fill=nbChato)) + 
  coord_equal() + 
  coord_quickmap() +
  theme_bw()
plotchato

  ##[RL] this methods can be fastest than geom_sf() but be careful with the map that relies on the graph dimension and the legend size (which is not the case with geom_sf()


# with sf
Loire_chato <- cbind(psfL, c(12,2,6,4,7))
Loire_chato
# rename column
Loire_chato.re <- rename(Loire_chato,
                         Chateaux="c.12..2..6..4..7.")
Loire_chato.re

Loire_chato.re %>%
  ggplot() +
  geom_sf(aes(fill = Chateaux)) +
  labs(title = 'Nos châteaux de la Loire')

# also with sf you can group department into region directly for example. There are several other handy functions like this one in sf relating to tidyverse which allows quick calculations on data frame.


# ----------------------------------------------------------------------------
    # 2° RASTER LAYERS ####

      # a - How to reach the values ----


clc
  # it is a big raster
  # thus as example we are going to subset a spatial area of this raster to ease and speed the next steps

?crop
bbox(p.Loire) # spatial extent on which we are going to subset the raster
proj4string(p.Loire) == proj4string(clc)


clc.Loire <- crop(clc,# raster to crop
                  extent(bbox(p.Loire))) # extent on which to subset it has  to be a extent type object
plot(clc.Loire)
plot(p.Loire,add=T)
unique(getValues(clc.Loire)) 
# normally we get all the number informing about the class of landscape cover in this area


      # b - How to change the values ----

# we want only the forest and we do not care about the type natural grassland
# in the legend of clc natural grassland == 26

new.clc.grass <- clc.Loire
new.clc.grass[] <- 0 #☺ empty raster
w2 <- which(clc.Loire[]==26)
new.clc.grass[w2] <- 1
w3 <- which(is.na(clc.Loire[]))
new.clc.grass[w3] <- NA
plot(new.clc.grass) # binary forest cover
unique(getValues(new.clc.grass))
new.clc.grass
extent(new.clc.grass) == extent(clc.Loire)



      # c - Ugly quick raster maps ----

# Transform rasters as data frame 
library(ggplot2)
# this function can be found at :
# https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r

gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

  #♠ change the raster in a frame for ggplot2 to handle
gplot_Loireclc <- gplot_data(clc.Loire)

ggplot() +
  geom_tile(data =gplot_Loireclc,aes(x=x,y=y,fill=Value)) +
  scale_fill_gradient("Value",
                      low = 'yellow', high = 'blue',
                      na.value = NA) +
  coord_quickmap() +
  coord_equal()


# with SF ####
extent(gplot_Loireclc)
extent(psfL)
ggplot() +
  geom_tile(data =gplot_Loireclc,aes(x=x,y=y,fill=Value)) +
  scale_fill_gradient("Value",
                      low = 'yellow', high = 'blue',
                      na.value = NA) +
  geom_sf(data=psfL,fill=NA,colour="black",size=0.2)  


# please note the difference background between the plot in sf and the classic ggplot. The background in sf displays the wgs projection (source projection) but gives the same picture as the classic ggplot




                                #### %%%%%%%%%%%%%%  ####

                          
                #### D - SPATIAL REQUEST and COMPUTATIONS ####


                               #### %%%%%%%%%%%%%%  ####

# ----------------------------------------------------------------------------
    # 1° Spatial computations on vectors ####

# compute how many points inside one polygon
# we are going to compute false data set of points

bbox(p.Loire)
p.Loire
p <- data.frame(runif(50,min=bbox(p.Loire)[1,1],max=bbox(p.Loire)[1,2]),runif(50,min=bbox(p.Loire)[2,1],max=bbox(p.Loire)[2,2]))
names(p)[1] <- "X"
names(p)[2] <- "Y"
coordinates(p) <- ~ X + Y
proj4string(p) <- CRS(proj4string(p.Loire))
class(p)

plot(p.Loire)
plot(p,add=T)

# Select only the points inside a Loire's department
library(rgeos)
u <- gIntersects(p,p.Loire,byid=T)
u
class(u)
  # u provides a matrix of T/F. it is true when the point j (column) is inside a department i (row)
nrow(u)
ncol(u)
nrow(p.Loire)
length(p)

stock <- c()
  #☺ we create a vector into which we are going to stock the ID of each point contained in each department
for(i in 1:nrow(u)){
  num.row.pix <- levels(as.factor(which(u[i,] == T)))
    # here we select for the point that are contained u[i,]==T
    # and we extract the ID levels(as.factor(which(u[i,]==T)))
  stock <- c(stock,num.row.pix)
    # thus for each row (i or department) we will stock the ID of each point contained in the department i
}
stock
  # thus, here we have the ID of all the points that are contained into the Loire's departments


pointloire <- p[as.numeric(stock),]
  # then, here we select the row (the point) that correspond to the ID extracted with the loop
length(pointloire)
length(stock)
plot(pointloire,add=T,lwd=2,cex=3,col="red",pch=20)

# * be aware that we might not have the same plot because we created our points from a random process *



# Now let's assign the number of points to each polygon containing them
# (we could have done both in one step)

u
nrow(u) # number of polygons
ncol(u) # number of points
length(p)

stock <- numeric(nrow(p.Loire))

for(i in 1:nrow(p.Loire)){
  num.pts <- levels(as.factor(which(u[i,] == T)))
    # as before, here we select for the point that are contained u[i,]==T
    # and we extract the ID levels(as.factor(which(u[i,]==T)))
  if(length(num.pts) == 0L){ # if the point is not contained
    stock[i] <- "0"
  }else{ # if there are points contained
      # we do not care about the ID now
      # but about how many IDs are contained
      # how many number of points
    stock[i] <- length(num.pts)
  }
}
stock
head(p.Loire@data)
p.Loire@data$Pts <- stock
fptot <- fortify(p.Loire) # can take a little time
data.tot <- p.Loire@data
data.tot$id <- rownames(data.tot)
fptot <- merge(fptot,data.tot,by="id") # can take a little time

head(fptot)

plotpts <- ggplot() + 
  geom_polygon(data=fptot,aes(x=long,y=lat,group=group,fill=Pts)) + 
  coord_equal() + 
  coord_quickmap() +
  theme_bw()
plotpts
  

# rather than gIntersects ou can also take a look at gContains in the same library


# ----------------------------------------------------------------------------
    # 2° Buffers around points around polygons ####

# let's take our points contained in the Loire departments
pointloire
plot(pointloire)
proj4string(pointloire)

?gBuffer
bufloirepts <- gBuffer(pointloire,byid=T,width=10000)

plot(bufloirepts,col="green")
plot(pointloire,add=T)

# around our departments
# it allows you in computations later in the script to avoid edge effect in some computations such as proportion of landscape in a buffer etc
p.Loire
bufloire <- gBuffer(p.Loire,byid=T,width=5000) # can take a little time

plot(bufloire,col="grey")
plot(p.Loire,add=T,col="light blue")




# ----------------------------------------------------------------------------
    # 3° Request on raster layer ####

# let's refresh
# before line 376
#• we cropped our raster on the extent of the Loire region

# now to ease our computations we would like to select only the raster inside the polygons and not just on the extent
# clc.Loire <- crop(clc,extent(bbox(p.Loire)))
# moreover, to avoid some edge effects later
# we are going to select a larger part than just the polygon
# thus the buffer of the polygons


clc.Loirebuf <- mask(crop(clc,extent(bufloire)), # raster cropped
                     bufloire) # buffer on which we subset our raster
clc.Loirebuf
extent(bufloire) == extent(clc.Loirebuf)
  #↓extent(clc.Loirebuf) <- extent(bufloire)
  # if we do that we change the resolution of our raster
  # the fact that the raster is in cells (that can differ in size) it creates a litlle difference in extent
plot(clc.Loirebuf)
plot(bufloire,add=T)
plot(p.Loire,add=T)

# ----------------------------------------------------------------------------
    # 4° Create a raster layer from a data frame ####

# from our points, we want to rasterize the results
# on our departments
plot(pointloire)

# create a raster ----
# create the basic raster on which we are going to create the point density
?raster
basic.raster <- raster(extent(bbox(bufloire)))
basic.raster
projection(basic.raster) <- CRS(proj4string(clc.Loirebuf))
basic.raster
# let's create the resolution at 10km square (for visibility of results)
res(basic.raster) <- 10000
basic.raster
basic.raster[] <- 0
basic.raster <- crop(basic.raster,extent(bbox(bufloire)))
basic.raster <- mask(basic.raster,bufloire)

extent(bufloire) == extent(basic.raster)

# count of observations ----
occ.ras <- rasterize(coordinates(pointloire),basic.raster,fun="count",field = 0,background=NA)
  # cells with 1 and 2
plot(occ.ras)
unique(getValues(occ.ras))

# change values > 1 for 1s to create a binary raster ----
occ.ras[] <- ifelse(occ.ras[]>1,1,occ.ras[])
plot(occ.ras)
  # value only 1
unique(getValues(occ.ras))

# b) change NAs for 0s ----
occ.ras[] <- ifelse(is.na(occ.ras[]),0,occ.ras[])
plot(occ.ras)
  # values 0 and 1
unique(getValues(occ.ras))

# c) change counts out of the area of interests to NAS ----
#occ.ras <- crop(occ.ras,basic.raster)
occ.ras <- mask(occ.ras,basic.raster)
plot(occ.ras)

unique(getValues(occ.ras))
# we have NAs
# wa have 0s
# we have 1s
extent(occ.ras) == extent(basic.raster)


########
# STOP #
########

# IN CASE WE WANT TO STOP THERE AND SAVE OUR FILES

writeRaster(occ.ras,
            #overwrite=T,
            "./Raster_binary_Loire_LAEA.tif")
writeOGR(bufloire,
         #overwrite=T,
         "./Buffer_Depart_Loire_LAEA.shp",
         layer="Buffer_Depart_Loire_LAEA",
         driver="ESRI Shapefile")





                      #### %%%%%%%%%%%%%%  ####


                #### E - SPATIAL COMPUTATIONS ####


                    #### %%%%%%%%%%%%%%  ####

# ----------------------------------------------------------------------------
    # 1° Inside vector layer or between vector layers ####

      # a) Distance neighbour and spatial correlation test ----

xy <- coordinates(pointloire)
xy

library(spdep)

?dnearneigh() # distance-based nearest neighbour
ahah <- dnearneigh(xy,d1=0,# from 0
                   d2=25000) # to 25 km
# create a neighbour type 
ahah
class(ahah)
plot(xy)
plot(ahah,xy)
 
wm2 <- nb2listw(ahah, style='B',zero.policy=TRUE) # create a weighted 0/1 neighbour matrix
class(wm2)

# let's add data to our points to do a spatial autocorrelation test
length(pointloire)
vec <- c(rnorm(length(pointloire),0,1))
class(pointloire)
pointloire$vec <- vec
pointloire@data
class(pointloire)

# Moran's test as spatial autocorrelation test ####
test <- moran.test(pointloire@data$vec, wm2, randomisation=T,zero.policy = TRUE)
test



       # b) Create distances between points and select for nearest neighbour ----

# with rgeos library
?gDistance
val <- gDistance(pointloire,byid=T)
val
summary(val)

# nearest neighbour nnb it is like dnearneigh ou knearneigh 
stock <- numeric(nrow(pointloire))
for(i in 1:nrow(pointloire)){
  valu <- val[i,]
  valu <- valu[which(valu != 0)]
  stock[i] <- min(valu)
}
stock

# elegant way ;)
sy <- knn2nb(knearneigh(xy,k=1)) # k nearest neighbour here k= 1
sy
plot(sy,xy)
# if you want the distances
dsts <- unlist(nbdists(sy,xy))
dsts
cbind(dsts,stock) # we find the same thing
summary(dsts)


# for example if we want 2 and then 4 nearest neighbours for each points
sy2 <- knn2nb(knearneigh(xy,k=2))
sy2
plot(sy2,xy)

sy4 <- knn2nb(knearneigh(xy,k=4))
sy4
plot(sy4,xy)

# please see the vignette of the spdep package:
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

# this type of analysis can be done on polygons using centroids
# but also distance between boudaries so linear type objects
# this can also be done between two vector layers

# as such we can calculate the distance between centroid of departments and points
val <- gDistance(pointloire,p.Loire,byid=T)
# thus here for each point (column) or each polygon (row) we can extract the nearest centroid or nearest point


# ----------------------------------------------------------------------------
    # 2° Inside raster layer or between raster layers ####

      # a) Distance raster ----

# let's refresh
# line 388
#◘ we created a binary raster of natural grassland
# we want only the forest now and we do not care about the type of forest
# in the legend of clc forest == 23.24.25
# and on the buffer of the departments if later we want to do calculations

# to ease and speed

clc.Loire2 <- mask(crop(clc,# raster to crop
                  extent(bbox(bufloire))),bufloire)
plot(clc.Loire2)

new.clc.foret2 <- clc.Loire2
new.clc.foret2[] <- 0 #☺ empty raster
w2 <- which(clc.Loire2[]==23 | clc.Loire2[]==24 | clc.Loire2[]==25)
new.clc.foret2[w2] <- 1
w3 <- which(is.na(clc.Loire2[]))
new.clc.foret2[w3] <- NA
plot(new.clc.foret2) # binary forest cover
unique(getValues(new.clc.foret2))

# prepare the raster to compute distances to forest
val <- getValues(new.clc.foret2)
unique(val)
val
new.clc.foret2[new.clc.foret2==0] <- NA
unique(values(new.clc.foret2))
class(as.numeric(values(new.clc.foret2)))

# get numeric values to be taken into account in the distance function

values(new.clc.foret2) <- as.numeric(values(new.clc.foret2))
plot(new.clc.foret2)

extent(new.clc.foret2) == extent(clc.Loire2)

# Change resolution to visualize better ####
?aggregate
d.foret <- aggregate(new.clc.foret2,fact=8)
extent(d.foret) == extent(new.clc.foret2)
plot(d.foret)
d.foret

#  distance to nearest cell that is not NA
?distance
distah <- distance(d.foret)
distah
val <- getValues(d.foret)
summary(val)
plot(distah)

distah.vrai <- mask(distah,bufloire)
plot(distah.vrai) # be carefull of edge effect in the sea here ! :)



      # b) Proportion raster ----

# prepare the forest raster as 0/1 and NAs

#  change NAs for 0s
d.foret[] <- ifelse(is.na(d.foret[]),0,d.foret[])
plot(d.foret)
unique(getValues(d.foret))

#change counts out of the area of interests to NAS 
d.foret <- crop(d.foret,bufloire)
d.foret <- mask(d.foret,bufloire)
plot(d.foret)

unique(getValues(d.foret))
# we have NAs
# wa have 0s
# we have 1s

# create weighted matrix
?focalWeight
fw<-focalWeight(d.foret,20000, "circle") 
# creates circular filter with a radius of 400m
# ce qui Ã©quivaut Ã  3cell*3cells soit 9 cells
##fw
# or
mat <- matrix(1, nrow = 3, ncol = 3)

# create proportion raster
?focal
ForestCover_prop <- focal(d.foret,w=fw,fun="sum",na.rm=T,pad=TRUE,padValue=F)  
# to compare
#ForestCover_prop <- focal(ra.new,w=fw,fun="sum")
ForestCover_prop.mat <- focal(d.foret,w=mat,fun="mean",na.rm=T,pad=TRUE,padValue=F)  
unique(getValues(ForestCover_prop))
# Attribute the proportion of forest cover
# in a radius of 20km for each cell
# thus between 0 and 1
plot(ForestCover_prop)
ForestCover_prop <- mask(ForestCover_prop,d.foret)
plot(ForestCover_prop)
plot(ForestCover_prop.mat)
ForestCover_prop
extent(ForestCover_prop) == extent(d.foret)

# be aware that gradient can be also evaluate with those kind of approaches




                      #### %%%%%%%%%%%%%%  ####


            #### F - DATA FRAME FOR ANALYSIS FROM RASTERS ####


                      #### %%%%%%%%%%%%%%  ####


#

# let's say you want to do an analysis from all those rasters we did before
# you need to be carefull about extent, resolution and projection
# thus all your rasters are supposed to overlap perfectly between each other

# let's refresh about our rasters

clc.Loire2 #resolution 250m
occ.ras # resolution 10 000m
ForestCover_prop # resolution 2000m
distah.vrai #resolution 2000m
extent(ForestCover_prop) == extent(distah.vrai)
extent(distah.vrai) == extent(occ.ras)
extent(distah.vrai) == extent(clc.Loire2)

plot(clc.Loire2)
plot(ForestCover_prop)
plot(distah.vrai)
plot(occ.ras)

# as the occurence raster occ.ras would be our rasponse variable
# we will put every raster at this resolution and extent

# here we lose a lot of information

# here we use resample but see before the function aggregate and the existence of disaggregate also
clc.Loire2.ras <- resample(clc.Loire2,occ.ras,method="ngb")
plot(clc.Loire2.ras)
unique(getValues(clc.Loire2.ras))
extent(clc.Loire2.ras) == extent(occ.ras)
ncell(clc.Loire2.ras)
ncell(occ.ras)

ForestCover_prop.ras1 <- aggregate(ForestCover_prop,fact=5,fun="mean")
plot(ForestCover_prop.ras1)
unique(getValues(ForestCover_prop.ras1))
summary(unique(getValues(ForestCover_prop.ras1)))
extent(ForestCover_prop.ras1) == extent(occ.ras)
ncell(ForestCover_prop.ras1)
ncell(occ.ras)
ForestCover_prop.ras2 <- resample(ForestCover_prop,occ.ras,method="bilinear")
plot(ForestCover_prop.ras2)
summary(unique(getValues(ForestCover_prop.ras2)))
extent(ForestCover_prop.ras2) == extent(occ.ras)
ncell(ForestCover_prop.ras2)
ncell(occ.ras)
# we can see that the information is practically the same
# however since we want to overlap perfectly
# we here have to use the function resample

distah.vrai
distah.vrai.ras <- resample(distah.vrai,occ.ras,method="bilinear")
plot(distah.vrai.ras)
summary(unique(getValues(distah.vrai.ras)))
extent(distah.vrai.ras) == extent(occ.ras)
ncell(distah.vrai.ras)
ncell(occ.ras)

distah.vrai.ras <- mask(distah.vrai.ras,occ.ras)
ForestCover_prop.ras2 <- mask(ForestCover_prop.ras2,occ.ras)
clc.Loire2.ras <- mask(clc.Loire2.ras,occ.ras)

# let's create a data frame for analysis
Occurence <- getValues(occ.ras) 
Dist.Foret <- getValues(distah.vrai.ras)
Prop.Foret <- getValues(ForestCover_prop.ras2)
Landuse <- getValues(clc.Loire2.ras)
montableau <- cbind(Occurence,Dist.Foret,Prop.Foret,Landuse)
head(montableau)
nrow(montableau)
class(montableau)
montableau <- data.frame(montableau)
head(montableau)
nrow(montableau)
is.na(montableau)
montableauVrai <- montableau[which(is.na(montableau$Occurence)==F),]
head(montableauVrai)
nrow(montableauVrai)



                          #### %%%%%%%%%%%%%%  ####


    #### G - LANDSCAPE METRICS: BETWEEN RASTER AND VECTOR LAYERS ####


                          #### %%%%%%%%%%%%%%  ####


# ----------------------------------------------------------------------------
    # 1° distance from points to forest ####

# let's refresh and take our points
pointloire
plot(pointloire)
# and our binary raster of forest cover
d.foret
plot(d.foret)
plot(pointloire,add=T,pch=20)

# distance from each points to nearest forest patch for each cell
distances <- distanceFromPoints(d.foret,pointloire)
plot(distances)
distances <- mask(distances,bufloire)
plot(distances)
extent(distances) == extent(d.foret)

# please see: https://rpubs.com/ricardo_ochoa/415839

    

# ----------------------------------------------------------------------------
    # 2° Landscape metrics (easy) ####


library(SDMTools)

# let's refresh and take our landscape raster
clc.Loirebuf
plot(clc.Loirebuf)

# compute the metric percentage of landscape for each class of the raster
?ClassStat
class.clc <- ClassStat(clc.Loirebuf,cellsize=250,bkgd=NA,latlon=F)
names(class.clc)
class.clc$class
length(class.clc$class)
vecnoms <- c("Urban conti","Urban disc","Industrial","Roads","Port","Airports",
             "Extraction sites","Dump","Construction","Green urban","Sport",
             "Arable land","Vineyards","Fruit trees","Pastures","Annual crops",
             "Complex cultivation","agriculture","Broad-leaved forest",
             "Coniferous forest","Mixed forest","Grasslands","Moors",
             "Woodland","Beaches","Inland","Peat bogs","Salt marshes",
             "Salines","Intertidal flats","Water courses","Water bodies",
             "Estuaries","Sea and ocean")
class.clc$prop.landscape
class.clc$landuse <- vecnoms[class.clc$class]
head(class.clc)
head(class.clc[,c("class","landuse","patch.density")])

# compute the patch metric for forest now
# let's take our forest binary raster
d.foret
  # binary matrix
plot(d.foret)
matra <- ConnCompLabel(d.foret)
  # patch value matrix
plot(matra)
patchforest <- PatchStat(matra,cellsize=res(matra)[1],latlon = F)
  # patch metrics
names(patchforest)
# how many patches there is for forest ?
patchforest$patchID
  # one row == one patch

# from this matrix/raster you can compute boundaries of forest patches and then compute indices and distances such as distance to forest edges. See function boundaries in the package raster

# let's rasterize the results of the patch metrics indexes
shapeind <- matra
w <- which(shapeind[]==0)
shapeind[w] <- NA

for(p in 2:(dim(patchforest)[1])) {
  w <- which(shapeind[]==patchforest$patchID[p])
  shapeind[w] <- patchforest$shape.index[p]
}

plot(shapeind) # raster of shape.index metrics for forest




# let's do that for point, buffer and/or polygon

plot(d.foret)
plot(pointloire[6,],add=T,pch=20)
bufr <- buffer(pointloire[6,],20000)
plot(d.foret)
plot(bufr,add=T)
plot(pointloire[6,],add=T)

ra <- mask(crop(d.foret,bufr),bufr)
plot(ra)
plot(bufr,add=T)
plot(pointloire[6,],add=T)
classbuf <- ClassStat(ra,cellsize=res(ra),bkgd=NA,latlon=F)
classbuf
  # we have a binary matrix 0/1 non forest/forest
  # thus it gives a value for each type (2 rows)
  # for this buffer we have a patch density index for forest of 1.025237e-08 



# ----------------------------------------------------------------------------
    # 3° Landscape metrics (not easy) for bufferS, pointS ####

# extraction and analysis of a buffers serie distributed on a regular grid
plot(bufloire)
sample <- spsample(bufloire,n=100,type="regular")
  # create a sample of point regulary spaced inside the shapefile
s <- as.data.frame(sample)
plot(clc.Loirebuf)
#plot(sample,add=T)
points(s[,1],s[,2],pch=3)

# extract the n buffers and create a regular grid
n <- nrow(s)
maille <- s[2,1] - s[1,1]
li <- list("vector", length=n)
for(i in 1:n){
  xmin <- s[i,1] - maille/2
  xmax <- s[i,1] + maille/2
  ymin <- s[i,2] - maille/2
  ymax <- s[i,2] + maille/2
  e <- extent(c(xmin, xmax, ymin, ymax))
  li[[i]] <- crop(clc.Loirebuf,e)
  plot(e,add=T)
}

# we could do the next step in a loop but a lapply is much more elegant ;)
# it means that for each point we compute the classstat metrics
li2 <- lapply(X=li, FUN=ClassStat, cellsize = res(clc.Loirebuf))
length(li2)
li2[[1]]
head(li2[[1]])
names(li2[[1]]) # this is just as before when we did it just for one point
                # this list gives for each point the class stat
                # thus n classstat in a list of length n
li2[[1]][,c("class","n.patches","patch.density")]  
                # if we want to access some particular columns

vecnoms <- c("Urban conti","Urban disc","Industrial","Roads","Port","Airports",
             "Extraction sites","Dump","Construction","Green urban","Sport",
             "Arable land","Vineyards","Fruit trees","Pastures","Annual crops",
             "Complex cultivation","agriculture","Broad-leaved forest",
             "Coniferous forest","Mixed forest","Grasslands","Moors",
             "Woodland","Beaches","Inland","Peat bogs","Salt marshes",
             "Salines","Intertidal flats","Water courses","Water bodies",
             "Estuaries","Sea and ocean")
# this is the vecnoms we used before
# but in our loop from 1 to 44 there will be a row
# thus we have to fill vecnoms as well 
vecnoms2 <- c("Urban conti","Urban disc","Industrial","Roads","Port","Airports",
             "Extraction sites","Dump","Construction","Green urban","Sport",
             "Arable land","","","Vineyards","Fruit trees","","Pastures",
             "Annual crops",
             "Complex cultivation","agriculture","","Broad-leaved forest",
             "Coniferous forest","Mixed forest","Grasslands","Moors","",
             "Woodland","Beaches","","","","","Inland","Peat bogs","Salt marshes",
             "Salines","Intertidal flats","Water courses","Water bodies","",
             "Estuaries","Sea and ocean")


li3 <- list("vector", length=n)
  # create a list of the same length of li2 our number of regular points
nbrow <- length(vecnoms2) # number of landuse
nosnoms <- names(li2[[1]])

for(i in 1:n){
  toto <- data.frame(matrix(NA,ncol= ncol(li2[[i]]),nrow=nbrow))
    # we create a matrix full of NAs
  names(toto) <- nosnoms
    # we give it the names of the ClassStat metrics
  toto[li2[[i]]$class,] <- li2[[i]][,]
    # for the ith point, we gave to toto when the class of the ith point is present
    # the metrics calculated
    # so at the end there is a NA when at this point i the class is not present
    # and there are values when this class is present at the point i
  toto$landuse <- vecnoms2  
    # we inform the landuse
  li3[[i]] <- toto
    # we create a newlist containing the classstat for each point and the landuse
}
length(li3)
  # our number of points
li3[[1]]
# our first point
li3[[1]][,c("class","landuse","patch.density","edge.density")]


#sum(li3[[1]][c(23,24,25),"patch.density"])


# extraction of the pathc density metrics for the type broad-leaved forest ----
pland<-function(x) return(x$patch.density[23])
                      #[23] pour broad-leaved forest
pland_foret <- lapply(X=li3, FUN=pland)
pland_foret
  # here we have a list of 1 == a point
  # containing for each the patch.density index for the borad-leaved forest
pland_foret <- unlist(pland_foret)
pland_foret
            # we have one point not containing any broad-leaved forest

w <- which(is.na(pland_foret)==TRUE)
w

pland_foret[w] <- 0 # replace the NA


# transform results into a raster ----
?rasterFromXYZ
pland_foret.r <- rasterFromXYZ(xyz=data.frame(s, # inform our coordinates
                                              pland_foret), # and what value goes with
                                 res=c(NA,NA), 
                               crs=NA,
                               digits=5)
plot(pland_foret.r, axes=F, box=F)
plot(bufloire,add=T)

projection(pland_foret.r) <- CRS(proj4string(bufloire))


# if we want to include this variable into our analysis data frame

Occurence <- getValues(occ.ras) 
Dist.Foret <- getValues(distah.vrai.ras)
Prop.Foret <- getValues(ForestCover_prop.ras2)
Landuse <- getValues(clc.Loire2.ras)


pland_foret.r
pland_foret.r.ras <- resample(pland_foret.r,occ.ras,method="bilinear")
plot(pland_foret.r.ras,axes=F, box=F)
summary(unique(getValues(pland_foret.r.ras)))
extent(pland_foret.r.ras) == extent(occ.ras)
ncell(pland_foret.r.ras)
ncell(occ.ras)

Patch.Dens <- getValues(pland_foret.r.ras)

montableau <- cbind(Occurence,Dist.Foret,Prop.Foret,Landuse,Patch.Dens)
head(montableau)
nrow(montableau)
montableau <- data.frame(montableau)

montableauVrai2 <- montableau[which(is.na(montableau$Occurence) == F),]



#○ please see : http://j.p.rossi.free.fr/rpackages/ecpaysage/TD_Analyse_quanti_paysage.html#(1)

# ----------------------------------------------------------------------------
    # 4° Landscape metrics for polygons ####

# from Zheng's script ####
#extract land cover data for each poly, given buffer area
bufloire
plot(bufloire)
landpercent <- raster::extract(d.foret, bufloire)
str(landpercent)
class(landpercent)
length(bufloire)
length(landpercent)
# for the 5 polygons 

plot(d.foret)
plot(bufloire,add=T,lwd=2)
# you have a list of values from the raster
unique(getValues(clc.Loire2))
# that are inside your buffer

# summarize each site's data by proportion of each cover type
landpercent2 <- lapply(landpercent, function(x){
  prop.table(table(x))
})

landpercent2[[1]] # for value 19 forest

# extraction de la métrique pour la forêt

pland2 <- function(x) return(x[2]) #[2] pour forest (value 1)
pland_foret <- lapply(landpercent2,FUN=pland2)
bufloire$Prop.Foret <- pland_foret
head(bufloire@data)


fptot <- fortify(bufloire) # can take a little time
data.tot <- bufloire@data
data.tot$id <- rownames(data.tot)
fptot <- merge(fptot,data.tot,by="id") # can take a little time

head(fptot)
fptot$Prop.Foret <- as.numeric(fptot$Prop.Foret)

plotgd <- ggplot() + 
  geom_polygon(data=fptot,aes(x=long,y=lat,group=group,fill=Prop.Foret)) + 
  coord_equal() + 
  coord_quickmap() +
  theme_bw() 
plotgd

plot(ForestCover_prop.ras2)  




                            #### %%%%%%%%%%%%%%  ####


                        #### H - SOME LAST REMARKS ####


                            #### %%%%%%%%%%%%%%  ####

# another package allows for landscape metrics -----
# it is 
# Get the stable version from CRAN
# install.packages("landscapemetrics")
# https://r-spatialecology.github.io/landscapemetrics/

# also spatialEco with the function land.metrics ----
# https://www.rdocumentation.org/packages/spatialEco/versions/1.1-1/topics/land.metrics

# A very interesting landscape metrics is the heterogeneity (and diversity) which is commonly calculated with the Shannon's index. But you might want to take a look
# at the Rao's index ----  
# that also take into account richness and diversity 
# spectralrao function please see the article Rocchini et al.2017
# https://www.sciencedirect.com/science/article/pii/S1470160X16304319
# and the function : 
# https://ars.els-cdn.com/content/image/1-s2.0-S1470160X16304319-mmc1.txt
# github repository:
# https://github.com/mattmar/spectralrao

# a possibility we did not show here is to work with a stack of rasters ----
# it allows to perform the same function only once on all your rasters contained in the stack

# for MNT computations please see the function raster::terrain()----
# elevation <- terrain(mnt.raster,opt="slope",unit="degrees",neighbors = 8)
# neighbors = 8 better for rough surfaces





                            #### %%%%%%%%%%%%%%  ####


                  #### I - REFERENCES AND ACKNOWLEDGMENTS ####


                            #### %%%%%%%%%%%%%%  ####

# Please acknowledge that this script has been constructed with a lot of outside help from a lot of great internet websites (the major ones are informed in the script) and from the help of a lot of people around or not (cesco and co)

# You would be interested to take a look at those sites also:
# https://r-spatial.github.io/sf/articles/
# https://github.com/oliviergimenez/intro_tidyverse/blob/master/tidyogv2.pdf
# https://oliviergimenez.github.io/introspatialR/#16
# https://geocompr.robinlovelace.net/attr.html
# https://statnmap.com/fr/2018-07-14-initiation-a-la-cartographie-avec-sf-et-compagnie/
# https://nceas.github.io/oss-lessons/spatial-data-gis-law/3-mon-intro-gis-in-r.html

# moreover, a lot of other packages are available such as cartography, maptools and others and others in preparation ;) (like stars a future sf package for raster)
# possibility of creating interactive maps as well
# and to open GRASS in R (https://cran.r-project.org/web/packages/rgrass7/index.html)

