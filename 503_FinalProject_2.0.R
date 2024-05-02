# Final Project - Seepage velocity maps for Alder Creek
# Hanna Szydlowski
# GEOL 503 Spring 2024
# 02/28/24

## Objective: Create 6 seepage velocity maps for Alder Creek (located in
#  Ontario, Canada). There are two sites of interest that will be referred to as
#  Site 1 and Site 2. The 6 maps will be split between 3 for Site 1 and 3 for 
#  Site 2. Maps will be: mini piezometer seepage velocities calculated using K
#  values from past literature (map 1 & 4), mini piezometer seepage velocities
#  using K values from field temperature readings (map 2 & 5), and SBPVP seepage
#  velocity (map 3 & 6).

#----------
## Set working directory
setwd("C:/Users/hszyd/Documents/KU Research/UW_Work/Alder_Creek/503FinalProject")

#----------
## List of .csv files to be read in
#  Note: Not all files will be used at once. The code is written using generic
#  names for the dataframes, so some files will need to be commented out depending
#  on the map being made.

MiniPiezData <- read.csv(file = "MiniPiezData_Site1.csv") # Site 1 (map 1)
LocationCoord <- read.csv(file = "James_Site1_MiniPiez_Coord.csv") # Site 1 (map 1 & 2)
#LocationCoord <- read.csv(file = "James_Site1_SBPVP_Coord.csv") # Site 1 (map 3)
#MiniPiezData <- read.csv(file = "_________") # Site 2 (map 4)
LocationCoord <- read.csv(file = "Jeff_Site2_MiniPiez_Coord.csv") # Site 2 (map 4 & 5)
#LocationCoord <- read.csv(file = "Jeff_Site2_SBPVP_Coord.csv") # Site 2 (map 6)

#----------
# Step 1: Calculate the seepage velocities for the mini piezometer locations
# Note: Step 1 only needs to be done for maps 1 and 4

# Step 2: Convert coordinates to UTM and create a dataframe with 3 columns:
# Easting_m, Northing_m, and SV_mday

# Step 3: Interpolate the data to get seepage velocities over the entire stream
# reach

#-------------------------------------------------------------------------------

## Step 1: Calculate the seepage velocities for the mini piezometer locations
#  Note: Step 1 only needs to be done for maps 1 and 4

# Header key for MiniPiezData_Site1 and ______________ --->
# Name = ID number of the mini piezometer
# WLS_ = Water Level of Stream from streambed to water level (cm)
# WLaS_ = Water Level Above Stream from stream level to level in mini piez (cm)
# Depth_ = depth of mini piez below the streambed (cm)
# 21, 23, 27 = day of the month (October) data is from

# Convert na characters to NA numeric
ncol(MiniPiezData)
MiniPiezData[2:10] <- as.numeric(unlist(MiniPiezData[2:10]))

# Create new columns that calculate the hydraulic gradient for each day
# hydraulic gradient (HG) = WL above stream / depth
MiniPiezData$HG_21 <- (MiniPiezData$WLaS_21)/(MiniPiezData$Depth_21)
MiniPiezData$HG_23 <- (MiniPiezData$WLaS_23)/(MiniPiezData$Depth_23)
MiniPiezData$HG_27 <- (MiniPiezData$WLaS_27)/(MiniPiezData$Depth_27)

# define hydraulic conductivity (K) and porosity (n)
K <- 7.4e-2  # cm/sec
n <- 0.3

# Create new columns that calculate the seepage velocity for each day (cm/sec)
# seepage velocity (SV) (cm/sec) = (K/n)*HG
MiniPiezData$SV_21_cmsec <- (K/n)*MiniPiezData$HG_21
MiniPiezData$SV_23_cmsec <- (K/n)*MiniPiezData$HG_23
MiniPiezData$SV_27_cmsec <- (K/n)*MiniPiezData$HG_27

# Create new columns that calculate the seepage velocity for each day (cm/day)
# seepage velocity (SV) (cm/day) = SV (cm/sec) * # of seconds in a day
# define how many seconds are in a day
secondsInDay <- 86400  # sec
MiniPiezData$SV_21_cmday <- MiniPiezData$SV_21_cmsec*secondsInDay
MiniPiezData$SV_23_cmday <- MiniPiezData$SV_23_cmsec*secondsInDay
MiniPiezData$SV_27_cmday <- MiniPiezData$SV_27_cmsec*secondsInDay

# Create new columns that calculate the seepage velocity for each day (m/day)
# seepage velocity (SV) (m/day) = SV (cm/day) / # of cm in a m (100cm)
MiniPiezData$SV_21_mday <- MiniPiezData$SV_21_cmday/100
MiniPiezData$SV_23_mday <- MiniPiezData$SV_23_cmday/100
MiniPiezData$SV_27_mday <- MiniPiezData$SV_27_cmday/100

# Create new columns that calculate the average seepage velocity over the three days
# do this for cm/sec, cm/day, and m/day

# first figure out what the number of each column is
match("SV_21_cmsec", names(MiniPiezData))
match("SV_21_cmday", names(MiniPiezData))
match("SV_21_mday", names(MiniPiezData))

MiniPiezData$AVG_SV_cmsec <- rowMeans(MiniPiezData[, c(14,15,16)], na.rm = TRUE)
MiniPiezData$AVG_SV_cmday <- rowMeans(MiniPiezData[, c(17,18,19)], na.rm = TRUE)
MiniPiezData$AVG_SV_mday <- rowMeans(MiniPiezData[, c(20,21,22)], na.rm = TRUE)

#-------------------------------------------------------------------------------

## Step 2: Convert coordinates to UTM and create a dataframe with 3 columns:
#          Easting_m, Northing_m, and SV_mday

# load in packages
library(sp)
library(gstat)
library(ggplot2)
library(lattice)

# plot the x/y points to view orientation
# first using Lat/Long
plot(LocationCoord$Lat..N., LocationCoord$Long..W.) # north pointing left
plot(LocationCoord$Long..W., LocationCoord$Lat..N.) # north pointing up

# Give a spatial extent going from lat/long to UTM
# First - set the coordinates (column names from LocationCoord)
coordinates(LocationCoord) <- ~ Long..W. + Lat..N.
# Second - view what the class of data is
class(LocationCoord) # using the sp package
str(LocationCoord) # view the string of data
# Third - give the data a CRS
proj4string(LocationCoord) <- CRS("+proj=longlat +datum=WGS84")
# Fourth - transform the data to the CRS you want
LocationCoord_UTM <- spTransform(LocationCoord, "+proj=utm +zone=17 +datum=WGS84")
LocationCoord_UTM

# Make a new dataframe with the UTM coordinates and seepage velocity
# This will be the dataframe to use for the interpolation
UTMCoords <- LocationCoord_UTM@coords
LocationCoord_UTM_SV <- as.data.frame(UTMCoords)
colnames(LocationCoord_UTM_SV) <- c("Easting_m", "Northing_m")
LocationCoord_UTM_SV$SV_mday <- MiniPiezData$AVG_SV_mday
LocationCoord_UTM_SV

# plot the UTM coordinates to view orientation
plot(LocationCoord_UTM_SV$Northing_m, LocationCoord_UTM_SV$Easting_m) # north pointing to the left
#plot(LocationCoord_UTM_SV$Easting_m, LocationCoord_UTM_SV$Northing_m) # north pointing upwards

# remove all NA values from the seepage velocity (SV_mday) if there are any
LocationCoord_UTM_SV <- LocationCoord_UTM_SV[- which(is.na(LocationCoord_UTM_SV$AVG_SV_mday)),]

#----------

# Mess around with plotting to see distribution before interpolation if desired

# use ggplot to plot the velocities indicated by size
ggplot(
  data = LocationCoord_UTM_SV,
  mapping = aes(x = Northing_m, y = Easting_m)
) +
  geom_point(mapping = aes(size = SV_mday)
) +
  labs(
    title = "FILL IN TITLE HERE",
    subtitle = "FILL IN SUBTITLE HERE",
    x = "Northing (m)", y = "Easting (m)",
    size = "Seepage Velocity (m/day)")

# use ggplot to plot the velocities indicated by color
ggplot(
  data = LocationCoord_UTM_SV,
  mapping = aes(x = Northing_m, y = Easting_m)
) +
  geom_point(mapping = aes(color = SV_mday)
) +
  labs(
    title = "FILL IN TITLE HERE",
    subtitle = "FILL IN SUBTITLE HERE",
    x = "Northing (m)", y = "Easting (m)",
    color = "Seepage Velocity (m/day)")


#--------------edit starting here

## Step 3: Interpolate the data to get seepage velocities over the entire stream
#          reach

# Create the interpolation grid
# Manually set ranges and make expanded grid based on min and max values
x.range <- as.integer(c(min(LocationCoord_UTM_SV$Northing_m),
                        max(LocationCoord_UTM_SV$Northing_m)))
y.range<-as.integer(c(min(LocationCoord_UTM_SV$Easting_m),
                      max(LocationCoord_UTM_SV$Easting_m)))


LocationGrid<-expand.grid(x=seq(x.range[1], x.range[2], 1), 
                 y=seq(y.range[1], y.range[2], 1))
plot(LocationGrid, col = "grey")
points(LocationCoord_UTM_SV$Northing_m, LocationCoord_UTM_SV$Easting_m,
       col = "black")

idp = 2

colnames(LocationCoord_UTM_SV)[1:2] <- c("x" , "y")

idwlogbtx<-idw(AVG_SV_mday~1, ~x + y , data = LocationCoord_UTM_SV,
               newdata=grd, idp=idp) #interpolate

idwoutput=as.data.frame(idwlogbtx) #output to dataframe
names(idwoutput)[1:3]<-c("x","y","logbtx") #give columns names
# this leaves 4th column with name var1.var - i.e. unnamed

xidw<-idwoutput$x #x vector on expanded grid
yidw<-idwoutput$y #y vector on expanded grid
xcontour<-x.range #vector of x values for contour plot
ycontour<-y.range #vector of y values for contour plot

plot(xidw,yidw,pch=1,col="transparent", #empty,scaled plot space
     ylim=rev(range(yidw)),xlim=range(xidw))
points(LocationCoord_UTM_SV$x,LocationCoord_UTM_SV$y,pch=1, col = "black")




