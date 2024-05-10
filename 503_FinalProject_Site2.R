# Final Project - Seepage velocity maps for Alder Creek Site 2
# Hanna Szydlowski
# GEOL 503 Spring 2024
# 02/28/24

## Objective: Create 3 seepage velocity maps for Site 2 at Alder Creek (located in
#  Ontario, Canada). Maps will be: mini piezometer seepage velocities calculated
#  using K values from past literature (map 4), mini piezometer seepage
#  velocities using K values from field temperature readings (map 5), and SBPVP
#  seepage velocity (map 6).

# Map 4 - Mini Piezometer seepage velocity using K from literature
# Map 5 - Mini Piezometer seepage velocity using K from UW temperature data
# Map 6 - SBPVP seepage velocity

## Lines of code that need to be commented in/out depending on map:
#  lines 33-38 (data files)
#  lines 66-67 (K values) 
#  lines 126-127 (for map 3 only - location coordinates to make mini piez)
#  lines 177-178 (for map 3 only - min and max values for x and y axes)

#----------
## Set working directory
setwd("C:/Users/hszyd/Documents/KU Research/UW_Work/Alder_Creek/503FinalProject")

#----------
## List of .csv files to be read in
#  Note: Not all files will be used at once. The code is written using generic
#  names for the dataframes, so some files will need to be commented out depending
#  on the map being made.
#  Determine which files to use based on the comment to the right

SVData <- read.csv(file = "MiniPiezData_Site2.csv")               # map 4 & 5
KvaluesUW <- read.csv(file = "MiniPiez_KvaluesUW_Site2.csv")      # map 5
#SVData <- read.csv(file = "SBPVP_SVmday_Site2.csv")               # map 6
#colnames(SVData)[2] <- c("AVG_SV_mday")                           # map 6
LocationCoord <- read.csv(file = "Jeff_Site2_MiniPiez_Coord.csv") # map 4 & 5
#LocationCoord <- read.csv(file = "Jeff_Site2_SBPVP_Coord.csv")    # map 6

#----------
# load in packages
library(sp)
library(gstat)
library(ggplot2)
library(lattice)
library(reshape)
library(latticeExtra)

#----------
# Step 1: Calculate the seepage velocities for the mini piezometer locations
# Note: Step 1 only needs to be done for maps 4 and 5
# Note again: Please make sure to use the correct K value as labeled 

# Step 2: Convert coordinates to UTM and create a dataframe with 3 columns:
# Easting_m, Northing_m, and SV_mday

# Step 3: Interpolate the data to get seepage velocities over the entire stream
# reach

#-------------------------------------------------------------------------------

## Step 1: Calculate the seepage velocities for the mini piezometer locations
#  Note: Step 1 only needs to be done for maps 4 and 5

# define hydraulic conductivity (K) and porosity (n)
#K <- 3.5e-2  # cm/sec from literature
K <- KvaluesUW$K_cmsec # cm/sec from UW temperature data
n <- 0.3

# Convert na characters to NA numeric
ncol(SVData)
SVData[2:4] <- as.numeric(unlist(SVData[2:4]))

# Create new columns that calculate the hydraulic gradient for each day
# hydraulic gradient (HG) = WL above stream / depth
SVData$HG_20 <- (SVData$WLaS_20)/(SVData$Depth_20)

# Create new columns that calculate the seepage velocity for each day (cm/sec)
# seepage velocity (SV) (cm/sec) = (K/n)*HG
SVData$SV_20_cmsec <- (K/n)*SVData$HG_20

# Create new columns that calculate the seepage velocity for each day (cm/day)
# seepage velocity (SV) (cm/day) = SV (cm/sec) * # of seconds in a day
# define how many seconds are in a day
secondsInDay <- 86400  # sec
SVData$SV_20_cmday <- SVData$SV_20_cmsec*secondsInDay

# Create new columns that calculate the seepage velocity for each day (m/day)
# seepage velocity (SV) (m/day) = SV (cm/day) / # of cm in a m (100cm)
SVData$AVG_SV_mday <- SVData$SV_20_cmday/100

#-------------------------------------------------------------------------------

## Step 2: Convert coordinates to UTM and create a dataframe with 3 columns:
#          Easting_m, Northing_m, and SV_mday

# plot the x/y points to view orientation
# first using Lat/Long
plot(LocationCoord$Long..W., LocationCoord$Lat..N.,
     xlab = "Longitude", ylab = "Latitude") # north pointing up

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
LocationCoord_UTM_SV_wNA <- as.data.frame(UTMCoords)
colnames(LocationCoord_UTM_SV_wNA) <- c("Easting_m", "Northing_m")
LocationCoord_UTM_SV_wNA$SV_mday <- SVData$AVG_SV_mday
LocationCoord_UTM_SV_wNA

# Change the coordinate system to be in local meter scale coordinates
# do this by subtracting the minimum values from easting and northing

#----FOR MAP 3 ONLY----#
LocationCoord_UTM_SV_wNA$x <- LocationCoord_UTM_SV_wNA$Easting_m - 538811.7
LocationCoord_UTM_SV_wNA$y <- LocationCoord_UTM_SV_wNA$Northing_m - 4797728
#----------------------#

# DO NOT RUN LINES 131-134 FOR MAP 3
LocationCoord_UTM_SV_wNA$x <- LocationCoord_UTM_SV_wNA$Easting_m - 
  min(LocationCoord_UTM_SV_wNA$Easting_m)
LocationCoord_UTM_SV_wNA$y <- LocationCoord_UTM_SV_wNA$Northing_m - 
  min(LocationCoord_UTM_SV_wNA$Northing_m)

# plot the UTM coordinates to view orientation
#plot(LocationCoord_UTM_SV$Northing_m, LocationCoord_UTM_SV$Easting_m,
#     xlab = "Northing (m)", ylab = "Easting (m)") # north pointing to the left
plot(LocationCoord_UTM_SV_wNA$Easting_m, LocationCoord_UTM_SV_wNA$Northing_m,
    xlab = "Easting (m)", ylab = "Northing (m)") # north pointing upwards

# remove all NA values from the seepage velocity (SV_mday) if there are any
LocationCoord_UTM_SV <- LocationCoord_UTM_SV_wNA[- which(is.na(LocationCoord_UTM_SV_wNA$SV_mday)),]

#-------------------------------------------------------------------------------

# Mess around with plotting to see distribution before interpolation if desired

# use ggplot to plot the velocities indicated by size and color
ggplot(
  data = LocationCoord_UTM_SV,
  mapping = aes(x = Easting_m, y = Northing_m)
) +
  geom_point(mapping = aes(size = SV_mday, color = SV_mday)
) +
  xlim(538811.7, 538825.3) +
  ylim(4797728, 4797764) +
  scale_size_continuous(breaks = c(0,5,10,15,20,25,30) # change as needed
) +
  labs(
    title = "FILL TITLE HERE",
    subtitle = "FILL SUBTITLE HERE",
    x = "Easting (m)", y = "Northing (m)",
    size = "Seepage Velocity (m/day)",
    color = "Seepage Velocity (m/day)"
)

#-------------------------------------------------------------------------------

## Step 3: Interpolate the data to get seepage velocities over the entire stream
#          reach

# Create the interpolation grid
# Manually set ranges and make expanded grid based on min and max values

#-----FOR MAP 3 ONLY----#
x.range <- as.integer(c(-1, 15))
y.range <- as.integer(c(-5, 41))
#-----------------------

# DO NOT RUN LINES 182-185 FOR MAP 3
x.range <- as.integer(c((min(LocationCoord_UTM_SV$x))-2,
                        (max(LocationCoord_UTM_SV$x))+2))
y.range <- as.integer(c((min(LocationCoord_UTM_SV$y))-5,
                      (max(LocationCoord_UTM_SV$y))+5))

# Generate the grid
LocationGrid <- expand.grid(x = seq(x.range[1], x.range[2], 0.1), 
                 y = seq(y.range[1], y.range[2], 0.1))
plot(LocationGrid$x, LocationGrid$y, col = "grey") # plot the gird
points(LocationCoord_UTM_SV$x, LocationCoord_UTM_SV$y, pch = 19,
       col = "black") # plot the points on the grid

# Set the idp for interpolation
idp = 2

# Begin the interpolation on the seepage velocities (SV_mday)
# use inverse distance weighting
idwSV <- idw(SV_mday~1, ~x + y , data = LocationCoord_UTM_SV,
               newdata = LocationGrid, idp = idp)
# Change the column name for column 3
colnames(idwSV)[3]<-c("idwSV_mday")

# plot the idw results to show seepage velocity distribution
contourplot(idwSV_mday~x+y,idwSV,
            main = paste("FILL TITLE HERE"),
            cuts=10,                  #factors defining panes
            region=T, 
            #aspect="iso",
            filled=T,
            contour = F, 
            col.regions=heat.colors(100,alpha=.75,rev=T))+
  as.layer(xyplot(y~x, LocationCoord_UTM_SV_wNA,
                  cex=2,pch=20,col="black"),
           panel= panel.polygon(x=c(-2,-2,2,5,10,13), # lower
                                y=c(-5,29,19,7,-1,-5),
                                col="darkgreen"),
           panel= panel.polygon(x=c(0,3,5,8,10,16,16), # upper
                                y=c(41,37,31,24,17,8,41),
                                col="darkgreen"))


#----------
# Export idvSV to .csv file (I need this to make a larger figure)

write.csv(idwSV,
          "C:/Users/hszyd/Documents/KU Research/UW_Work/Alder_Creek/503FinalProject\\Site2_SBPVP_idwTable.csv" )



#------this code did not work but maybe someday I'll come back to it ¯\_(ツ)_/¯
# It was the original attempt to make the contour plot as done in the above code

# Reference data to grid
xidw <- idwSV$x # x vector on expanded grid
yidw <- idwSV$y # y vector on expanded grid
# vector of x values for contour plot
xcontour <- seq(from = as.integer((min(LocationCoord_UTM_SV$x))-5), 
                to = as.integer((max(LocationCoord_UTM_SV$x))+5), by = 1)
# vector of y values for contour plot
ycontour <- seq(from = as.integer((min(LocationCoord_UTM_SV$y))-5),
                to = as.integer((max(LocationCoord_UTM_SV$y))+5), by = 1)

# Plot the empty scaled area
plot(xidw, yidw, pch=1, col="transparent", ylim = range(yidw),
     xlim = range(xidw), xlab = "Easting (m)", ylab = "Northing (m)")
# Plot the points over the empty area
points(LocationCoord_UTM_SV$x, LocationCoord_UTM_SV$y, pch = 19, col = "black")

par(new = T)

xlength <- as.numeric(length(xcontour))
ylength <- as.numeric(length(ycontour))

# Reshape as wide view
mtrxSV <- matrix(rep(rep(NA, xlength), ylength), ncol = ylength)

jj=0
for(i in 1:xlength){ # inserting interpolation results into expanded grid in wide-view
  for(j in 1:ylength){
    jj <- jj + 1
    mtrxSV[i,j] <- idwSV[jj, 3]
  }
}

# Transpose
# x-direction is now down columns
# y-direction is now along rows
mtrxSV_t <- t(mtrxSV)

contour(xcontour, ycontour, mtrxSV, nlevels = 3, xlim = range(xcontour),
        ylim = range(ycontour))

# Use lattice to fill in plot with color using lattice::contourplot
# contourplot requires data to be entered in long-view

# first step is to convert mtrxSV to long view
# this means entering x-coordinates as row names and
# y-coordinates as column names
mtrxSV_2.0 <- mtrxSV #create new dataframe for lattice work
colnames(mtrxSV_2.0) <- ycontour #moving along row is y direction
rownames(mtrxSV_2.0) <- xcontour #moving down columns (row by row) is x direction

# Melt to long-view form; this creates columns for x, y, SV_mday
mtrxSV_2.0 <- melt(mtrxSV_2.0, na.rm = FALSE, value.name = "SV_mday",
                   id = c("x","y","SV_mday"))
colnames(mtrxSV_2.0)<-c("x","y","SV_mday") #name columns

# Make contourplot with aspect preserved and selected colors
# alpha gives some transparency to colors, lowering intensity
# rev in colors makes red the 'hot' zone (default puts yellow there)
contourplot(mtrxSV_2.0$SV_mday~mtrxSV_2.0$x+mtrxSV_2.0$y,
            main = paste("FILL TITLE HERE"),
            cuts = 21, 
            region = T, 
            aspect = "iso",
            xlim = range(xcontour),
            xlab = "Northing (m)",
            ylim = range(ycontour),
            ylab = "Easting (m)",
            contour = F,
            col.regions = heat.colors(20, alpha = 1, rev = T)) +
  as.layer(xyplot(LocationCoord_UTM_SV$y ~ LocationCoord_UTM_SV$x, pch = 19,
                  col = "black"))
