# Final Project - Seepage velocity Darcy calculations for Alder Creek
# Hanna Szydlowski
# GEOL 503 Spring 2024
# 02/28/24

## Step 1: Calculate the seepage velocities for the mini piezometer locations 
#          based on K values from past studies (for Site 1)
## Step 2: Interpolate the data to get seepage velocities over the entire stream reach
#          based on the past studies K calculation
## Step 3: Repeat the same steps for Site 2 for K values from past studies

## Step 4: Repeat the same steps for Waterloo K values for Site 1 and 2

## Step 5: Repeat the same step for SBPVP values for Site 1 and 2

#-------------------------------------------------------------------------------

## Step 1: Calculate the seepage velocities for the mini piezometer locations 
#          based on K values from past studies

# Site 1 - James

# Set working directory
setwd("C:/Users/hszyd/Documents/KU Research/UW_Work/Alder_Creek/Corey_data")

# Read in .csv file with field data
MiniPiezData_S1 <- read.csv(file = "MiniPiezData.csv")
MiniPiezData_S1

# Key for what each header represents --->
# Name = ID number of the mini piezometer
# WLS_ = Water Level of Stream from streambed to water level (cm)
# WLaS_ = Water Level Above Stream from stream level to level in mini piez (cm)
# Depth_ = depth of mini piez below the streambed (cm)
# 21, 23, 27 = day of the month (October) data is from

# Convert na characters to NA numeric
ncol(MiniPiezData_S1)
MiniPiezData_S1[2:10] <- as.numeric(unlist(MiniPiezData_S1[2:10]))

# Create new columns that calculate the hydraulic gradient for each day
# hydraulic gradient (HG) = WL above stream / depth
MiniPiezData_S1$HG_21 <- (MiniPiezData_S1$WLaS_21)/(MiniPiezData_S1$Depth_21)
MiniPiezData_S1$HG_23 <- (MiniPiezData_S1$WLaS_23)/(MiniPiezData_S1$Depth_23)
MiniPiezData_S1$HG_27 <- (MiniPiezData_S1$WLaS_27)/(MiniPiezData_S1$Depth_27)

# define hydraulic conductivity (K) and porosity (n)
K <- 7.4e-2  # cm/sec
n <- 0.3

# Create new columns that calculate the seepage velocity for each day (cm/sec)
# seepage velocity (SV) (cm/sec) = (K/n)*HG
MiniPiezData_S1$SV_21_cmsec <- (K/n)*MiniPiezData_S1$HG_21
MiniPiezData_S1$SV_23_cmsec <- (K/n)*MiniPiezData_S1$HG_23
MiniPiezData_S1$SV_27_cmsec <- (K/n)*MiniPiezData_S1$HG_27

# Create new columns that calculate the seepage velocity for each day (cm/day)
# seepage velocity (SV) (cm/day) = SV (cm/sec) * # of seconds in a day
# define how many seconds are in a day
secondsInDay <- 86400  # sec
MiniPiezData_S1$SV_21_cmday <- MiniPiezData_S1$SV_21_cmsec*secondsInDay
MiniPiezData_S1$SV_23_cmday <- MiniPiezData_S1$SV_23_cmsec*secondsInDay
MiniPiezData_S1$SV_27_cmday <- MiniPiezData_S1$SV_27_cmsec*secondsInDay

# Create new columns that calculate the seepage velocity for each day (m/day)
# seepage velocity (SV) (m/day) = SV (cm/day) / # of cm in a m (100cm)
MiniPiezData_S1$SV_21_mday <- MiniPiezData_S1$SV_21_cmday/100
MiniPiezData_S1$SV_23_mday <- MiniPiezData_S1$SV_23_cmday/100
MiniPiezData_S1$SV_27_mday <- MiniPiezData_S1$SV_27_cmday/100

# Create new columns that calculate the average seepage velocity over the three days
# do this for cm/sec, cm/day, and m/day

# first figure out what the number of each column is
match("SV_21_cmsec", names(MiniPiezData_S1))
match("SV_21_cmday", names(MiniPiezData_S1))
match("SV_21_mday", names(MiniPiezData_S1))

MiniPiezData_S1$AVG_SV_cmsec <- rowMeans(MiniPiezData_S1[, c(14,15,16)], na.rm = TRUE)
MiniPiezData_S1$AVG_SV_cmday <- rowMeans(MiniPiezData_S1[, c(17,18,19)], na.rm = TRUE)
MiniPiezData_S1$AVG_SV_mday <- rowMeans(MiniPiezData_S1[, c(20,21,22)], na.rm = TRUE)

#-------------------------------------------------------------------------------


## Step 2: Interpolate the data to get seepage velocities over the entire stream
#          reach based on the past studies K calculation

# install and load in packages
#install.packages("sp")
library(sp)
#install.packages("gstat")
library(gstat)
library(ggplot2)
#install.packages("akima")
library(akima)
library(lattice)

# Reset working directory for coordinates data
setwd("C:/Users/hszyd/Documents/KU Research/UW_Work/Alder_Creek/Coordinates")

# Read in .csv file with coordinate data
MPCoord <- read.csv(file = "James_Site1_MiniPiez_Coord.csv")
MPCoord

# plot the x/y points to view orientation
# first using Lat/Long
plot(MPCoord$Lat..N., MPCoord$Long..W.) # north pointing left
plot(MPCoord$Long..W., MPCoord$Lat..N.) # north pointing up

# Give a spatial extent going from lat/long to UTM
# First - set the coordinates (column names in MPCoord)
coordinates(MPCoord) <- ~ Long..W. + Lat..N.
# Second - view what the class of data is
class(MPCoord) # using the sp package
str(MPCoord) # view the string of data
# Third - give the data a CRS
proj4string(MPCoord) <- CRS("+proj=longlat +datum=WGS84")
# Fourth - transform the data to the CRS you want
MPCoord_UTM <- spTransform(MPCoord, "+proj=utm +zone=17 +datum=WGS84")
MPCoord_UTM

# convert MPCoord back to a dataframe
# MPCoord <- as.data.frame(MPCoord)

# make a new dataframe with the UTM coordinates and seepage velocity
# this will be the dataframe to use for the interpolation
UTMCoords <- MPCoord_UTM@coords
MPCoord_UTM_SV <- as.data.frame(UTMCoords)
colnames(MPCoord_UTM_SV) <- c("Easting_m", "Northing_m")
MPCoord_UTM_SV$AVG_SV_mday <- MiniPiezData_S1$AVG_SV_mday
MPCoord_UTM_SV

# plot the UTM coordinates to view orientation
plot(MPCoord_UTM_SV$Northing_m, MPCoord_UTM_SV$Easting_m) # north pointing to the left
plot(MPCoord_UTM_SV$Easting_m, MPCoord_UTM_SV$Northing_m) # north pointing upwards

# remove all NA values
MPCoord_UTM_SV <- MPCoord_UTM_SV[- which(is.na(MPCoord_UTM_SV$AVG_SV_mday)),]

#------------------

#mess around with plotting to see distribution before interpolation

# use ggplot to plot the velocities indicated by size
ggplot(
  data = MPCoord_UTM_SV,
  mapping = aes(x = Northing_m, y = Easting_m)
) +
  geom_point(mapping = aes(size = AVG_SV_mday)
) +
  labs(
    title = "Seepage Velocity Distribution - Site 1",
    subtitle = "Using K values from the literature",
    x = "Northing (m)", y = "Easting (m)",
    size = "Seepage Velocity (m/day)")

# use ggplot to plot the velocities indicated by color
ggplot(
  data = MPCoord_UTM_SV,
  mapping = aes(x = Northing_m, y = Easting_m)
) +
  geom_point(mapping = aes(color = AVG_SV_mday)
) +
  labs(
    title = "Seepage Velocity Distribution - Site 1",
    subtitle = "Using K values from the literature",
    x = "Northing (m)", y = "Easting (m)",
    color = "Seepage Velocity (m/day)")






gridV<- akima::interp(MPCoord_UTM_SV$Northing_m,
                      MPCoord_UTM_SV$Easting_m,
                      MPCoord_UTM_SV$AVG_SV_mday)

gridVdf <- data.frame(x = rep(gridV$x), ncol(gridV$z),
                      y = rep(gridV$y, each = nrow(gridV$z)),
                      z = as.numeric(gridV$z))

ggplot(data = gridVdf, aes(x = x, y = y, z = z) +
         geom_contour_filled(size = 0.3, bins = 20) )











#-------------------------------------------------------------------------------
### FAILED CODE BELOW vvvvv




#--------------------
## There are 3 steps to implement kriging (https://rpubs.com/nabilabd/118172)
  # 1. Convert the dataframe to a spatial points dataframe
  # 2. Fit a variogram model to the data
  # 3. Krige the data according to the variogram
#--------------------

## 1. Convert the dataframe to a spatial points dataframe

# all values need to be numerical values not NA values
# drop row 13 which has an NA
MiniPiezCoord <- MiniPiezCoord[- which(is.na(MiniPiezCoord$AVG_SV_mday)),]

# view the MiniPiezCoord dataframe
str(MiniPiezCoord)
# 'data.frame':	24 obs. of  6 variables:
#   $ Mini.Piez  : chr  "1" "1b" "2" "3" ...
# $ Lat..N.    : num  43.3 43.3 43.3 43.3 43.3 ...
# $ Long..W.   : num  -80.5 -80.5 -80.5 -80.5 -80.5 ...
# $ AVG_SV_mday: num  41 64.8 51.5 63.9 50.3 ...

# Specify which columns are the coordinates
# long first and then lat
coordinates(MiniPiezCoord) <- ~ Long..W. + Lat..N.
class(MiniPiezCoord) # using the sp package
str(MiniPiezCoord)

# access the corner square values
bbox(MiniPiezCoord)
#                min       max
# Long..W. -80.52694 -80.52653
# Lat..N.   43.34653  43.34730

# Check the spatial extent of the data
# not specified yet so the output will be NA
proj4string(MiniPiezCoord)

# Give a spatial extent going from lat/long to UTM
proj4string(MiniPiezCoord) <- CRS("+proj=longlat +datum=WGS84")
MiniPiezCoord_UTM <- spTransform(MiniPiezCoord, "+proj=utm +zone=17 +datum=WGS84")
MiniPiezCoord_UTM
# make a new dataframe with the UTM coordinates
UTMCoords <- MiniPiezCoord_UTM@coords
colnames(UTMCoords) <- c("Easting (m)", "Northing (m)")
# plot the UTM coordinates
plot(y = UTMCoords[,1], x = UTMCoords[,2]) # north pointing to the left
plot(x = UTMCoords[,1], y = UTMCoords[,2]) # north pointing upwards


#------------------
## 2. Fit a variogram model to the data

# AVG_SV_mday is the variable being used to find missing seepage velocities
MPCoord_vgm <- variogram(MiniPiezCoord$AVG_SV_mday~1, MiniPiezCoord) # calculates sample variogram values

MPCoord_fit <- fit.variogram(MPCoord_vgm, model=vgm(1, "Sph", 200, 1)) # fit model

# plot to see how well the fit is
plot(MPCoord_vgm, MPCoord_fit) # plot the sample values, along with the fit model




