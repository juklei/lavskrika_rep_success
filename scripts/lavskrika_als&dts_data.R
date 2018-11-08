## This settlript aims at importing relevant raster files derived from ALS data
## create habitat metrics and calculate the distance to the closest settlement
## for Siberian jay territories at different radiuses around the nest.

## First edit: 20171208
## Last edit: 20181107

## Author: Julian Klein

## Clear environment and load libraries ----------------------------------------

rm(list = ls())

library(data.table)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(dplyr)

## Define or source functions used in this settlript ------------------------------

#extract_by_nest <- function(x) {AALS_by_yearALS_by_year[[paste0("y_", x$year)]], x[, c("X", "Y")],
                 buffer = s_rad)
  all <- all[[1]]
  
  # set.seed(3)
  # rrows <- sample(length(all[, 1]), s_size) 
  # all <- all[rrows, ]
  
  nd_ratio <- sum(na.omit(all[, 1] == -9999))/sum(!is.na(all[, 1]))
  
  if(nd_ratio < nd_thresh) {
    
    all[all[] == -9999] <-  NA
    aval_dense <- sum(all[!is.na(all[, 2]), 2] > 8.8)/sum(!is.na(all[, 2]))
    mean <- colMeans(all, na.rm = TRUE)
    sd <- apply(all, 2, FUN = sd, na.rm = TRUE)
    vals <- cbind(rbind(mean, sd, aval_dense), var = c("mean", "sd", "aval_dense"))
    
  } else {
    
    all[] <- NA
    vals <- cbind(unique(all), var = "excluded")
    
  }
  
  return(cbind(x, vals))
  
}

## Load and explore data -------------------------------------------------------

dir("data")

## Text files
nest_pos <- read.csv("data/nest_positions.csv")
settl <- read.csv("data/settlement_positions.csv")

## Raster files
ALS <- merge(stack(c("data/unmanaged_ElevP95.asc",
                     "data/unmanaged_5_Perc_above.asc",
                     "data/unmanaged_0.5_Perc_above.asc")),
              stack(c("data/managed_ElevP95.asc",
                      "data/managed_5_Perc_above.asc",
                      "data/managed_0.5_Perc_above.asc")))
names(ALS) <- c("height", "vd_5to", "vd_0to")

## Shape files
study_area <- shapefile("data/120ha_buffer_study.shp")
forestry <- shapefile("data/forestry.shp")

## -----------------------------------------------------------------------------

## Process ALS data and adjust it to forestry interventions --------------------

## Calculate percentage returns between 0.5 and 5m
ALS$vd_0to5 <- ALS$vd_0to - ALS$vd_5to

## Drop percentage returns above 0.5
ALS <- ALS[[-3]]

## Reduce ALS data set to study area
ALS <- mask(ALS, study_area)

## Add forest area layer: All pixels with forest height below 2 m or NA(water)
ALS$area <- !(ALS[["height"]][] < 2 | is.na(ALS[["height"]][]))

## Look at ALS metrics correlations
layerStats(ALS, stat = "pearson", na.rm = TRUE)

## Create a ALS data sample from all clear cuts in the study area. We use this 
## sample from later on in the settlript.
sample_cc <- mask(ALS, 
                  forestry[forestry@data$GRIDCODE %in% as.character(1997:2010) &
                           forestry@data$action == "cut", ])

## Duplicate ALS data from collection year in 2010 onwards for all years we have 
## nest data for

## Create list for yearly ALS data for which nest data exists and store ALS data
## in it
ALS_by_year <- vector("list", length(c(1998:2004, 2011:2013)))
ALS_by_year[1:length(ALS_by_year)] <- ALS ## Ignore warning!
names(ALS_by_year) <- paste0("y_", c(1998:2004, 2011:2013))

## In the following loop we replace pixels that were cut after collection year 
## 2010 with random samples from clear cuts from 1997 to 2009.

temp <- NULL

for(i in 2011:2013) {
  
  ## Store ALS data used in this loop in temp to make understanding easier
  temp <- ALS_by_year[[paste0("y_", i)]]
  
  ## Which forestry shapes were cut between collection year 2010 and year i
  B1 <- forestry@data$GRIDCODE %in% as.character(2010:i) &
        forestry@data$action == "cut"
  
  ## Replace all pixels not NA which were cut after collection year 2010 with 
  ## random samples from clear cuts from 1997 to 2009.
  if(sum(B1) > 0) {
    
    B2 <- !is.na(mask(temp, forestry[B1, ])[[1]][])
    temp[B2] <- sampleRandom(sample_cc, sum(B2))
    
  } 
  
  ## Update forest are layer with new height data on modelled clear cuts
  temp$area <- !(temp[["height"]][] < 2 | is.na(temp[["height"]][]))
  
  ## Store temp in respective place in ALS_by_year list
  ALS_by_year[[paste0("y_", i)]] <- temp
  
  temp <- NULL
  
} 

## -----------------------------------------------------------------------------

## add no data (-9999) to ALS layers -------------------------------------------
## Replace pixels with -9999 where no data exists. These pixels are forests 
## older than a clear cut or thinnings that occured before data acquisition in 
## 2010 as well as forests after thinnings that happened after 2010. 
## No data pixels bear the value -9999 because NA is used for water.  

## layer names needed for after loop
names <- names(ALS_by_year[[1]])

for(i in unique(nest_pos$year)) {
  
  if(i < 2010) {
    
    ## All cc and thinnings after year i until 2010 are no data
    ## i+1 because it is the interventions in the years after i that matter
    B3 <- forestry@data$GRIDCODE %in% as.character((i+1):2010)
    
  } else {
    
    ## All thinnings after 2010 become no data
    B3 <- forestry@data$GRIDCODE %in% as.character(2010:i) &
          forestry@data$action == "thinned" ## Not correct if th/cc for 2010
    
  }  
  
  if(sum(B3) > 0) {
    
    for(j in 1:length(ALS[1])) {
      
      ALS_by_year[[paste0("y_", i)]][[j]] <- 
        rasterize(forestry[B3, ],
                  ALS_by_year[[paste0("y_", i)]][[j]],
                  update = TRUE, 
                  field = "value")
    
    }
    
  }
  
  names(ALS_by_year[[paste0("y_", i)]]) <- names
  
}

## -----------------------------------------------------------------------------

## Create and add settlement layer to the ALS_by_year --------------------------

## Until 2004, because one settlement disappeard after 2004
dts_u04 <- distanceFromPoints(ALS[[1]], as.matrix(settl[, 1:2]))
names(dts_u04) <- "dts"

dts_a04 <- distanceFromPoints(
  ALS[[1]], as.matrix(settl[settl$Settlement != "Fika_until_2004", 1:2]))
names(dts_a04) <- "dts"

## Add to ALS_by_year
ALS_by_year[paste0("y_", 1998:2004)] <- 
  lapply(ALS_by_year[paste0("y_", 1998:2004)], 
         FUN = function(x) stack(x, dts_u04))
ALS_by_year[paste0("y_", 2011:2013)] <- 
  lapply(ALS_by_year[paste0("y_", 2011:2013)], 
         FUN = function(x) stack(x, dts_a04))

## -----------------------------------------------------------------------------

## Extract ALS data from ALS_by_year for different radiuses around the nest

nest <- as.data.table(nest)

rad <- c(15, seq(50, 450, 50))

around <- nest[, fun_all(.SD), by = 1:nrow(nest)] ## must stand first
#at <- nest[, fun_nest(.SD), by = 1:nrow(nest)]

#at$sample_rad <- paste0("rad_", s_rad)
#at$focal_rad <- paste0("rad_", rad)

around$sample_rad <- paste0("rad_", s_rad)