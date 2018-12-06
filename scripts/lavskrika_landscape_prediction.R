## Make predictions of the probability of reproductive success for
## Arvidsjaur, Lycksele and Åsele kommun and compare it with bird survey data
## provided by bradter et al. 2018.
## In this script the prediction formula used in QGIS is verified.
##
## First edit: 20181114
## Last edit: 20181204
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(ggplot2)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(dplyr)
library(Metrics)

## 2. Define or source functions used in this script ---------------------------

## Define the formula for predicting nest success: 
## Results from rep_succ_analysis with 15m rad.

prob_success <- function(x, y) {
  
  p <- 0.1598 + 
    0.8465 * y + ## y = dts_low
    0.3624 * x - ## x = centered(log(vd))
    1.4957 * y * x
  
  return(exp(p)/(1+exp(p)))
  
}

## Define dts categorisation distance:
R <- 1450

## Define focal rad:
rad <- 15

## 3. Load and explore data ----------------------------------------------------

m.fit <- load("data/m.vd0t5_abs_log_c_15.rda")

m.fit.data <- read.csv("data/p.vd0t5_abs_log_c.csv")
head(m.fit.data)

vd_lp <- raster("data/Dens_mosaic_part.tif")*100
## *100 because input between 0 and 1 but formula expects between 0 and 100
vd_lp

dts <- raster("data/dts.tif")
dts

sj_occ <- raster("data/SFBestCellNbhMainland_99TM.tif")
sj_occ

## 4. Compare the predictions with the formula with the predictions from -------
##    predict()

m.fit.data$dts_low[m.fit.data$dts_cat == "high_ca"] <- 0
m.fit.data$dts_low[m.fit.data$dts_cat == "low_ca"] <- 1

## C is needed in prob_success
C <-  mean(log(m.fit.data$vd_0to5_abs), na.rm = TRUE)

D_pred <- cbind("fit_formula" = 
                  prob_success(x = log(m.fit.data$vd_0to5_abs) - C,
                               y = m.fit.data$dts_low),
                m.fit.data)

plot(D_pred$fit, D_pred$fit_formula); abline(0, 1, col = "red")
cor(D_pred$fit, D_pred$fit_formula)

## 5. Manipulate the raster files and predict ns for all Arvidsjaur, Malå,
##    Lycksele and Åsele

## Average vd_lp within radius used for the prediction in the study:

## Define raster with value = 1 for all cells where understorey is not NA:
num_vd_lp <- vd_lp 
num_vd_lp[!is.na(num_vd_lp)] <- 1 

#Define focal filter:
f <- focalWeight(vd_lp, d = rad, type = "circle") 
f[f > 0] <- 1 

#Calculate sum of raster cells in focal filter:
vd_lp_sum <- focal(vd_lp, f, fun = sum, na.rm = TRUE) 

#Calculate number of raster cells in focal filter: 
vd_lp_num <- focal(num_vd_lp, f, fun = sum, na.rm = TRUE)

## Calculate mean within rad:
vd_lp_mean <- vd_lp_sum / vd_lp_num

## Calculate centered log of vd_lp_mean:
vd_lp_log_cent <- log(vd_lp_mean)-C 

## Change the resolution of dts:
dts <- disaggregate(dts, fact = 4)

## Reduce dts to extent of vd_lp:
dts_red <- crop(dts, extent(vd_lp))

## Categorise dts:
dts_red_cat <- dts_red >= R

## Predict p(successful reprduction):
## Overlay fastest. Tested against calc and manual calculation. Results of all
## three the same

t1 <- Sys.time() 
lp_out <- overlay(vd_lp_log_cent, dts_red_cat, fun = prob_success)
t2 <- Sys.time()
t2-t1

dir.create("temp")
writeRaster(lp_out, "temp/lp_out_raw.tif")

## All NA's in lp_out become 0 as they are non-habitat and no reproduction is 
## expected there
lp_out[is.na(lp_out)] <- 0

## Reduce lp_out to range found in study:

range_dts_low <- c(min(m.fit.data$vd_0to5_abs[m.fit.data$dts_low == 1]),
                   max(m.fit.data$vd_0to5_abs[m.fit.data$dts_low == 1]))

range_dts_high <- c(min(m.fit.data$vd_0to5_abs[m.fit.data$dts_low == 0]),
                    max(m.fit.data$vd_0to5_abs[m.fit.data$dts_low == 0]))

## Define rasters that will become NA after the prediction:
NA_low <- vd_lp_mean <= range_dts_low[1] | vd_lp_mean >= range_dts_low[2]
NA_high <- vd_lp_mean <= range_dts_high[1] | vd_lp_mean >= range_dts_high[2]

## Set non study range predictions to 0. Not selected => no reproduction:
lp_out[dts_red_cat == 1 & NA_low] <- 0
lp_out[dts_red_cat == 0 & NA_high] <- 0

## 6. Compare P(succ_repr) with p(Occurence) -----------------------------------

## Resample sj_occ to resolution of lp_out:
sj_occ_red <- crop(sj_occ, extent(dts))
lp_out_resamp <- resample(lp_out, sj_occ_red)

## Compare lp_out with sj_occ:
plot(lp_out_resamp, sj_occ_red); abline(0, 1, col = "red")
layerStats(stack(sj_occ_red, lp_out_resamp), 
           stat = "pearson", 
           na.rm = TRUE)

lm(values(sj_occ_red) ~ values(lp_out_resamp), na.action = "na.exclude") %>% 
  summary(.) %>% capture.output(.) %>% write(., "results/lm_lp.txt")

## -------------------------------END-------------------------------------------
