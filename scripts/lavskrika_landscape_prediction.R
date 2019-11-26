## Make predictions of the probability of reproductive success for
## Arvidsjaur, Lycksele and Åsele kommun and compare it with bird survey data
## provided by bradter et al. 2018.
## In this script the prediction formula used in QGIS is verified.
##
## First edit: 20181114
## Last edit: 20191120
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
  
  p <- 0.2411 + 
    1.0427 * y + ## y = dts_low
    0.3965 * x - ## x = centered(log(vd))
    0.5196/2 - ## /2 for the mean between unmanaged/managed
    1.5703 * y * x
  
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

vd_lp <- raster("data/Dens_05_5_mosaic.tif")*100
## *100 because input between 0 and 1 but formula expects between 0 and 100
vd_lp

dts <- raster("data/dts.tif")
dts

sj_occ <- raster("data/SFBestCellNbhMainland_99TM.tif")
sj_occ

lp_shape <- shapefile("data/lp_shape.shp")

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

plot(D_pred$fit_formula, D_pred$fit); abline(0, 1, col = "red")
cor(D_pred$fit, D_pred$fit_formula) ## The dicrepancy is because the formula 
                                    ## takes the mean of managed/unmanaged, 
                                    ## while the model fit data has predictions 
                                    ## once for managed and once for unmanaged. 

## 5. Manipulate the raster files and predict ns for all Arvidsjaur, Malå,
##    Lycksele and Åsele

## Reduce vd_lp to the extent of dts:
vd_lp <- crop(vd_lp, extent(dts))

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

## Categorise dts:
dts_cat <- dts >= R

## Predict p(successful reprduction):
## Overlay fastest. Tested against calc and manual calculation. Results of all
## three the same

gc()
t1 <- Sys.time() 
lp_out <- overlay(vd_lp_log_cent, dts_cat, fun = prob_success)
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
lp_out[dts_cat == 1 & NA_low] <- 0
lp_out[dts_cat == 0 & NA_high] <- 0

## 6. Compare P(succ_repr) with p(Occurence) -----------------------------------

## Resample sj_occ to resolution of lp_out:
lp_out_resamp <- resample(lp_out, sj_occ)

## Reduce both to the Kommun borders of Arvidsjaur, Åsele. Malå and Lycksele
lp_out_resamp <- mask(lp_out_resamp, lp_shape)
sj_occ <- mask(sj_occ, lp_shape)

## Temporary_export:
dir.create("temp")
writeRaster(lp_out_resamp, "temp/lp_out_resamp.tif")

## Compare lp_out with sj_occ:

dir.create("results")
capture.output(
  
  "Av" = cor.test(mask(sj_occ, lp_shape[lp_shape$KnNamn == "Arvidsjaur", ])[],
                  mask(lp_out_resamp, 
                       lp_shape[lp_shape$KnNamn == "Arvidsjaur", ])[],
                  method = "spearman"),
  "Ma" = cor.test(mask(sj_occ, lp_shape[lp_shape$KnNamn == "Malå", ])[],
                  mask(lp_out_resamp, lp_shape[lp_shape$KnNamn == "Malå", ])[],
                  method = "spearman"),
  "Ly" = cor.test(mask(sj_occ, lp_shape[lp_shape$KnNamn == "Lycksele", ])[],
                  mask(lp_out_resamp, 
                       lp_shape[lp_shape$KnNamn == "Lycksele", ])[],
                  method = "spearman"),
  "As" = cor.test(mask(sj_occ, lp_shape[lp_shape$KnNamn == "Åsele", ])[],
                  mask(lp_out_resamp, lp_shape[lp_shape$KnNamn == "Åsele", ])[],
                  method = "spearman")
  
) %>% write(., "results/lm_lp_rho.txt")

## Look at disagreement:

lp_rank <- (lp_out_resamp - min(lp_out_resamp[], na.rm = TRUE))/
           (max(lp_out_resamp[], na.rm = TRUE) - min(lp_out_resamp[], na.rm = TRUE))
sjo_rank <- (sj_occ - min(sj_occ[], na.rm = TRUE))/
            (max(sj_occ[], na.rm = TRUE) - min(sj_occ[], na.rm = TRUE))

disag <- sjo_rank - lp_rank ## Use direction of difference

dir.create("results")
capture.output(
  print("Arvidsjaur"),
  mean(mask(disag, lp_shape[lp_shape$KnNamn == "Arvidsjaur", ])[], 
       na.rm = TRUE),
  print("Malå"),
  mean(mask(disag, lp_shape[lp_shape$KnNamn == "Malå", ])[], na.rm = TRUE),
  print("Lycksele"),
  mean(mask(disag, lp_shape[lp_shape$KnNamn == "Lycksele", ])[], na.rm = TRUE),
  print("Åsele"),
  mean(mask(disag, lp_shape[lp_shape$KnNamn == "Åsele", ])[], na.rm = TRUE)

) %>% write(., "results/disagreement_15.txt")

dir.create("temp")
writeRaster(disag, "temp/disag_15.tif")

## Compare 15 m with 83 m:

cor.test(na.omit(disag15[]), na.omit(disag83[]))

## Bin disagree:

bdg <- rbind(expand.grid("disag" = na.omit(mask(disag, lp_shape[lp_shape$KnNamn == "Arvidsjaur", ])[]), 
                         "county" = "Arvidsjaur"),
             expand.grid("disag" = na.omit(mask(disag, lp_shape[lp_shape$KnNamn == "Malå", ])[]), 
                         "county" = "Malå"),
             expand.grid("disag" = na.omit(mask(disag, lp_shape[lp_shape$KnNamn == "Lycksele", ])[]), 
                         "county" = "Lycksele"),
             expand.grid("disag" = na.omit(mask(disag, lp_shape[lp_shape$KnNamn == "Åsele", ])[]), 
                         "county" = "Åsele"))

png("figures/figure_lsp_dens.png", 10000/4, 6000/4, "px", res = 600/4)

ggplot(bdg) + 
  geom_vline(xintercept = c(0.25, 0, -0.25, -0.5, -0.75), linetype = "dashed", color = "grey") +
  # geom_histogram(aes(x = disag, color = county, fill = county, alpha = 0.5)) +
  # geom_histogram(aes(x = disag), color = "black", fill = "white") +
  geom_density(aes(x = disag, color = county, linetype = county), size = 2) +
  # facet_grid(county ~ ., scales = "free_y") +
  scale_x_continuous(breaks = c(0.25, 0, -0.25, -0.5, -0.75)) +
  xlab("disagreement values") + ylab("") + 
  labs(color = "counties from north to south", 
       linetype = "counties from north to south") +
  theme_classic(40) +
  # theme(legend.position = "none")
  theme(legend.position = c(0.75, 0.8), legend.key.size = unit(3, 'lines'))

dev.off()

## -------------------------------END-------------------------------------------
