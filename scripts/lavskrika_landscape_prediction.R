## Make predictions of the probability of reproductive success for
## Arvidsjaur, Lycksele and Åsele kommun and compare it with bird survey data
## provided by bradter et al. 2018.
## In this script this is done for a smaller area for trial. 
## The actual calculation is done in a different software.
##
## First edit: 20181114
## Last edit: 20181119
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

## 2. Define or source functions used in this script ---------------------------

## Define the formula for predicting nest success: 
## Results from rep_succ_analysis

prob_succsess <- function(data) {
  
  p <- 0.1598 + 
    0.8465 * data$dts_low + 
    0.3624 * (log(data$vd_0to5_abs)-C) - 
    1.4957 * data$dts_low * (log(data$vd_0to5_abs)-C)
  
  return(exp(p)/(1+exp(p)))
  
}

## 3. Load and explore data ----------------------------------------------------

p.vd0t5_abs_log_c <- read.csv("data/p.vd0t5_abs_log_c.csv")

## 4. Compare the predictions with the formula with the predictions from -------
##    predict()

p.vd0t5_abs_log_c$dts_low[p.vd0t5_abs_log_c$dts_cat == "high_ca"] <- 0
p.vd0t5_abs_log_c$dts_low[p.vd0t5_abs_log_c$dts_cat == "low_ca"] <- 1

## C is needed in prob_succsess
C <-  mean(log(p.vd0t5_abs_log_c$vd_0to5_abs), na.rm = TRUE)

D_pred <- cbind("fit_formula" = prob_succsess(data = p.vd0t5_abs_log_c),
                p.vd0t5_abs_log_c)

plot(D_pred$fit, D_pred$fit_formula); cor(D_pred$fit, D_pred$fit_formula)

## -------------------------------END-------------------------------------------