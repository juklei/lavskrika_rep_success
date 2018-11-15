## We want to know if SJ chose nests strategically or opportunistically.
## For this we calculate the difference in the understorey density
## between 2ha around the nest with the understorey density within 66ha,
## the territory size. If positive, denser nests than available are chosen,
## If negative, more open. The explanatory variables are availability 
## (mean understorey density at territory level) and corvid activity.
## 
##
## First edit: 20180810
## Last edit: 20181115

## Author: Julian Klein

## Clear environment and load libraries ----------------------------------------

rm(list = ls())

library(nlme)
library(lme4)
library(blmeco)
library(pscl)
library(boot)
library(ltm)
library(lmerTest)
library(car)
library(dplyr)

## Define or source functions used in this script ------------------------------

## dts categorisation radius:
R <- 1500

## Radius of 2 ha area relevant around the nest:
A <- 80

## Load and explore data -------------------------------------------------------

ALS <- read.csv("data/ALS_rep_succ.csv")
nest <- read.csv("data/nest_reproduction.csv")

head(ALS)
head(nest)

## 4. Merge and prepare data set for analysis ----------------------------------

## Reduce ALS data to the radiuse relevant around the nest and territory radius:
ALS_red <- merge(ALS[ALS$sample_rad == A, ], 
                 ALS[ALS$sample_rad == max(ALS$sample_rad), ],
                 by = c("name", "territory", "year"))

## Calculate the difference between vd_0to5_abs within A and within territory:
ALS_red$vd_diff <- ALS_red$vd_0to5_abs.x - ALS_red$vd_0to5_abs.y  

## Exlcude experimental data from nest:
nest <- nest[nest$excl_nest_pl == "no", 
             c("name", "female_ring", "male_ring", "hab_qual")]

## Merge reproduction with ALS and nest height data:
np_ALS <- na.omit(merge(ALS_red, nest, by = "name"))

## We want to categorise vd_0to5_abs within the territory into open and dense:
np_ALS$vd_terr_cat <- ifelse(np_ALS$vd_0to5_abs.y > mean(np_ALS$vd_0to5_abs.y),
                             "dense", "open")

## Categorise dts according to findings in rep_success_analysis:
np_ALS$dts_cat <- ifelse(np_ALS$dts.x > R, "low_ca", "high_ca")

## Make factors:
np_ALS$female_ring <- as.factor(np_ALS$female_ring)
np_ALS$male_ring <- as.factor(np_ALS$male_ring)
np_ALS$year <- as.factor(np_ALS$year)

## Center vd_0to5_abs.y:
np_ALS$vd_terr_cont_cent <- np_ALS$vd_0to5_abs.y - mean(np_ALS$vd_0to5_abs.y)

## 5. Make a model, export results and make predictions for figures:

m.nest_pos <- lmer(vd_diff ~  vd_terr_cont_cent * dts_cat +
                (1|male_ring) +
                (1|female_ring) +
                (1|year) +
                (1|hab_qual),
              data = np_ALS)

## Test model assumptions and export results:
qqnorm(summary(m.nest_pos)$residuals); qqline(summary(m.nest_pos)$residuals)

## Store results:
capture.output(summary(m.nest_pos),
               r.squaredGLMM(m.nest_pos),
               "vif:" =vif(m.nest_pos)) %>% 
  write(., "results/nest_pos_terr_cont.txt")

## Export data set with prediction for figures to data:
cbind(predict(m.nest_pos, 
              re.form = NA, 
              se.fit = TRUE),
              np_ALS[, c("vd_diff", "dts_cat", "vd_0to5_abs.y")]) %>% 
  write.csv(., "data/p.nest_pos.csv", row.names = FALSE)

## ------------------------------END--------------------------------------------