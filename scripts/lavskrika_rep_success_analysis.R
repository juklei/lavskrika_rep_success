## Describe the purpose of the script
## 

## First edit: 20180602
## Last edit: 20181113

## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(ggplot2)
library(data.table)
library(reshape2)
library(nlme)
library(lme4)
library(MuMIn)
library(blmeco)
library(DHARMa)
library(pscl)
library(boot)
library(ltm)
library(lmerTest)
library(dplyr)

## 2. Define or source functions used in this script ---------------------------

source()

## 3. Load and explore data ----------------------------------------------------

dir("data")

ALS <- read.csv("data/ALS_rep_succ.csv")
nest_height <- na.omit(read.csv("data/nest_heights.csv"))
nest <- read.csv("data/nest_reproduction.csv")

head(ALS)
head(nest_heights)
head(nest)

## 4. Merge and prepare data set for analysis ----------------------------------

## Exlcude experimental data from nest and add rep success yes/no because nest 
## predation by corvids is absolute :
nest <- nest[nest$excl_rep_suc == "no", ]
nest$rep_succ <- ifelse(nest$fledl > 0, 1, 0)

## Merge reproduction with ALS and nest height data:
nest_ALS <- merge(nest, ALS, by = c("name", "territory", "year"))
nest_height <- merge(nest, nest_height, all.x = TRUE, by = "name")

## Make factors:
nest_ALS$female_ring <- as.factor(nest_ALS$female_ring)
nest_ALS$male_ring <- as.factor(nest_ALS$male_ring)
nest_ALS$year <- as.factor(nest_ALS$year)

## 5. Make models for different dts categorisation -----------------------------

## Reduce nest_ALS to one sample radius:
D1 <- nest_ALS[nest_ALS$sample_rad == 15, ]

## Make loop through different categorisation distances and store results:

r.dts_cat <- NULL
for(i in seq(1000, 2000, 100)) {

  ## Categorise dts
  D1$dts_cat <- ifelse(D1$dts > i, "low_ca", "high_ca")

  m.dts_cat <- glmer(rep_succ ~  dts_cat +
                       (1|female_ring) +
                       (1|male_ring) +
                       (1|hab_qual) +
                       (1|year),
                     family = binomial,  
                     data = D1,
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 100000)))

  ## Store model output for all categorisation distances:
  r.dts_cat <- rbind(r.dts_cat, 
                     cbind(summary(m.dts_cat)$coefficients, 
                           "AIC" = summary(m.dts_cat)$AIC["AIC"],
                           "dts_cat" = i))

}

## Test model assumptions with DHARMa for chosen dts_cat:

## Select dts_cat with lowest p value:
R <- r.dts_cat[r.dts_cat[,"Pr(>|z|)"] == min(r.dts_cat[,"Pr(>|z|)"]), "dts_cat"]

D1$dts_cat <- ifelse(D1$dts > R, "low_ca", "high_ca")

m.dts_cat <- glmer(rep_succ ~  dts_cat +
                     (1|female_ring) +
                     (1|male_ring) +
                     (1|hab_qual) +
                     (1|year),
                   family = binomial,  
                   data = D1,
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 100000)))

sim <- simulateResiduals(m.dts_cat)
plot(sim)
testUniformity(sim)
testZeroInflation(sim)

dir.create("results")
write.csv(r.dts_cat, "results/dts_cat_results.csv")

## 6. Make the main model testing ALS on reproductive success ------------------

## Categorise dts according to the results above in nest_ALS:
nest_ALS$dts_cat <- ifelse(nest_ALS$dts > R, "low_ca", "high_ca")





## -------------------------------END-------------------------------------------
