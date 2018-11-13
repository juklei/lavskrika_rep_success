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
library(car)

## 2. Define or source functions used in this script ---------------------------

## Specify control values for all models in this script:
cont_spec <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

## Test model assumptions:
test_my_model <- function(m.out) {
  
  sim <- simulateResiduals(m.out)
  plot(sim)
  print(testUniformity(sim))
  print(testZeroInflation(sim))
  print(testDispersion(sim))
  print("Variance Inflation Factor:")
  print(vif(m.out))
  
}

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

  m.dts_cat <- glmer(rep_succ ~ dts_cat +
                       (1|female_ring) +
                       (1|male_ring) +
                       (1|hab_qual) +
                       (1|year),
                     family = binomial,  
                     data = D1,
                     control = cont_spec)

  ## Store model output for all categorisation distances:
  r.dts_cat <- rbind(r.dts_cat, 
                     cbind(summary(m.dts_cat)$coefficients, 
                           "AIC" = AIC(m.dts_cat),
                           "dts_cat" = i))

}

## Test model assumptions with DHARMa for chosen dts_cat:

## Select dts_cat with lowest p value:
R <- r.dts_cat[r.dts_cat[,"Pr(>|z|)"] == min(r.dts_cat[,"Pr(>|z|)"]), "dts_cat"]

D1$dts_cat <- ifelse(D1$dts > R, "low_ca", "high_ca")

m.dts_cat <- glmer(rep_succ ~ dts_cat +
                     (1|female_ring) +
                     (1|male_ring) +
                     (1|hab_qual) +
                     (1|year),
                   family = binomial,  
                   data = D1,
                   control = cont_spec)

test_my_model(m.dts_cat)

dir.create("results")
write.csv(r.dts_cat, "results/dts_cat_results.csv")

## Categorise dts according to the results above in nest_ALS:
nest_ALS$dts_cat <- ifelse(nest_ALS$dts > R, "low_ca", "high_ca")

## 6. Make the main model testing ALS on reproductive success ------------------
##    for 15m radius around the nest, beacuse most data for that one.

D15 <- nest_ALS[nest_ALS$sample_rad == 15 & !is.na(nest_ALS$height), ]

## Add logarithmic version of vd_0to5
D15$vd0t5_rel_log <- log(D15$vd_0to5_rel)
D15$vd0t5_abs_log <- log(D15$vd_0to5_abs)
D15$vd_5to_log <- log(D15$vd_5to)

## Center all continuous variables to avoid covariate correlations:
D15$vd0t5_rel_log_c <- D15$vd0t5_rel_log - mean(D15$vd0t5_rel_log)
D15$vd0t5_abs_log_c <- D15$vd0t5_abs_log - mean(D15$vd0t5_abs_log)
D15$vd0t5_abs_c <- D15$vd_0to5_abs - mean(D15$vd_0to5_abs)
D15$vd0t5_rel_c <- D15$vd_0to5_rel - mean(D15$vd_0to5_rel)
D15$vd5t_log_c <- D15$vd_5to_log - mean(D15$vd_5to_log)
D15$height_c <- D15$height - mean(D15$height)

## 6a) Test log model with absolute veg density as a logarithmic predictor; 

m.vd0t5_abs_log_c <- glmer(rep_succ ~ dts_cat * vd0t5_abs_log_c + 
                             (1|female_ring) +
                             (1|male_ring) +
                             (1|hab_qual) +
                             (1|year) +
                             offset(area),
                           family = binomial,
                           data = D15,
                           control = cont_spec)

## Store results:
capture.output(summary(m.vd0t5_abs_log_c),
               r.squaredGLMM(m.vd0t5_abs_log_c),
               test_my_model(m.vd0t5_abs_log_c)) %>% 
  write(., "results/rep_succ_vd0t5_abs_log_c.txt")

## Store the model object for predictions for figures:
save(m.vd0t5_abs_log_c, file = "data/m.vd0t5_abs_log_c.rda")

## 6b) and compare to quadratic, linear and intercep only model; 

## Only the intercept:
m.int <- glmer(rep_succ ~ 
                 (1|female_ring) +
                 (1|male_ring) +
                 (1|hab_qual) +
                 (1|year),
               family = binomial,
               data = D15,
               control = cont_spec)

## Corvid activity only:
m.dts_cat_15 <- glmer(rep_succ ~ dts_cat +
                        (1|female_ring) +
                        (1|male_ring) +
                        (1|hab_qual) +
                        (1|year),
                      family = binomial,
                      data = D15,
                      control = cont_spec)

## Test model assumptions:
test_my_model(m.dts_cat_15)

## With absolute veg density as a linear predictor:
m.vd0t5_abs_c_15 <- glmer(rep_succ ~ dts_cat * vd0t5_abs_c +
                            (1|female_ring) +
                            (1|male_ring) +
                            (1|hab_qual) +
                            (1|year) +
                            offset(area),
                          family = binomial,
                          data = D15,
                          control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_abs_c_15)

## With absolute veg density as a quadratic predictor:
m.vd0t5_abs_c_poly_15 <- glmer(rep_succ ~ dts_cat * poly(vd0t5_abs_c, 2) +
                                 (1|female_ring) +
                                 (1|male_ring) +
                                 (1|hab_qual) +
                                 (1|year) +
                                 offset(area),
                               family = binomial,
                               data = D15,
                               control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_abs_c_15)

## List all models above in a model selection table and export results:
model.sel(m.vd0t5_abs_log_c,
          m.int, 
          m.dts_cat_15, 
          m.vd0t5_abs_c_15, 
          m.vd0t5_abs_c_poly_15) %>% capture.output(.) %>% 
  write(., "results/rep_succ_mod_sel_formula.txt")

## 6c) and compare vd_0to5_rel; 

m.vd0t5_rel_log_c <- glmer(rep_succ ~ dts_cat * vd0t5_rel_log_c + 
                             (1|female_ring) +
                             (1|male_ring) +
                             (1|hab_qual) +
                             (1|year) +
                             offset(area),
                           family = binomial,
                           data = D15,
                           control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_rel_log_c)

## List both models in a model selection table and export results:
model.sel(m.vd0t5_abs_log_c, m.vd0t5_rel_log_c) %>% capture.output(.) %>%  
  write(., "results/rep_succ_mod_sel_vd0to5.txt")

## 6d) and together with forest height in the same model:

## NAs need to be removed from data set if dredge() is used
D15_red <- na.omit(D15[, -which(colnames(D15) %in% c("eggs", "hatched"))])

m.all_forest <- glmer(rep_succ ~ dts_cat * 
                        (vd0t5_abs_log_c + height_c + vd5t_log_c) +
                        (1|female_ring) +
                        (1|male_ring) +
                        (1|hab_qual) +
                        (1|year) +
                        offset(area),
                      family = binomial,
                      na.action = "na.fail", 
                      data = D15_red,
                      control = cont_spec)

## Test model assumptions:
test_my_model(m.all_forest)

## Average all models < deltaAICc and export results:
dredge(m.all_forest) %>% model.avg(., subset = delta <= 2) %>% summary(.) %>% 
  capture.output(.) %>% write(., "results/rep_succ_mod_sel_forest_metrics.txt")
  
## 7. Test the logarithmic model with absolute density for all radiuses --------
##    around the nest and extract AICc values. Also calculate the correlation 
##    of the mean vd_0to5_abs of all radiuses with 15m around the nest and test 
##    the correlation of those correlations with the AICc of the radiuses.

## Reduce data set so the same nests are used for all radiuses:

#...

## -------------------------------END-------------------------------------------
