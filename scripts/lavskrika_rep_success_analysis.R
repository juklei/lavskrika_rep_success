## In this script we want to explain the probability of succesful reproduction
## with ALS metrics representing understorey and overstorey density as well as
## forest height in interaction with corvid activity, which in turn is 
## representet by the distance from settlements.
## 
## We also want to test at which radius around the nest the forest composition
## matters for reproductive success.

## First edit: 20180602
## Last edit: 20191119

## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(nlme)
library(lme4)
library(MuMIn)
library(blmeco)
library(DHARMa)
library(pscl)
library(boot)
library(ltm)
library(dplyr)
library(car)

## 2. Define or source functions used in this script ---------------------------

## Specify control values for all models in this script:
cont_spec <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

## Calculate absolute difference in rep_succ relative to categorisation of 
## vd_0to5_abs:
rs_thresh <- function(data, indices) {
  
  d <- data[indices,] # allows boot to select sample 
  return(abs(mean(d$rep_succ[d$dts_cat == "high_ca"], na.rm = TRUE)-
               mean(d$rep_succ[d$dts_cat == "low_ca"], na.rm = TRUE)))
  
}

## Calculate absolute difference in rep_succ relative to categorisation of 
## vd_0to5_abs:
vd_thresh <- function(data, indices) {
  
  dd <- data[indices,] # allows boot to select sample 
  return(abs(mean(d$rep_succ[dd$vd_cat == "dense"])-
               mean(d$rep_succ[dd$vd_cat == "open"])))
  
}

## bootstrapping function for the radius comparison:
rsq <- function(data, indices) {
  
  ddd <- data[indices, ]
  m.T <- glmer(rep_succ ~ dts_cat * vd0t5_abs_log_c + hab_qual + area +
                 (1|female_ring) +
                 (1|male_ring) +
                 (1|year),
               family = binomial,
               data = ddd,
               control = cont_spec)
  
  return(r.squaredGLMM(m.T)[1, 1])
  
}

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
head(nest_height)
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

## 5. Make models for different dts categorisation and categorise --------------
##    dts according to result 

## Reduce nest_ALS to one sample radius:
D <- nest_ALS[nest_ALS$sample_rad == 15, ]

## Make loop through different categorisation distances and store results:

r.dts_cat <- NULL
for(i in seq(500, 3500, 50)) {

  ## Categorise dts
  D$dts_cat <- ifelse(D$dts > i, "low_ca", "high_ca")

  m.dts_cat <- glmer(rep_succ ~ dts_cat + hab_qual +
                       (1|female_ring) +
                       (1|male_ring) +
                       (1|year),
                     family = binomial,  
                     data = D,
                     control = cont_spec)
  
  # ## Calculate threshold in rep_succ between high_ca and low_ca: 
  # T1 <- boot(data = D, statistic = rs_thresh, R = 1000)

  ## Store model output for all categorisation distances:
  r.dts_cat <- rbind(r.dts_cat, 
                     cbind(summary(m.dts_cat)$coefficients, 
                           "neg_logLik" = summary(m.dts_cat)$logLik[1],
                           "R2m" = r.squaredGLMM(m.dts_cat)[1, 1],
                           "AIC" = AIC(m.dts_cat),
                           # summary(T1),
                           "dts_cat" = i))

}

dir.create("results")
write.csv(r.dts_cat, "results/dts_cat_results.csv")

## Test model assumptions with DHARMa for chosen dts_cat:

## Select dts_cat with lowest Estimate:
R <- r.dts_cat[r.dts_cat[,"AIC"] == min(r.dts_cat[,"AIC"]), "dts_cat"]

D$dts_cat <- ifelse(D$dts > R, "low_ca", "high_ca")

m.dts_cat <- glmer(rep_succ ~ dts_cat + hab_qual +
                     (1|female_ring) +
                     (1|male_ring) +
                     (1|year),
                   family = binomial,  
                   data = D,
                   control = cont_spec)

test_my_model(m.dts_cat)

## 6. Make the main model testing ALS on reproductive success ------------------
##    for 15m radius around the nest, beacuse most data for that one.

## Categorise dts according to the results above in nest_ALS:
nest_ALS$dts_cat <- ifelse(nest_ALS$dts > R, "low_ca", "high_ca")

## Exclude NA's in data set for use in models below:
DD <- droplevels(na.omit(
  nest_ALS[, -which(colnames(nest_ALS) %in% c("eggs", "hatched"))]))

## Reduce DD to D15

D15 <- DD[DD$sample_rad == 15, ]

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

m.vd0t5_abs_log_c_15 <- glmer(rep_succ ~ dts_cat * vd0t5_abs_log_c + hab_qual +
                                (1|female_ring) +
                                (1|male_ring) +
                                (1|year),
                              family = binomial,
                              data = D15,
                              control = cont_spec)

## Store results:
capture.output(summary(m.vd0t5_abs_log_c_15),
               confint.merMod(m.vd0t5_abs_log_c_15),
               r.squaredGLMM(m.vd0t5_abs_log_c_15),
               test_my_model(m.vd0t5_abs_log_c_15)) %>% 
  write(., "results/rep_succ_vd0t5_abs_log_c.txt")

## Store the model object for predictions:
save(m.vd0t5_abs_log_c_15, file = "data/m.vd0t5_abs_log_c_15.rda")

## Export data set with prediction for figures to data:

newdata <- expand.grid.df(D15[, c("name", "rep_succ", "dts_cat", 
                                  "vd0t5_abs_log_c", "vd_0to5_abs")],
                          data.frame("hab_qual" = levels(D15$hab_qual)))

## You should make a new data frame with all possible combinations and then 
## plot only for unmanaged or managed:
cbind(predict(m.vd0t5_abs_log_c_15,
              newdata = newdata,
              re.form = NA, 
              se.fit = TRUE, 
              type = "response"), 
      newdata) %>% write.csv(., "data/p.vd0t5_abs_log_c.csv", row.names = FALSE)

## 6b) and compare to quadratic, linear and intercep only model; 

## Only the intercept:
m.int <- glmer(rep_succ ~ (1|female_ring) + (1|male_ring) + (1|year),
               family = binomial,
               data = D15,
               control = cont_spec)

## Only the habitat quality:
m.hq <- glmer(rep_succ ~ hab_qual +
                (1|female_ring) +
                (1|male_ring) +
                (1|year),
              family = binomial,
              data = D15,
              control = cont_spec)

## Corvid activity only:
m.dts_cat_15 <- glmer(rep_succ ~ dts_cat + hab_qual +
                        (1|female_ring) +
                        (1|male_ring) +
                        (1|year),
                      family = binomial,
                      data = D15,
                      control = cont_spec)

## Test model assumptions:
test_my_model(m.dts_cat_15)

## With absolute veg density only as a linear predictor:
m.vd0t5_abs_c_15_only <- glmer(rep_succ ~ vd0t5_abs_c + hab_qual +
                                 (1|female_ring) +
                                 (1|male_ring) +
                                 (1|year),
                               family = binomial,
                               data = D15,
                               control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_abs_c_15_only)

## With absolute veg density as a linear predictor in interaction with dts:
m.vd0t5_abs_c_15 <- glmer(rep_succ ~ dts_cat * vd0t5_abs_c + hab_qual +
                            (1|female_ring) +
                            (1|male_ring) +
                            (1|year),
                          family = binomial,
                          data = D15,
                          control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_abs_c_15)

## With absolute veg density as a quadratic predictor:
m.vd0t5_abs_c_poly_15 <- glmer(rep_succ ~ dts_cat * poly(vd0t5_abs_c, 2) +
                                 hab_qual +
                                 (1|female_ring) +
                                 (1|male_ring) +
                                 (1|year),
                               family = binomial,
                               data = D15,
                               control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_abs_c_15)

## List all models above in a model selection table and export results:
model.sel(m.vd0t5_abs_log_c_15,
          m.hq, 
          m.int,
          m.dts_cat_15, 
          m.vd0t5_abs_c_15_only,
          m.vd0t5_abs_c_15, 
          m.vd0t5_abs_c_poly_15) %>% capture.output(.) %>% 
  write(., "results/rep_succ_mod_sel_formula.txt")

## 6c) and compare vd_0to5_rel; 

m.vd0t5_rel_log_c_15 <- glmer(rep_succ ~ dts_cat * vd0t5_rel_log_c + hab_qual +
                                (1|female_ring) +
                                (1|male_ring) +
                                (1|year),
                              family = binomial,
                              data = D15,
                              control = cont_spec)

## Test model assumptions:
test_my_model(m.vd0t5_rel_log_c_15)

## List both models in a model selection table and export results:
model.sel(m.vd0t5_abs_log_c_15, m.vd0t5_rel_log_c_15) %>% capture.output(.) %>%  
  write(., "results/rep_succ_mod_sel_vd0to5.txt")

## 6d) and together with forest height in the same model:

m.all_forest <- glmer(rep_succ ~ dts_cat * 
                        (vd0t5_abs_log_c + height_c + vd5t_log_c) + hab_qual +
                        (1|female_ring) +
                        (1|male_ring) +
                        (1|year),
                      family = binomial,
                      na.action = "na.fail", 
                      data = D15,
                      control = cont_spec)

## Test model assumptions:
test_my_model(m.all_forest)

## Average all models < deltaAICc and export results:
dredge(m.all_forest) %>% model.avg(., subset = delta <= 2) %>% summary(.) %>% 
  capture.output(.) %>% write(., "results/rep_succ_mod_sel_forest_metrics.txt")
  
## 6f) Check the sensitivity of the results relative to dts_cat:

## Make loop through different categorisation distances and store results:

r.sensitivity <- NULL
for(i in seq(500, 3500, 50)) {
  
  ## Categorise dts
  D15$dts_cat <- ifelse(D15$dts > i, "low_ca", "high_ca")
  
  m.sensitivity <- glmer(rep_succ ~ dts_cat * vd0t5_abs_log_c + hab_qual +
                           (1|female_ring) +
                           (1|male_ring) +
                           (1|year),
                         family = binomial,  
                         data = D15,
                         control = cont_spec)
  
  ## Store model output for all categorisation distances:
  r.sensitivity <- rbind(r.sensitivity, 
                     cbind(summary(m.sensitivity)$coefficients, 
                           "R2m" = r.squaredGLMM(m.sensitivity)[1, 1],
                           "AIC" = AIC(m.sensitivity),
                           "dts_cat" = i))
  
}

dir.create("results")
write.csv(r.sensitivity, "results/sensitivity.csv")

## 7. Assess if a natural threshold can be identified for high corvid ----------
##    activity. For this we classify vd_to5_abs into two categories in steps 
##    low to high density and test the effect on rep_success. The threshold 
##    value is the were the effect size ist highest. 

## We want to do this only for high corvid activity because no threshold value
## is expected in low corvid activity:

DD_t <- D15[D15$dts_cat == "low_ca", ]

r.rep_succ_t <- NULL
for(i in seq(5, 20, 0.1)) {
  
  ## Categorise vd_0to5_abs
  DD_t$vd_cat <- ifelse(DD_t$vd_0to5_abs > i, "dense", "open")
  
  m.rep_succ_t <- glmer(rep_succ ~ vd_cat + hab_qual +
                          (1|female_ring) +
                          (1|male_ring) +
                          (1|year),
                        family = binomial,  
                        data = DD_t,
                        control = cont_spec)

  ## Calculate threshold in rep_succ between open and dense: 
  T2 <- boot(data = DD_t, statistic = vd_thresh, R = 1000)

  ## Store model output for all categorisation distances:
  r.rep_succ_t <- rbind(r.rep_succ_t, 
                             cbind(summary(m.rep_succ_t)$coefficients, 
                                   "AIC" = AIC(m.rep_succ_t),
                                   "R2m" = r.squaredGLMM(m.rep_succ_t)[1, 1],
                                   summary(T2),
                                   "vd_cat" = i))
  
}

dir.create("results")
write.csv(r.rep_succ_t, "results/vd_cat_results.csv")

## Test model assumptions with DHARMa for chosen vd_cat:

## Select vd_cat with highest Estimate:
S <- r.rep_succ_t[r.rep_succ_t[, "Estimate"] == max(r.rep_succ_t[, "Estimate"]),
                  "vd_cat"]

## Categorise vd_0to5_abs
DD_t$vd_cat <- ifelse(DD_t$vd_0to5_abs > S, "dense", "open")

m.rep_succ_t <- glmer(rep_succ ~ vd_cat + hab_qual +
                        (1|female_ring) +
                        (1|male_ring) +
                        (1|year),
                      family = binomial,  
                      data = DD_t,
                      control = cont_spec)

test_my_model(m.rep_succ_t)

## 8. Test the logarithmic model with absolute density for all radiuses --------
##    around the nest and extract AICc values. Also calculate the correlation 
##    of the mean vd_0to5_abs of all radiuses with 15m around the nest and test 
##    the correlation of those correlations with the AICc of the radiuses.

## Reduce data set so the same nests are used for all radiuses:

## Which radius has fewest nests?
N1 <- names(sort(table(DD$sample_rad)))[1]

## Which nests are those?
N2 <- DD[DD$sample_rad == N1 & !is.na(DD$height), "name"]

## Select only those nests:
DD_all_rad <- DD[DD$name %in% N2, ]

## Add centered log(vd_0to5_abs) by radius:
DD_all_rad <- as.data.table(DD_all_rad)
DD_all_rad[, "vd0t5_abs_log_c" := log(vd_0to5_abs) - mean(log(vd_0to5_abs)),
          by = "sample_rad"]

## Make a loop through all radiuses, run the log model and store results:

r.all_rad <- NULL
for(i in unique(DD_all_rad$sample_rad)) {

  tryCatch({
    
    print(i)
  
    ## Calculate correlation of vd_0to5_abs at all radiuses with 15m:
    r.cor <- cor(DD_all_rad[DD_all_rad$sample_rad == i, "vd_0to5_abs"],
                 DD_all_rad[DD_all_rad$sample_rad == 15, "vd_0to5_abs"])
  
    m.all_rad <- glmer(rep_succ ~ dts_cat * vd0t5_abs_log_c + area + hab_qual +
                         (1|female_ring) +
                         (1|male_ring) +
                         (1|year),
                       family = binomial,
                       data = DD_all_rad[DD_all_rad$sample_rad == i, ],
                       control = cont_spec)
 
    # ## Make bootstrapping SE's for R2:
    # T3 <- boot(data = as.data.frame(DD_all_rad[DD_all_rad$sample_rad == i, ]),
    #            statistic = rsq,
    #            R = 100)
  
    ## Store model output for all radiuses:
    r.all_rad <- rbind(r.all_rad, 
                       cbind("cor" = r.cor[1],
                             "estimate" = summary(m.all_rad)$coefficients[5, 1],
                             "SE" = summary(m.all_rad)$coefficients[5, 2],
                             "pvalue" = summary(m.all_rad)$coefficients[5, 4],
                             "AIC" = AIC(m.all_rad), 
                             "R2m" = r.squaredGLMM(m.all_rad)[1, 1],
                             # "bootBias_R2m" = summary(T3)$bootBias,
                             # "bootSE_R2m" = summary(T3)$bootSE,
                             "radius" = i))
  
  
  }, error = function(e) {cat("ERROR :",conditionMessage(e), "\n")})
  
}

## Add deltaAIC to results:
r.all_rad <- as.data.frame(r.all_rad)
r.all_rad$deltaAIC <- r.all_rad$AIC - min(r.all_rad$AIC)

## Export r.all_rad to results and data for making figures:
write.csv(r.all_rad, "results/all_rad_results.csv", row.names = FALSE)
write.csv(r.all_rad, "data/all_rad_results.csv", row.names = FALSE)

## -------------------------------END-------------------------------------------
