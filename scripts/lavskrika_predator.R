## Make a model for CJ activity
## 
## First edit: 20191113
## Last edit: 20191113
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(lme4)
library(MuMIn)
library(DHARMa)
library(ltm)
library(ggplot2)
library(magrittr)
library(MASS)

## 2. Define or source functions used in this script ---------------------------

## Test model assumptions:
test_my_model <- function(m.out) {
  sim <- simulateResiduals(m.out)
  plot(sim)
  print(testUniformity(sim))
  print(testZeroInflation(sim))
  print(testDispersion(sim))
}

## 3. Load and explore data ----------------------------------------------------

dir("data")
pred <- read.csv("data/corvid_data_2005.csv")
head(pred)

## 4. Plot data and build a GLM ------------------------------------------------

m <- glm(cj_presence ~ dts, data = pred, family = binomial) 
capture.output(summary(m), r.squaredGLMM(m), dose.p(m, p = 0.5)) %>% 
  write(., "results/predator.txt")
test_my_model(m)

pred_pred <- predict(m, re.form = NA, se.fit = TRUE, type = "response")

## Calculate LD50:
LD50 <- dose.p(m, p = 0.5)

## Make a plot:
png("figures/figure_pred.png", 10000/4, 7000/4, "px", res = 600/4)

ggplot() + 
geom_line(aes(x = pred$dts, y = pred_pred$fit), size = 2) +
geom_ribbon(aes(x = pred$dts, 
                ymin = pred_pred$fit - pred_pred$se.fit,
                ymax = pred_pred$fit + pred_pred$se.fit),
            alpha = 0.2) +
geom_jitter(aes(x = pred$dts, y = pred$cj_presence), size = 6, height = 0.05) + 
scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2500)) +
geom_vline(xintercept = 1316) + geom_vline(xintercept = c(1165, 1467), 
                                           linetype = "dashed") +
xlab("distance of the nest to the cosest settlement") +
ylab("probability of Eurasian jay presence") +
theme_classic(40)

dev.off()

## -------------------------------END-------------------------------------------
