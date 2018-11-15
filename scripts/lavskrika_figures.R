## In this script all figures for the SJ rep_sucess paper are produced
## 
##
## First edit: 
## Last edit: 
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(ggplot2)
library(data.table)
library(dplyr)

## 2. Define or source functions used in this script ---------------------------

source()

## 3. Load and explore data ----------------------------------------------------

dir("data")
D_fig1 <- read.csv("data/p.vd0t5_abs_log_c.csv")
D_fig2 <- read.csv("data/all_rad_results.csv")
#D_fig3 <- read.csv("data/")
head(D_fig1)
head(D_fig2)

## 4. Produce Figure 1 showing the main results for the 15m radius -------------

## Change dts_cat level names:
levels(D_fig1$dts_cat) <- c("high corvid activity", "low corvid activity") 

p1 <- ggplot(D_fig1, aes(x = vd_0to5_abs, 
                         y = rep_succ, 
                         fill = dts_cat, 
                         color = dts_cat))
p2 <- geom_jitter(size = 4, na.rm = TRUE, width = 0, height = 0.01)
p3 <- geom_line(aes(x = vd_0to5_abs, y = fit, lty = dts_cat), size = 3)
p4 <- geom_ribbon(aes(ymin = fit - (1.96*se.fit), ymax = fit + (1.96*se.fit)),
#                  colour = NA, 
                  alpha = 0.2)

p1 + p2 + p3 + p4 +
  scale_color_grey(start = 0.1, end = 0.6) + 
  scale_fill_grey(start = 0.1, end = 0.6) +
  ylab("probability of succesful reproduction") + 
  xlab("understorey density") + 
  theme_classic(45) +                  
  theme(legend.position = c(0.8, 0.1), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

## 5. Produce Figure 2 showing the relationship between vegetation -------------
##    correlation and AIC values for all areas around the nest the 15m radius 

x <- read.csv("rad_P.csv")

c1 <- ggplot(x, aes(x = radius))
c2 <- geom_line(size = 2, aes(y = estimate, colour = "estimate"))
c2b <- geom_ribbon(aes(ymin = estimate - std_error, 
                       ymax = estimate + std_error,
                       colour = "estimate"), 
                   alpha = 0.1)
c2c <- geom_point(size = 5, aes(y = estimate, colour = "estimate"))
c3 <- geom_line(size = 2, 
                linetype = "dashed", 
                aes(y = ((p.value*10)-2.5), colour = "p.value"))
c3b <- geom_point(size = 5, aes(y = ((p.value*10)-2.5), colour = "p.value"))
c4 <- scale_y_continuous(
  sec.axis = sec_axis(~(.+2.5)/10, 
                      name = "p.value of the interaction"))
c5 <- geom_hline(yintercept = (0.05*10)-2.5, linetype = "dashed", color = "grey")

c1 + c2 + c2c + c3 + c3b + c4 + c5 + 
  scale_color_grey() + 
  scale_fill_grey() +
  xlab("radius(m) around the nest") + 
  ylab("effect size of the interaction") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.2, 0.8), 
        legend.title = element_blank(),
        legend.key.size = unit(2.5, 'lines'))


## 6. Produce Figure 3 showing the nest site selection -------------------------

p1 <- ggplot(data = prd_diff, aes(x = dts, 
                                  y = vd_0to5_diff, 
                                  fill = terr_cat,
                                  colour = terr_cat))
p2 <- geom_point(size = 4, na.rm = TRUE, alpha = 0.8)
p2b <- geom_jitter(size = 4, na.rm = TRUE, width = 0, height = 0.01)
p3 <- geom_line(aes(x = dts, y = prd_diff, lty = terr_cat), size = 3)
p4 <- geom_ribbon(aes(ymin = prd_diff - se.fit, 
                      ymax = prd_diff + se.fit),
                  colour = NA, 
                  alpha = 0.2)
p5 <- geom_vline(xintercept = 1500, color = "black", linetype = "dashed")
p6 <- geom_hline(yintercept = 0, color = "black", linetype = "dashed")
# p7 <- geom_abline(intercept = 0, slope = 1, linetype = "dashed")
# p8 <- scale_y_continuous(limits = c(0, 25))

p1 + p5 + p6 + p2 + p3 + p4 +
  scale_color_grey() + 
  scale_fill_grey() +
  ylab("relative understorey density") + 
  xlab("distance(m) to closest settlement") + 
  theme_classic(45) +                  
  theme(legend.position = c(0.75, 0.1), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))


## -----------------------------------------------------------------------------
