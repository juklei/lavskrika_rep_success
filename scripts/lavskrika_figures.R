## In this script all figures for the SJ rep_sucess paper are produced
## 
##
## First edit: 20181101
## Last edit: 20190420
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(ggplot2)
library(data.table)
library(dplyr)

## 2. Define or source functions used in this script ---------------------------

#source()

## 3. Load and explore data ----------------------------------------------------

dir("data")

D_fig1 <- read.csv("data/p.vd0t5_abs_log_c.csv")
D_fig2 <- read.csv("data/all_rad_results.csv")
D_fig3 <- read.csv("data/p.nest_pos.csv")
head(D_fig1)
head(D_fig2)
head(D_fig3)

## 4. Produce Figure 1 showing the main results for the 15m radius -------------

## Change dts_cat level names:
levels(D_fig1$dts_cat) <- c("close to human settlement", 
                            "far from human settlement") 

p1 <- ggplot(D_fig1, aes(x = vd_0to5_abs, 
                         y = rep_succ, 
                         fill = dts_cat, 
                         color = dts_cat))
p2 <- geom_jitter(size = 4, na.rm = TRUE, width = 0, height = 0.01)
p3 <- geom_line(aes(x = vd_0to5_abs, y = fit, lty = dts_cat), size = 3)
p4 <- geom_ribbon(aes(ymin = pmax(fit - (1*se.fit), 0), 
                      ymax = pmin(fit + (1*se.fit), 1)),
                  colour = NA, 
                  alpha = 0.1)
p5 <- annotate("text", 10, 0.15, label = "R2m = 0.08 / R2c = 0.17", size = 10)

## Store figures:

dir.create("figures")

png("figures/lavskrika_F1.png", 10000/4, 7000/4, "px", res = 600/4)

p1 + p2 + p3 + p4 +
  scale_color_manual(breaks = c("close to human settlement", 
                                "far from human settlement"),
                     values=c("red", "blue")) + 
  scale_fill_manual(breaks = c("close to human settlement", 
                               "far from human settlement"),
                    values=c("red", "blue")) +
  ylab("probability of succesful reproduction") + 
  xlab("understory density") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.8, 0.2), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## 5. Produce Figure 2 showing the relationship between vegetation -------------
##    correlation and AIC values for all areas around the nest the 15m radius 

## Add actual area (ha) to data set:
D_fig2$area <- round(D_fig2$radius^2*pi/10000, digits = 2)

## Calculate cor values used for vertical lines below:
vert_lin <- D_fig2$cor[D_fig2$radius %in% c(78, 202, 450)]

c1 <- ggplot(D_fig2, aes(x = cor, y = R2m))
c2a <- geom_point(size = 4)
c2b <- geom_errorbar(aes(ymax = R2m + bootSE_R2m, ymin = R2m - bootSE_R2m),
                     color = "grey")
c2c <- geom_ribbon(aes(ymax = R2m + (1.96*bootSE_R2m),
                       ymin = pmax(R2m - (1.96*bootSE_R2m), 0)),
                   colour = NA, 
                   alpha = 0.1)
c3 <- scale_x_reverse()
c4a <- geom_vline(xintercept = vert_lin, color = "black", linetype = "dashed") 
c4b <- geom_text(aes(x = vert_lin[2] + 0.025, label = "450 m", y = 0.115), 
                 angle = 0,
                 size = 10)
c4c <- geom_text(aes(x = vert_lin[3] + 0.021, label = "80 m", y = 0.115), 
                 angle = 0,
                 size = 10)
c4d <- geom_text(aes(x = vert_lin[1] + 0.025, label = "200 m", y = 0.115), 
                 angle = 0,
                 size = 10)

png("figures/lavskrika_F2.png", 10000, 6000, "px", res = 600)

c1 + c4a + c4b + c4c + c4d + c2a + c3 +
  scale_color_grey(start = 0.1, end = 0.5) + 
  scale_fill_grey(start = 0.1, end = 0.5) +
  xlab("correlation of ud at nest with ud at radius i around the nest") + 
  ylab("r-squared with ud at radius i") + 
  theme_classic(40)

dev.off()

## 6. Produce Figure 3 showing the nest site selection -------------------------

## Change dts_cat level names:
levels(D_fig3$dts_cat) <- c("high corvid activity", "low corvid activity") 

q1 <- ggplot(D_fig3, aes(x = vd_0to5_abs.y, 
                         y = vd_diff, 
                         fill = dts_cat, 
                         color = dts_cat))
q2 <- geom_point(size = 4, na.rm = TRUE)
q3 <- geom_line(aes(x = vd_0to5_abs.y, y = fit, lty = dts_cat), size = 3)
q4 <- geom_ribbon(aes(ymin = fit - (1.96*se.fit), 
                      ymax = fit + (1.96*se.fit)),
                  colour = NA, 
                  alpha = 0.1)

## Store figures:

dir.create("figures")

png("figures/lavskrika_F3.png", 10000, 7000, "px", res = 600)

q1 + q2 + q3 + q4 +
  scale_color_grey(start = 0.1, end = 0.5) + 
  scale_fill_grey(start = 0.1, end = 0.5) +
  ylab("realtive understorey density nest") + 
  xlab("understorey density territory") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.8, 0.2), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()


## -----------------------------------------------------------------------------
