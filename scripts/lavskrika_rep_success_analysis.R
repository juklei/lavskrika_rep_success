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
require(lmerTest)

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


## -------------------------------END-------------------------------------------
