# the goal of this script is to generate simulation result from mrgsolve
# the simulation result will be compared to Julia simulation 

rm(list=ls())  #Clear out existing objects
gc()

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script

# original model
mod <- mread("../model/human_erythroid_lymphoid_myeloid") %>% param(epsilon_spl = 0.032, delta_dp = 1e-6) 

# pre-equlibrium
pre_sim <- mod  %>% init(LT0 = 1000, LT1 = 100) %>%  mrgsim(end = 600) %>% as_tibble() 

write.csv(pre_sim%>% select(time, LT0, LT1, RBC0, RBC1, GM0, GM1, BMrec0, BMrec1, cd4rec0, cd8rec0, cd4rec1, cd8rec1), "rbc.csv", row.names=FALSE)
