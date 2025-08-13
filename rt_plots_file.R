# Code to accompany the manuscript "Estimation of the effective reproduction
# number Rt using dynamical survival analysis" by Micaela Richter, 
# Grzegorz Rempala, and Eben Kenah.
# last modified: August 13, 2025
# This code generates all figures in the manuscript

## Required R packages
require(deSolve)  # version 1.40
require(MASS)     # version 7.3-60.2
require(survival) # version 3.6
require(dplyr)    # version 1.1.4
require(tidyr)    # version 1.3.1
require(ggplot2)  # version 3.5.1
require(data.table) # version 1.16.0
require(EpiEstim) # version 2.2.4

set.seed(918) # fixed seed for reproducibility
source("rt_est_functions.R") # source backend functions

## Fig 1: intervention with initial R0 = 2 and tE = tI = 4 ---------------------
c1 <- read.csv("2M_int2_44_Cori.csv")
d1 <- read.csv("2M_int2_44_DSA.csv")
pars1 <- list(delta = 1/4, gamma = 1/4, width = 8, end = 150)

fig1 <- list(Rt_out = d1,
             Cori_out = c1,
             parameters = pars1)

rt_plot(fig1, ymax = 3) + coord_cartesian(ylim = c(0, 3))

## Fig 2: intervention with initial R0 = 2 and tE = tI = 4, window = 4 ---------
c2 <- read.csv("2M_int2_44_w4_Cori.csv")
d2 <- read.csv("2M_int2_44_w4_DSA.csv")
pars2 <- list(delta = 1/4, gamma = 1/4, width = 4, end = 150)

fig2 <- list(Rt_out = d2,
             Cori_out = c2,
             parameters = pars2)

rt_plot(fig2, ymax = 3) + coord_cartesian(ylim = c(0, 3))

## Fig 3: no intervention, initial R0 = 2 and tE = 1, tI = 7 -------------------
c3 <- read.csv("2M_r02_17_Cori.csv")
d3 <- read.csv("2M_r02_17_DSA.csv")
pars3 <- list(delta = 1, gamma = 1/7, width = 8, end = 150)

fig3 <- list(Rt_out = d3,
             Cori_out = c3,
             parameters = pars3)

rt_plot(fig3, ymax = 3) + coord_cartesian(ylim = c(0, 3))

## Fig 4: no intervention, initial R0 = 4 and tE = tI = 4 ----------------------
c4 <- read.csv("2M_r04_44_Cori.csv")
d4 <- read.csv("2M_r04_44_DSA.csv")
pars4 <- list(delta = 1/4, gamma = 1/4, width = 8, end = 100)

fig4 <- list(Rt_out = d4,
             Cori_out = c4,
             parameters = pars4)

rt_plot(fig4, ymax = 4.5) + coord_cartesian(ylim = c(0, 4.5))

### Appendix figures

## Fig 5: intervention with initial R0 = 2 and tE = 1, I = 7 -------------------
c5 <- read.csv("2M_int2_17_Cori.csv")
d5 <- read.csv("2M_int2_17_DSA.csv")
pars5 <- list(delta = 1, gamma = 1/7, width = 8, end = 150)

fig5 <- list(Rt_out = d5,
             Cori_out = c5,
             parameters = pars5)

rt_plot(fig5, ymax = 3) + coord_cartesian(ylim = c(0, 3))

## Fig 6: intervention with initial R0 = 2 and tE = 1, I = 7, window = 4 -------
c6 <- read.csv("2M_int2_17_w4_Cori.csv")
d6 <- read.csv("2M_int2_17_w4_DSA.csv")
pars6 <- list(delta = 1, gamma = 1/7, width = 4, end = 150)

fig6 <- list(Rt_out = d6,
             Cori_out = c6,
             parameters = pars6)

rt_plot(fig6, ymax = 3) + coord_cartesian(ylim = c(0, 3))

## Fig 7: intervention with initial R0 = 4 and tE = tI = 4 ---------------------
c7 <- read.csv("2M_int4_44_Cori.csv")
d7 <- read.csv("2M_int4_44_DSA.csv")
pars7 <- list(delta = 1/4, gamma = 1/4, width = 8, end = 100)

fig7 <- list(Rt_out = d7,
             Cori_out = c7,
             parameters = pars7)

rt_plot(fig7, ymax = 5) + coord_cartesian(ylim = c(0, 5))

## Fig 8: intervention with initial R0 = 4 and tE = tI = 4, window = 4 ---------
c8 <- read.csv("2M_int4_44_w4_Cori.csv")
d8 <- read.csv("2M_int4_44_w4_DSA.csv")
pars8 <- list(delta = 1/4, gamma = 1/4, width = 4, end = 100)

fig8 <- list(Rt_out = d8,
             Cori_out = c8,
             parameters = pars8)

rt_plot(fig8, ymax = 5) + coord_cartesian(ylim = c(0, 5))

## Fig 9: no intervention with initial R0 = 2 and tE = tI = 4 ------------------
c9 <- read.csv("2M_r02_44_Cori.csv")
d9 <- read.csv("2M_r02_44_DSA.csv")
pars9 <- list(delta = 1/4, gamma = 1/4, width = 8, end = 150)

fig9 <- list(Rt_out = d9,
             Cori_out = c9,
             parameters = pars9)

rt_plot(fig9, ymax = 3) + coord_cartesian(ylim = c(0, 3))

