# Code to create csvs in the data folder
# Note that running this file will take several days

# required R packages
require(deSolve)  # version 1.40
require(MASS)     # version 7.3-60.2
require(survival) # version 3.6.4
require(dplyr)    # version 1.1.4
require(tidyr)    # version 1.3.1
require(ggplot2)  # version 3.5.1
require(data.table) # version 1.16.0
require(EpiEstim) # version 2.2.4

source("localDSA_functions.R")
source("rt_est_functions.R")
source("util2.R")
source("simulation.R")

set.seed(918)
pop <- 2e6
name <- "2M"

# scenario 1: intervention with initial R0 = 2, delta = 1/4, gamma = 1/4
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1/4, gamma = 1/4, R0 = 2\n")
  system.time(rt_int2_44 <- rt_est(n = pop, begin = 0, end = 150, width = 8, 
                                   step = 1, R0 = 2, R0_int = 0.8, R0_end = 1.15, 
                                   delta = 1/4, gamma = 1/4, I_init = 3e-5 * pop, 
                                   int_time1 = 60, int_time2 = 90, days1 = 7, days2 = 7))

  write.csv(rt_int2_44$Rt_out, paste0(name, "_int2_44_DSA.csv"))
  write.csv(rt_int2_44$Cori_out, paste0(name, "_int2_44_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 2: intervention with  R0 = 2, delta = 1/4, gamma = 1/4, window = 4
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1/4, gamma = 1/4, R0 = 2, w = 4\n")
  system.time(rt_int2_44w4 <- rt_est(n = pop, begin = 0, end = 150, width = 4, step = 1, 
                                 R0 = 2, R0_int = 0.8, R0_end = 1.15, delta = 1/4, 
                                 gamma = 1/4, I_init = 3e-5 * pop, int_time1 = 60, 
                                 int_time2 = 90, days1 = 7, days2 = 7))

  write.csv(rt_int2_44w4$Rt_out, paste0(name, "_int2_44_w4_DSA.csv"))
  write.csv(rt_int2_44w4$Cori_out, paste0(name, "_int2_44_w4_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 3: delta = 1, gamma = 1/7, R0 = 2
tryCatch({
  cat("\nrunning scenario: int = N, delta = 1, gamma = 1/7, R0 = 2\n")
  system.time(rt_2_17 <- rt_est(n = pop, begin = 0, end = 150, width = 8, step = 1,
                                R0 = 2, delta = 1, gamma = 1/7, I_init = 3e-5 * pop))
  
  write.csv(rt_2_17$Rt_out, paste0(name, "_r02_17_DSA.csv"))
  write.csv(rt_2_17$Cori_out, paste0(name, "_r02_17_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})


# scenario 4: delta = 1/4, gamma = 1/4, R0 = 4
tryCatch({
  cat("\nrunning scenario: int = N, delta = 1/4, gamma = 1/4, R0 = 4\n")
  system.time(rt_4 <- rt_est(n = pop, begin = 0, end = 100, width = 8, step = 1,
                             R0 = 4, delta = 1/4, gamma = 1/4, I_init = 3e-5 * pop))
  
  write.csv(rt_4$Rt_out, paste0(name, "_r04_44_DSA.csv"))
  write.csv(rt_4$Cori_out, paste0(name, "_r04_44_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 5: intervention with delta = 1, gamma = 1/7, R0 = 2
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1, gamma = 1/7, R0 = 2\n")
  system.time(rt_int2 <- rt_est(n = pop, begin = 0, end = 150, width = 8, step = 1,
                                R0 = 2, R0_int = 0.8, R0_end = 1.15, delta = 1,
                                gamma = 1/7, I_init = 3e-5 * pop, int_time1 = 60,
                                int_time2 = 90, days1 = 7, days2 = 7))
  
  write.csv(rt_int2$Rt_out, paste0(name, "_int2_17_DSA.csv"))
  write.csv(rt_int2$Cori_out, paste0(name, "_int2_17_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 6: intervention with delta = 1, gamma = 1/7, R0 = 2, window = 4
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1, gamma = 1/7, R0 = 2, w = 4\n")
  system.time(rt_int2w4 <- rt_est(n = pop, begin = 0, end = 150, width = 4, step = 1,
                                R0 = 2, R0_int = 0.8, R0_end = 1.15, delta = 1,
                                gamma = 1/7, I_init = 3e-5 * pop, int_time1 = 60,
                                int_time2 = 90, days1 = 7, days2 = 7))

  write.csv(rt_int2w4$Rt_out, paste0(name, "_int2_17_w4_DSA.csv"))
  write.csv(rt_int2w4$Cori_out, paste0(name, "_int2_17_w4_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})


# scenario 7: intervention with delta = 1/4, gamma = 1/4, R0 = 4
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1/4, gamma = 1/4, R0 = 4\n")
  system.time(rt_int444 <- rt_est(n = pop, begin = 0, end = 100, width = 8, step = 1, 
                                  R0 = 4, R0_int = 1.6, R0_end = 2.3, delta = 1/4, 
                                  gamma = 1/4, I_init = 3e-5 * pop, int_time1 = 20, 
                                  int_time2 = 30, days1 = 7, days2 = 7))

  write.csv(rt_int444$Rt_out, paste0(name, "_int4_44_DSA.csv"))
  write.csv(rt_int444$Cori_out, paste0(name, "_int4_44_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 8: intervention with delta = 1/4, gamma = 1/4, R0 = 4, window = 4
tryCatch({
  cat("\nrunning scenario: int = Y, delta = 1/4, gamma = 1/4, R0 = 4, w = 4\n")
system.time(rt_int444w4 <- rt_est(n = pop, begin = 0, end = 100, width = 4, step = 1, 
                                R0 = 4, R0_int = 1.6, R0_end = 2.3, delta = 1/4, 
                                gamma = 1/4, I_init = 3e-5 * pop, int_time1 = 20, 
                                int_time2 = 30, days1 = 7, days2 = 7))

  write.csv(rt_int444w4$Rt_out, paste0(name, "_int4_44_w4_DSA.csv"))
  write.csv(rt_int444w4$Cori_out, paste0(name, "_int4_44_w4_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})

# scenario 9: delta = 1/4, gamma = 1/4, R0 = 2
tryCatch({
  cat("\nrunning scenario: int = N, delta = 1/4, gamma = 1/4, R0 = 2\n")
  system.time(rt_2 <- rt_est(n = pop, begin = 0, end = 150, width = 8, step = 1,
                             R0 = 2, delta = 1/4, gamma = 1/4, I_init = 3e-5 * pop))
  
  p3 <- rt_plot(rt_2$Rt_out, rt_2$Cori_out, n = pop, end = 150, R0 = 2, delta = 1/4,
                gamma = 1/4, I_init = 3e-5 * pop) + coord_cartesian(ylim = c(0, 10))
  
  ggsave(paste0(name, "_r02_44.png"), plot = p3, width = 6, height = 4, dpi = 300,
         bg = "white")
  write.csv(rt_2$Rt_out, paste0(name, "_r02_44_DSA.csv"))
  write.csv(rt_2$Cori_out, paste0(name, "_r02_44_Cori.csv"))
}, error = function(e) {
  message("oops!:", conditionMessage(e))
})




