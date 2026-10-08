# EoN analysis
source("eon_analysis.R")

# number of human sensors 
n_sens <- 2000

# human sensor network from EoN
eon_dsa <- read.csv("dat_deg100_r02.csv")

# tidy the data:
# fix initial cases: etimes = 0
eon_dsa$etime[eon_dsa$itime == 0] <- 0
start <- 0
begin <- 0
end <- max(eon_dsa$rtime, na.rm = TRUE) + 1
n_tot <- dim(eon_dsa)[1]

# everyone else with etime = NA is censored; they stayed in S forever
# these people should Inf event times
eon_dsa$etime <- ifelse(is.na(eon_dsa$etime), Inf, eon_dsa$etime)
eon_dsa$itime <- ifelse(eon_dsa$etime == Inf, Inf, eon_dsa$itime)
eon_dsa$rtime <- ifelse(eon_dsa$etime == Inf, Inf, eon_dsa$rtime)

# add indicators
eon_dsa$estat <- ifelse(eon_dsa$etime < end, 1, 0)
eon_dsa$istat <- ifelse(eon_dsa$itime < end, 1, 0)
eon_dsa$rtsat <- ifelse(eon_dsa$rtime < end, 1, 0)

names(eon_dsa) <- c("X", "id", "eTime", "iTime", "rTime", "estat", "istat", "rstat")

# take a sample of size = n_sens (or read in previously selected sample)
eon_sample <- eon_dsa[sample(nrow(eon_dsa), n_sens), ]

# true Rt data
true_sim <- read.csv("true_deg100_r02.csv")

# Cori estimates
cori_gt <- read.csv("cori_pairs_deg100_r02.csv")
cori_incidence <- read.csv("cori_incidence_deg100_r02.csv")
names(cori_incidence) <- c("time", "incidence")

# find point at which S = 0.9
# s_prop <- 0.9
# c_cutoff <- min(true_sim$time[true_sim$S/n_tot == s_prop])
# cori_res <- cori_est(cori_gt, cori_incidence, end)

# Cori oracle estimator
# delta_true <- 1
# gamma_true <- 1/7
# parlist <- {
#   list(
#     t_E = 1 / delta_true,
#     t_I = 1 / gamma_true
#   )
# }
# parlist$true_mean_GI = (parlist$t_E + parlist$t_I)
# parlist$true_var_GI = parlist$t_E^2 + parlist$t_I^2

#cori_res1 <- get_cori(cori_incidence, icol_name = 'incidence', window = 1)

# Cori parametric 
#cori_gt$SI <- cori_gt$etime_infectee - cori_gt$etime_infector
cori_para <- coriC(5, end, 1, cori_gt, cori_incidence)

# one window size
system.time(res2 <- eon_est(dat = eon_sample, begin = 0, end = 150, width = 10, 
                            step = 1, obs_end = TRUE, use_empEIR = TRUE,
                            CIs = TRUE))                                                                                                                  
# system.time(res2a <- eon_est(dat = eon_sample, begin = 10, end = 14, width = 2, 
#                             step = 1, obs_end = TRUE, use_empEIR = TRUE,
#                             CIs = TRUE))   

plot(res2a$time, res2a$estimate, type = "l", ylim = c(0, 3), xlim = c(0, 150))
lines(true_sim$time, true_sim$true_rt)
lines(res2$time, res2$upperRt, col = "green")
lines(res2$time, res2$lowerRt, col = "green")

# ggplot(data = res2) +
#   geom_line(aes(x = time, y = estimate), color = "#DA8210") +
#   geom_ribbon(aes(x = time, ymin = lowerRt, ymax = upperRt), fill = "#FFDBB5", alpha = 0.4) +
#   geom_line(data = true_sim, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, 5), xlim = c(0, 100))

# several windows smoothed 

# run windows in parallel
library(progressr)
handlers(global = TRUE) 
handlers("cli")

plan(multisession, workers = 4)
windows_deg100_2bf <- adaptive_smooth1(c(4, 6, 8, 10), CIs = TRUE)
plan(sequential)

plot(windows_deg100_2bf[[3]]$time, windows_deg100_2bf[[3]]$estimate, type = "l", na.rm = T)

adaptive_smooth_plot(DSA = windows_deg100_2bf[[3]], Cori = cori_para, truth = true_sim, R0 = 2,
                     pop = n_sens, ymax = 3)
