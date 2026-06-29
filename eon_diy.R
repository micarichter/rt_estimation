# EoN analysis
source("eon_analysis.R")

# number of human sensors 
n_sens <- 1500

# human sensor network from EoN
eon_dsa <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/configeon_dsa_deg100_r02.csv")

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
true_sim <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/config_deg100_r02.csv")

# Cori estimates
cori_gt <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/cori_config_pairs_deg100_r02.csv")
cori_incidence <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/config_incidence_deg100_r02.csv")
names(cori_incidence) <- c("time", "incidence")

# find point at which S = 0.9
s_prop <- 0.9
c_cutoff <- min(true_sim$time[true_sim$S/n_tot == s_prop])
cori_res <- cori_est(cori_gt, cori_incidence, c_cutoff)


# one window size
system.time(res2 <- eon_est(dat = eon_sample, begin = 40, end = 100, width = 6, 
                            step = 1, obs_end = TRUE, use_empEIRsurv = TRUE))

# several windows smoothed
short <- adaptive_smooth1(c(2, 4, 6))
long <- adaptive_smooth1(c(8, 10, 12))


# plot
adaptive_smooth_plot(DSA = windows_deg10_4a[[3]], Cori = cori_res, truth = true_sim, R0 = 4,
                     pop = n_sens, ymax = 5)


windows_deg100_2 <- adaptive_smooth1(c(4, 6, 8, 10))

adaptive_smooth_plot(DSA = windows_deg100_2[[3]], Cori = cori_res, truth = true_sim, R0 = 3,
                     pop = n_sens, ymax = 4)

