library(cmdstanr)
set.seed(918)

# Compile the Stan model
mod <- cmdstan_model("rt_estimation.stan")
# Check that it compiled successfully
mod$print()

# number of human sensors 
n_sens <- 2000

# input and format data
eon_dsa <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/configeon_dsa_deg100_r02.csv")

# tidy the data:
# fix initial cases: etimes = 0
eon_dsa$etime[eon_dsa$itime == 0] <- 0
start <- 0
begin <- 0
end <- ceiling(max(eon_dsa$rtime, na.rm = TRUE))
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
full_dat <- EPIdat(eon_sample, begin, end)

# make the times mapping 
width = 6
step = 1
tstarts <- seq(begin, end - width, by = step)
times <- tstarts
unique_times <- sort(unique(c(full_dat$Etime, full_dat$Itime, full_dat$Rtime)))
Etime_idx <- match(full_dat$Etime, full_dat$Etime)

# empirical ch
empsurv <- HSsurv(full_dat)
EIRcumhaz <- data.frame(Esurv = summary(empsurv$Esurv, times = tstarts,
                                        data.frame = TRUE)$surv,
                        Ecumhaz = summary(empsurv$Esurv, times = tstarts,
                                          data.frame = TRUE)$cumhaz,
                        Esurv_se = summary(empsurv$Esurv, times = tstarts,
                                           data.frame = TRUE)$std.err,
                        Isurv = summary(empsurv$Isurv, times = tstarts,
                                        data.frame = TRUE)$surv,
                        Icumhaz = summary(empsurv$Isurv, times = tstarts,
                                          data.frame = TRUE)$cumhaz,
                        Isurv_se = summary(empsurv$Isurv, times = tstarts,
                                           data.frame = TRUE)$std.err,
                        Rsurv = summary(empsurv$Rsurv, times = tstarts,
                                        data.frame = TRUE)$surv,
                        Rcumhaz = summary(empsurv$Rsurv, times = tstarts,
                                          data.frame = TRUE)$cumhaz,
                        Rsurv_se = summary(empsurv$Rsurv, times = tstarts,
                                           data.frame = TRUE)$std.err)


rt_data <- list(N = n_sens, start = start, end = end, width = width, t0 = 0, 
                ltimes = length(times),
                Etime = full_dat$Etime, Itime = full_dat$Itime, 
                Rtime = full_dat$Rtime, Estat = full_dat$Estat,
                Istat = full_dat$Istat, Rstat = full_dat$Rstat,
                Erisk = full_dat$Erisk, Irisk = full_dat$Irisk,
                Rrisk = full_dat$Rrisk, Etime_idx = Etime_idx, Tmax = end, 
                times = times, use_empEIR = TRUE, emp_Esurv = EIRcumhaz$Esurv, 
                emp_Isurv = EIRcumhaz$Isurv, emp_Rsurv = EIRcumhaz$Rsurv,
                emp_Esurv_se = EIRcumhaz$Esurv_se, emp_Isurv_se = EIRcumhaz$Isurv_se,
                emp_Rsurv_se = EIRcumhaz$Rsurv_se)

fit_mcmc <- mod$sample(
  data = rt_data,
  seed = 123,
  chains = 2,
  parallel_chains = 2
)
