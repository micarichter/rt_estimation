# MLE in stan to improve behavior & speed 
# script to accompany rt_window.stan 
# last edited 11.5.2026

library(cmdstanr)
set.seed(918)

# number of human sensors 
n_sens <- 2000

# read in true data
true_sim <- read.csv("true_deg100_r02.csv")

# input and format data
eon_dsa <- read.csv("dat_deg100_r02.csv")

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


# compile the stan model
mod <- cmdstan_model("rt_window.stan")

eon_est_stan <- function(dat, begin, end, width, step, EIRcumhaz) {
  
  full_dat <- EPIdat(dat, begin, end)
  tstarts <- seq(begin, end - width, by = step)
  n_times <- length(tstarts)
  
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

  Rt_out <- data.frame(
    time = numeric(n_times),
    estimate = numeric(n_times),
    lowerRt = numeric(n_times),
    upperRt = numeric(n_times),
    beta = numeric(n_times),
    gamma = numeric(n_times),
    S_est = numeric(n_times),
    rt_var = numeric(n_times)
  )
  
  for (i in seq_along(tstarts)) {
    tstart <- tstarts[i]
    dat_i <- HSsubset(full_dat, tstart = tstart, tstop = tstart + width)
    
    t0 <- tstart - 0.001
    times <- sort(unique(c(dat_i$Etime, dat_i$Itime, dat_i$Rtime)))
    times <- times[is.finite(times) & times > t0]
    times <- sort(unique(c(times, tstart + width)))
    
    Etime_idx <- sapply(dat_i$Etime, function(x) {
      if (is.na(x) || is.infinite(x)) return(length(times))
      if (x <= t0) return(1)
      idx <- which.min(abs(times - x))
      if (length(idx) == 0) return(length(times)) 
      return(idx)
    })

    stan_data <- list(
      N = nrow(dat_i),
      t0 = t0,
      ltimes = length(times),
      times = times,
      Etime = dat_i$Etime, Itime = dat_i$Itime, Rtime = dat_i$Rtime,
      Estat = dat_i$Estat, Istat = dat_i$Istat, Rstat = dat_i$Rstat,
      Erisk = dat_i$Erisk, Irisk = dat_i$Irisk, Rrisk = dat_i$Rrisk,
      Etime_idx = Etime_idx,
      use_empEIR = 0,
      emp_Esurv = EIRcumhaz$Esurv[i], emp_Esurv_se = EIRcumhaz$Esurv_se[i],
      emp_Isurv = EIRcumhaz$Isurv[i], emp_Isurv_se = EIRcumhaz$Isurv_se[i],
      emp_Rsurv = EIRcumhaz$Rsurv[i], emp_Rsurv_se = EIRcumhaz$Rsurv_se[i]
    )
    
    # Stan Optimizer + Laplace Approximation
    fit_draws <- tryCatch({
      # Step A: Fast MLE Optimization
      fit_mle <- mod$optimize(
        data = stan_data, 
        seed = 123, 
        jacobian = TRUE,
        algorithm = "lbfgs",
        refresh = 0 # Suppresses terminal output
      )
      # step b: laplace approx
      lap <- mod$laplace(
        data = stan_data,
        mode = fit_mle,
        draws = 4000, # Number of draws to estimate variance/CIs
        refresh = 0
      )
      dummy_check <- lap$draws("Rt_end")
      lap # return laplace object if successful
    }, error = function(e) { 
    
    # step c: if laplace fails use pathfinder (variational inference)
    tryCatch({
      pf <- mod$pathfinder(
        data = stan_data,
        draws = 4000,
        seed = 123,
        refresh = 0
      )
      return(pf)
    }, error = function(e) {
    return(NULL)
  })
})
    if (!is.null(fit_draws)) {
      # laplace object supports the exact same $summary() calls as MCMC
      sum_rt <- fit_draws$summary("Rt_end")
      sum_beta <- fit_draws$summary("beta")
      sum_gamma <- fit_draws$summary("gamma")
      sum_S <- fit_draws$summary("S_end")
      
      # Calculate variance
      q_low <- if ("q5" %in% names(sum_rt)) sum_rt$q5 else sum_rt$`5%`
      q_high <- if ("q95" %in% names(sum_rt)) sum_rt$q95 else sum_rt$`95%`
      rt_var <- (sum_rt$sd)^2
      
      Rt_out[i, ] <- c(
        time = tstart + width,
        estimate = sum_rt$median,
        lowerRt = sum_rt$q5,
        upperRt = sum_rt$q95,
        beta = sum_beta$median,
        gamma = sum_gamma$median,
        S_est = sum_S$median,
        rt_var = rt_var
      )
    } else {
      Rt_out[i, ] <- c(tstart + width, rep(NA, 7))
    }
  }
  return(Rt_out)
}

res <- eon_est_stan(eon_sample, 0, 100, 8, 1, EIRcumhaz = EIRcumhaz)


plot(res$time, res$estimate, type = "l")
lines(true_sim$time, true_sim$true_rt, col = "red")
grid()



