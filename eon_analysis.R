# Estimate Rt using data generated from Epidemics on Networks package in Python

# required R packages and functions
require(deSolve)  # version 1.40
require(MASS)     # version 7.3-60.2
require(survival) # version 3.6.4
require(dplyr)    # version 1.1.4
require(tidyr)    # version 1.3.1
require(ggplot2)  # version 3.5.1
require(data.table) # version 1.16.0
require(EpiEstim) # version 2.2.4

source("localDSA_functions.R")
source("util2.R")
source("simulation.R")

# Function to estimate Rt with network data
eon_est <- function(dat, begin, end, width, step, obs_end, use_empEIRsurv = FALSE) {

  full_dat <- EPIdat(dat, begin, end)

  tstarts <- seq(begin, end - width, by = step)
  Rt_out <- data.frame("time" = rep(NA, length(tstarts)),
                       "estimate" = rep(NA, length(tstarts)),
                       "upperRt" = rep(NA, length(tstarts)),
                       "lowerRt" = rep(NA, length(tstarts)),
                       "beta" = rep(NA, length(tstarts)),
                       "delta" = rep(NA, length(tstarts)),
                       "gamma" = rep(NA, length(tstarts)),
                       "est_S" = rep(NA, length(tstarts)),
                       "R0" = rep(NA, length(tstarts)),
                       "rt_var" = rep(NA, length(tstarts)))

  names(dat) <- c("X", "id", "Etime", "Itime", "Rtime", "Estat", "Istat", "Rstat")
  
  if (use_empEIRsurv == TRUE) {
    empsurv <- HSsurv(dat)
    EIRsurv <- data.frame(Esurv = summary(empsurv$Esurv, times = tstarts,
                                          data.frame = TRUE)$surv,
                          Ecumhaz = summary(empsurv$Esurv, times = tstarts,
                                            data.frame = TRUE)$cumhaz,
                          Ecumhaz_se = summary(empsurv$Esurv, times = tstarts,
                                               data.frame = TRUE)$std.err,
                          Isurv = summary(empsurv$Isurv, times = tstarts,
                                          data.frame = TRUE)$surv,
                          Icumhaz = summary(empsurv$Isurv, times = tstarts,
                                            data.frame = TRUE)$cumhaz,
                          Icumhaz_se = summary(empsurv$Isurv, times = tstarts,
                                               data.frame = TRUE)$std.err,
                          Rsurv = summary(empsurv$Rsurv, times = tstarts,
                                          data.frame = TRUE)$surv,
                          Rcumhaz = summary(empsurv$Rsurv, times = tstarts,
                                            data.frame = TRUE)$cumhaz,
                          Rcumhaz_se = summary(empsurv$Rsurv, times = tstarts,
                                              data.frame = TRUE)$std.err)
  }
  for (i in 1:length(tstarts)) {
    try({
      tstart <- tstarts[i]
      print(tstart + width)
      dat <- HSsubset(full_dat, tstart = tstart, tstop = tstart + width)
      if (use_empEIRsurv == FALSE) {
        # add init
        DSAest <- DSAmle(dat, empEIRsurv = NULL, method = "L-BFGS-B")
        pvec <- as.numeric(exp(DSAest$point$point))
        beta <- pvec[1]
        lbeta <- as.numeric(DSAest$point$point)[1]
        delta <- pvec[2]
        gamma <- pvec[3]
        lgamma <- as.numeric(DSAest$point$point)[3]
        xrhoE <- pvec[4]
        xrhoI <- pvec[5]
        xrhoR <- pvec[6]
        rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
        rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
        rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
      } else {
        empEIRsurv <- EIRsurv[i, ]
        DSAest <- DSAmle(dat, empEIRsurv = empEIRsurv, method = "L-BFGS-B")
        pvec <- as.numeric(exp(DSAest$point$point))
        beta <- pvec[1]
        lbeta <- as.numeric(DSAest$point$point)[1]
        delta <- pvec[2]
        gamma <- pvec[3]
        lgamma <- as.numeric(DSAest$point$point)[3]
        #inits <- c(log(beta), log(delta), log(gamma), log(rhoE), log(rhoI), log(rhoR))
        rhoE <- as.numeric(DSAest$empEIR$rhoE)
        rhoI <- as.numeric(DSAest$empEIR$rhoI)
        rhoR <- as.numeric(DSAest$empEIR$rhoR)
      }
      
      if (obs_end == TRUE) {
        S_est <- last(SEIRepidemic(beta = beta, delta = delta, gamma = gamma,
                                   rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, tmin = tstart,
                                   tmax = tstart + width, tstep = 0.01)[, "S"])
      } else {
        epi <- SEIRepidemic(beta = beta, delta = delta, gamma = gamma,
                              rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, 
                              tmin = tstart, tmax = tstart + width, 
                              tstep = 0.01)
        S_est <- epi$S[epi$time == (tstart + tstart + width) / 2]
      }
    
      Rt_out[i, "time"] <- tstart + width
      Rt_out[i, "est_S"] <- S_est
      Rt_out[i, c("beta", "delta", "gamma")] <- c(beta, delta, gamma)
      Rt_out[i, "estimate"] <- exp(lbeta - lgamma + log(S_est))
      Rt_out[i, "R0"] <- exp(lbeta - lgamma)

      mlesamp <- DSApred_mlesamp(DSAest, empEIRsurv = empEIRsurv)
      cis <- DSApred_ci(mlesamp, times = seq(tstart, tstart + width, 0.1))
      
      Rt_out[i, c("upperRt", "lowerRt", "rt_var")] <- 
        last(cis$bounds[, c("upperRt", "lowerRt", "Rt_var")])
    })
  }
    Rt_out = Rt_out
}


# Cori estimates
cori_est <- function(cori_gt, cori_incidence, c_cutoff) {
  pairs_sample <- cori_gt[cori_gt$etime_infectee < c_cutoff, ]
  pairs_sample <- pairs_sample[sample(nrow(pairs_sample), n_sens), ]
  
  gt_df <- data.frame("EL" = floor(pairs_sample$etime_infector),
                      "ER" = ceiling(pairs_sample$etime_infector),
                      "SL" = floor(pairs_sample$etime_infectee),
                      "SR" = ceiling(pairs_sample$etime_infectee))
  gt_df[] <- lapply(gt_df, as.integer)
  
  mcmc_control <- make_mcmc_control(burnin = 1000, thin = 10, seed = 13)
  
  config <- make_config(incid = cori_incidence$incidence,
                        method = "si_from_data",
                        si_parametric_distr = "G",
                        mcmc_control = mcmc_control,
                        n1 = 500,
                        n2 = 50,
                        seed = 918)
  
  estimate_R(incid = cori_incidence$incidence, method = "si_from_data",
                         si_data = gt_df, 
                         config = config)
  
}

# # Plotting function
# eon_plots <- function(DSA, Cori, truth, R0, width, pop, ymax) {
#  dsa_dat <- DSA
#  cori_dat <- Cori$R
#  true_dat <- truth
#  
#  annot <- data.frame(
#    x = end - 40,
#    y = ymax - c(0.15, 0.2, 0.25) * (ymax - 0),
#    label = c(
#      paste0("R0 = ", R0),
#      paste0("DSA window = ", width, " days"),
#      paste0("sample = ", pop, " human sensors")
#    )
#  )
#  
#  dsa_dat %>%
#    ggplot() +
#    theme_minimal() +
#    
#    # Cori CI
#    geom_ribbon(data = cori_dat, aes(x = (t_start + t_end) / 2, ymin = `Quantile.0.025(R)`, 
#                                     ymax = `Quantile.0.975(R)`), fill = "gray", alpha = 0.3) +
#    
#    # DSA CI
#    geom_ribbon(aes(x = time, ymin = lowerRt, ymax = upperRt), fill = "#FFDBB5", 
#                 alpha = 0.4) +
#    
#    # true Rt
#    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
#    
#    # Cori estimate
#    geom_smooth(data = cori_dat, aes(x = t_start, y = `Median(R)`, color = "Cori estimate"),
#                 se = FALSE) +
#    geom_line(data = cori_dat, aes(x = (t_start + t_end) / 2, y = `Median(R)`, color = "Cori estimate")) +
#    #geom_vline(xintercept = c_cutoff, linetype = "dotted") +
#    # DSA estimate 
#    #geom_smooth(aes(x = time, y = estimate, color = "DSA estimate"), se = FALSE) +
#    geom_line(aes(x = time, y = estimate, color = "DSA estimate")) +
#    
#    labs(x = "Time (days)", y = expression(R[t]), color = NULL) +
#    scale_color_manual(
#      values = c(
#        "DSA estimate" = "#DA8210",
#        "Cori estimate" = "#006C67"
#      )
#    ) +
#    geom_text(data = annot, aes(x = x, y = y, label = label), hjust = 0) +
#    coord_cartesian(ylim = c(0, ymax)) 
#  }
# 
# eon_plots(DSA = res4, Cori = cori_res, truth = true_sim, width = 8, R0 = 4,
#           pop = n_sens, ymax = 5)
# 
# eon_plots(DSA = res2, Cori = NULL, truth = true_sim2, width = 6, R0 = 2,
#           pop = n_sens, ymax = 3)

# Smoothed estimate using many windows
adaptive_smooth1 <- function(windows_vec = c(2, 4, 6, 8)) {
  
  # estimate Rt at different window sizes
  window_list <- lapply(windows_vec, function(x) {
    eon_est(dat = eon_sample, begin = begin, end = end, width = x, 
            step = 1, obs_end = TRUE, use_empEIRsurv = TRUE) %>%
    dplyr::select(time, estimate, beta, gamma, R0, rt_var) %>%
    dplyr::rename(!!paste0("rt", x) := estimate,
                  !!paste0("beta", x) := beta,
                  !!paste0("gamma", x) := gamma,
                  !!paste0("R0", x) := R0,
                  !!paste0("var", x) := rt_var)
  })
  
  # housekeeping
  rt_all <- Reduce(function(x, y) 
    inner_join(x, y, by = "time", na_matches = "never"), window_list)
  stable <- apply(rt_all[ , grepl("rt|var", names(rt_all))], 1,
                  function(row) all(is.finite(row)))
  idx <- which(stable)
  rt_trim <- rt_all[idx, ]
  
  # inverse variance weights
  rt_cols <- grep("^rt", names(rt_trim), value = TRUE)
  beta_cols <- grep("^beta", names(rt_trim), value = TRUE)
  gamma_cols <- grep("^gamma", names(rt_trim), value = TRUE)
  var_cols <- grep("^var", names(rt_trim), value = TRUE)
  rt_mat <- as.matrix(rt_trim[ , rt_cols])
  var_mat <- as.matrix(rt_trim[ , var_cols])
  
  precision <- rowSums(1/var_mat) # na.rm = T
  rt_comb <- rowSums(rt_mat/var_mat) / precision
  var_comb <- 1/precision
  
  # CIs
  ntimes <- nrow(rt_trim)
  lower <- numeric(ntimes)
  upper <- numeric(ntimes)
  
  for(i in 1:ntimes){
    
    mu <- as.numeric(rt_trim[i, rt_cols])
    lmu <- log(mu)
    vars <- as.numeric(rt_trim[i, var_cols])
    lvars <- vars / mu^2
    Sigma <- diag(vars)
    lSigma <- diag(lvars)
    #browser()
    lsamples <- mvrnorm(4000, mu = lmu, Sigma = lSigma)
    
    weights <- 1/vars
    weights <- weights/sum(weights)
    lrt_sim <- apply(lsamples, 1, function(x) sample(x, size = 1, prob = weights))
    
    lower[i] <- exp(quantile(lrt_sim, 0.025))
    upper[i] <- exp(quantile(lrt_sim, 0.975))
  }
  
  res_df <- 
    data.frame(
    time = rt_trim$time,
    estimate = rt_comb,
    lower = lower,
    upper = upper
  )
  return(list(rt_trim, rt_all, res_df))
}


# Adaptive window plot
adaptive_smooth_plot <- function(DSAsmooth, Cori, truth, R0, pop, ymax) {
  
  dsa_dat <- DSAsmooth
  cori_dat <- Cori
  cori_dat <- Cori$R
  true_dat <- truth
  
  annot <- data.frame(
    x = end - 40,
    y = ymax - c(0.15, 0.2) * (ymax - 0),
    label = c(
      paste0("R0 = ", R0),
      paste0("sample = ", pop, " human sensors")
    )
  )
  
  cori_dat %>%
    ggplot() +
    theme_minimal() +
    
    # DSA CI
    geom_ribbon(data = dsa_dat, aes(x = time, ymin = lower, ymax = upper),
                fill = "#FFDBB5", alpha = 0.4) +
    # Cori CI
    geom_ribbon(aes(x = (t_start + t_end) / 2, ymin = `Quantile.0.025(R)`, 
                ymax = `Quantile.0.975(R)`), fill = "gray", alpha = 0.3) +
    
    # Cori estimate
    #geom_smooth(data = cori_dat, aes(x = t_start, y = `Median(R)`, color = "Cori estimate"),
    #             se = FALSE) +
    geom_line(aes(x = (t_start + t_end) / 2, y = `Median(R)`, color = "Cori estimate")) +
    geom_vline(xintercept = c_cutoff, linetype = "dotted") +
    
    # true Rt
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
    
    # DSA estimate 
    #geom_smooth(aes(x = time, y = estimate, color = "DSA estimate"), se = FALSE) +
    geom_line(data = dsa_dat, aes(x = time, y = estimate, color = "DSA estimate")) +
    
    labs(x = "Time (days)", y = expression(R[t]), color = NULL) +
    scale_color_manual(
      values = c(
        "DSA estimate" = "#DA8210",
        "Cori estimate" = "#006C67"
      )
    ) +
    geom_text(data = annot, aes(x = x, y = y, label = label), hjust = 0) +
    coord_cartesian(ylim = c(0, ymax), xlim = c(0, 200))
}






