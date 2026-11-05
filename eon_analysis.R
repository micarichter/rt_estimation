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
require(future.apply) # version 4.5.2
require(numDeriv)

source("localDSA_functions.R")
source("util2.R")
source("simulation.R")

# Function to estimate Rt with network data
eon_est <-function(dat, begin, end, width, step, obs_end, use_empEIR = FALSE, 
                   CIs = TRUE, maxit_sann = 5000) {

  full_dat <- EPIdat(dat, begin, end)

  tstarts <- seq(begin, end - width, by = step)
  Rt_out <- data.frame("time" = rep(NA, length(tstarts)),
                       "estimate" = rep(NA, length(tstarts)),
                       "upperRt" = rep(NA, length(tstarts)),
                       "lowerRt" = rep(NA, length(tstarts)),
                       "beta" = rep(NA, length(tstarts)),
                       "delta" = rep(NA, length(tstarts)),
                       "gamma" = rep(NA, length(tstarts)),
                       "rhoE" = rep(NA, length(tstarts)),
                       "rhoI" = rep(NA, length(tstarts)),
                       "rhoR" = rep(NA, length(tstarts)),
                       "est_S" = rep(NA, length(tstarts)),
                       "R0" = rep(NA, length(tstarts)),
                       "rt_var" = rep(NA, length(tstarts)))

  names(dat) <- c("X", "id", "Etime", "Itime", "Rtime", "Estat", "Istat", "Rstat")
  
  if (use_empEIR == TRUE) {
    empsurv <- HSsurv(dat)
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
  }
  init2 <- c(0, 0, 0, 0, 0, 0)
   mle_sol <- FALSE
  for (i in 1:length(tstarts)) {
    tryCatch({
      tstart <- tstarts[i]
      print(tstart + width)
      dat <- HSsubset(full_dat, tstart = tstart, tstop = tstart + width)
      if (use_empEIR == FALSE) {
        DSAest <- DSAmle(dat, empEIR = NULL, method = "SANN",
                         control = list(maxit = maxit_sann, temp = 30, 
                                        tmax = 100))
      } else {
        empEIR <- EIRcumhaz[i, ]
        
        # two step optimization
        # step 1: SANN
        
        #print("step 1")
        if(!mle_sol) {
          # use empEIR until sann finds a solution
          DSAest_sann <- DSAmle(dat, empEIR = empEIR, method = "SANN",
                                control = list(maxit = maxit_sann,
                                               temp = 100, tmax = 30))
        } else {
          DSAest_sann <- DSAmle(dat, method = "SANN", init = init2,
                                control = list(maxit = maxit_sann,
                                               temp = 100, tmax = 30))
        }
          sann_est <- DSAest_sann$point$point
          #print(DSAest_sann$loglik)

        # step 2: Nelder-Mead or L-BFGS-B
        #print("step 2")
        lwr <- c(lnbeta = -4, lndelta = -4, lngamma = -4, 
                 lnxrhoE = -8, lnxrhoI = -8, lnxrhoR = -8)
        upr <- c(lnbeta = 4, lndelta = 4, lngamma = 4, 
                          lnxrhoE = 8, lnxrhoI = 8, lnxrhoR = 8)
        DSAest <- DSAmle(dat, method = "L-BFGS-B", init = sann_est, lower = lwr,
                         upper = upr)
        #print(DSAest$loglik)
        
        # check if BFGS found a solution
        if (!is.null(DSAest$point$point) && any(!is.na(DSAest$point$point))) {
          init2 <- DSAest$point$point
          mle_sol <- TRUE
        }
      }
      
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
    
      if (obs_end == TRUE) {
        S_est <- last(SEIRepidemic(beta = beta, delta = delta, gamma = gamma,
                                   rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, 
                                   tmin = tstart, tmax = tstart + width, 
                                   tstep = 0.01)[, "S"])
      } else {
        epi <- SEIRepidemic(beta = beta, delta = delta, gamma = gamma,
                              rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, 
                              tmin = tstart, tmax = tstart + width, 
                              tstep = 0.01)
        S_est <- epi$S[epi$time == (tstart + tstart + width) / 2]
      }
      
      rt_est <- exp(lbeta - lgamma + log(S_est))
      r0_est <- exp(lbeta - lgamma)
      
      bounds <- c(NA, NA, NA)

      if (CIs == TRUE) {
        mlesamp <- DSApred_mlesamp(DSAest, empEIR = empEIR)
        cis <- DSApred_ci(mlesamp, times = seq(tstart, tstart + width, 0.1))
      }
      Rt_out[i, ] <- list(tstart + width, rt_est, bounds[1], bounds[2], 
                          beta, delta, gamma, rhoE, rhoI, rhoR, 
                          S_est, r0_est, bounds[3])
    }, error = function(e) {
      
    })
  }
    return(Rt_out)
}


# Cori estimates
cori_B <- function(cori_gt, cori_incidence, c_cutoff) {
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


# Cori estimate using parametric SI
# first, calculate the SI 
true_sim <- read.csv("eon_simulations/true_deg100_r02.csv")

# Cori estimates
cori_gt <- read.csv("eon_simulations/cori_pairs_deg100_r02.csv")
cori_incidence <- read.csv("eon_simulations/cori_incidence_deg100_r02 .csv")
names(cori_incidence) <- c("time", "incidence")
cori_gt$SI <- cori_gt$etime_infectee - cori_gt$etime_infector

coriC <- function(begin, end, step, cori_gt, cori_incidence) {
  tstarts <- seq(begin, end, by = step)
  tryCatch({
    cori_est <- data.frame("time" = rep(NA, length(tstarts)),
                         "estimate" = rep(NA, length(tstarts)),
                         "upperRt" = rep(NA, length(tstarts)),
                         "lowerRt" = rep(NA, length(tstarts)))
    for (i in begin:length(tstarts)) {
      tstart <- tstarts[i]
      mean_si <- mean(cori_gt$SI[cori_gt$etime_infectee < i])
      std_si <- sqrt(var(cori_gt$SI[cori_gt$etime_infectee < i]))
    
      coriC <- estimate_R(cori_incidence$incidence, method = "parametric_si",
                          config = make_config(list(
                          mean_si = mean_si, std_si = std_si)))
      cori_est$time[i] <- coriC$R$t_start[i]
      cori_est$estimate[i] <- coriC$R$`Median(R)`[i]
      cori_est$upperRt[i] <- coriC$R$`Quantile.0.025(R)`[i]
      cori_est$lowerRt[i] <- coriC$R$`Quantile.0.975(R)`[i]
    }
  })
    return(cori_est)
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
adaptive_smooth1 <- function(windows_vec = c(2, 4, 6, 8), CIs = TRUE) {
  windows_vec <- sort(windows_vec)
  
  # progress bar
  p <- progressor(along = windows_vec)
  
  # estimate Rt at different window sizes
  window_list <- future_lapply(windows_vec, function(x) {
    res <- eon_est(dat = eon_sample, begin = begin, end = end, width = x, 
            step = 1, obs_end = TRUE, use_empEIR = TRUE, CIs = CIs) |>
    dplyr::select(time, estimate, beta, gamma, R0, rt_var) |>
    dplyr::rename(!!paste0("rt", x) := estimate,
                  !!paste0("beta", x) := beta,
                  !!paste0("gamma", x) := gamma,
                  !!paste0("R0", x) := R0,
                  !!paste0("var", x) := rt_var)
    
    p(sprintf("finished window size %d", x))
    return(res)
  }, future.seed = TRUE, future.packages = c("dplyr", "deSolve"))
  
  # housekeeping
  rt_all <- Reduce(function(x, y)
    left_join(x, y, by = "time", na_matches = "never"), window_list)

  # keep rows that have at least one pair of rt and var
  rt_cols <- grep("^rt", names(rt_all), value = TRUE)
  var_cols <- grep("^var", names(rt_all), value = TRUE)

  rt_mat <- as.matrix(rt_all[, rt_cols])
  var_mat <- as.matrix(rt_all[, var_cols])
  
  stable_rows <- rowSums(is.finite(rt_mat) & is.finite(var_mat)) > 0
  rt_trim <- rt_all[stable_rows, ]
  
  rt_mat <- as.matrix(rt_trim[, rt_cols])
  var_mat <- as.matrix(rt_trim[, var_cols])
  
  # inverse variance weights
  precision <- rowSums(1/var_mat, na.rm = TRUE) # na.rm = T
  rt_comb <- rowSums(rt_mat/var_mat, na.rm = TRUE) / precision
  var_comb <- 1/precision
  
  res_df <- data.frame(
    time = rt_trim$time,
    estimate = rt_comb
  )
  
  # CIs
  if (CIs == TRUE) {
    ntimes <- nrow(rt_trim)
    lower <- numeric(ntimes)
    upper <- numeric(ntimes)
    
    for(i in 1:ntimes){
      
      mu <- as.numeric(rt_trim[i, rt_cols])
      vars <- as.numeric(rt_trim[i, var_cols])
      
      # keep only finite and not NA estimates 
      stable <- is.finite(mu) & is.finite(vars) & mu > 0 & vars > 0
      
      mu_stable <- mu[stable]
      vars_stable <- vars[stable]
      
      lmu <- log(mu_stable)
      lvars <- vars_stable / mu_stable^2
      
      if (length(mu_stable) == 1) {
        lrt_sim <- rnorm(4000, mean = lmu, sd = sqrt(lvars))
      } else {
        weights <- 1 / vars_stable
        weights <- weights / sum(weights)

        # sample mixture component indices then draw from rnorm
        idx <- sample(seq_along(mu_stable), size = 4000, replace = TRUE, prob = weights)
        lrt_sim <- rnorm(4000, mean = lmu[idx], sd = sqrt(lvars[idx]))
      }
      
      lower[i] <- exp(quantile(lrt_sim, 0.025))
      upper[i] <- exp(quantile(lrt_sim, 0.975))
    }
    res_df$lower <- lower
    res_df$upper <- upper
  }
  res_df <- subset(res_df, time >= max(windows_vec))
  return(list(rt_trim, rt_all, res_df))
}


# Adaptive window plot
adaptive_smooth_plot <- function(DSAsmooth, Cori = NULL, truth, ymax) {
  
  dsa_dat <- DSAsmooth
  cori_dat <- Cori
  #cori_dat <- Cori$R
  true_dat <- truth
  
  # annot <- data.frame(
  #   x = end - 40,
  #   y = ymax - c(0.15, 0.2) * (ymax - 0),
  #   label = c(
  #     paste0("R0 = ", R0),
  #     paste0("sample = ", pop, " human sensors")
  #   )
  # )
  
  cori_dat %>%
    ggplot() +
    theme_minimal() +
    
    # DSA CI
    geom_ribbon(data = dsa_dat, aes(x = time, ymin = lower, ymax = upper),
                fill = "#FFDBB5", alpha = 0.4) +
    # Cori CI
    geom_ribbon(aes(x = time, ymin = lowerRt, 
                    ymax = upperRt), fill = "gray", alpha = 0.3) +
    # geom_ribbon(aes(x = (t_start + t_end) / 2, ymin = `Quantile.0.025(R)`, 
    #             ymax = `Quantile.0.975(R)`), fill = "gray", alpha = 0.3) +
    
    # Cori estimate
    #geom_smooth(data = cori_dat, aes(x = t_start, y = `Median(R)`, color = "Cori estimate"),
    #             se = FALSE) +
    #geom_line(aes(x = (t_start + t_end) / 2, y = `Median(R)`, color = "Cori estimate")) +
    geom_line(aes(x = time, y = estimate, color = "Cori estimate")) +
    #geom_vline(xintercept = c_cutoff, linetype = "dotted") +
    
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
    # geom_text(data = annot, aes(x = x, y = y, label = label), hjust = 0) +
    coord_cartesian(ylim = c(0, ymax), xlim = c(0, 150))
}


# Adaptive smooth plot 2
adaptive_smooth_plot1 <- function(DSAsmooth, Cori = NULL, truth, ymax) {
  
  dsa_dat <- DSAsmooth
  true_dat <- truth
  
  plot <- dsa_dat %>%
    ggplot() +
    theme_minimal() +
    
    # DSA CI
    geom_ribbon(data = dsa_dat, aes(x = time, ymin = lower, ymax = upper),
                fill = "#FFDBB5", alpha = 0.4) +
    # DSA estimate 
    geom_line(data = dsa_dat, aes(x = time, y = estimate, color = "DSA estimate")) +
    # true Rt
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", linetype = "dashed")
  
  if (!is.null(Cori)) {
    cori_dat <- Cori
    # Cori CI
    plot <- plot + 
      geom_ribbon(data = cori_dat, aes(x = time, ymin = lowerRt, ymax = upperRt), 
                  fill = "gray", alpha = 0.3) +
      # Cori point estimate
      geom_line(data = cori_dat, aes(x = time, y = estimate, color = "Cori estimate")) +
      scale_color_manual(
        values = c(
          "DSA estimate" = "#DA8210",
          "Cori estimate" = "#006C67"))
  } else {
    plot <- plot +
      scale_color_manual(
      values = c("DSA estimate" = "#DA8210"))
  }
    plot <- plot + labs(x = "Time (days)", y = expression(R[t]), color = NULL) +
      coord_cartesian(ylim = c(0, ymax), xlim = c(0, 150))
    return(plot)
}


# Cori variations plot
cori_plots <- function(coriA, coriB, coriC, true_dat, R0, pop, ymax) {
  
  pairs_res <- coriB$R
  # annot <- data.frame(
  #   x = end - 40,
  #   y = ymax - c(0.15, 0.2) * (ymax - 0),
  #   label = c(
  #     paste0("R0 = ", R0),
  #     paste0("sample = ", pop, " human sensors")
  #   )
  #)
  
  coriA %>%
    ggplot() +
    theme_minimal() +
  
    # CoriA (oracle estimator)
    geom_ribbon(aes(x = time, ymin = Cori.025, 
                    ymax = Cori.975), fill = "gray", alpha = 0.3) +
    geom_line(aes(x = time, y = Cori.mean, color = "Oracle estimate")) +
    
    # true Rt
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
    
    # CoriB (all pairs)
    geom_ribbon(data = pairs_res, aes(x = t_start, ymin = `Quantile.0.025(R)`, 
                                ymax = `Quantile.0.975(R)`), fill = "#FFDBB5", alpha = 0.5) +
    geom_line(data = pairs_res, aes(x = t_start, y = `Median(R)`, color = "Cori pairs data")) +
    
    # CoriC (most fair)
    geom_ribbon(data = coriC, aes(x = time, ymin = lowerRt, ymax = upperRt),
                fill = "palegreen", alpha = 0.3) +
    geom_line(data = coriC, aes(x = time, y = estimate, color = "Chronological Cori")) +
    
    labs(x = "Time (days)", y = expression(R[t]), color = NULL) +
    scale_color_manual(
      values = c(
        "Cori pairs data" = "#DA8210",
        "Oracle estimate" = "#00008B",
        "Chronological Cori" = "#06402B" 
      )
    ) + coord_cartesian(ylim = c(0, ymax), xlim = c(0, 150))
}

# Cori panel plot
cori_plots1 <- function(coriA, coriB, coriC, true_dat, ymax) {
  
  pairs_res <- coriB$R
  common_layers <- list(
    theme_minimal(),
    coord_cartesian(ylim = c(0, ymax), xlim = c(0, 150)),
    theme(legend.position = "none")
  )
  
  # Panel 1: CoriA (Oracle estimate)
  p1 <- ggplot() +
    common_layers +
    geom_ribbon(data = coriA, aes(x = time, ymin = Cori.025, ymax = Cori.975), 
                fill = "gray", alpha = 0.3) +
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", 
              linetype = "dashed") +
    geom_line(data = coriA, aes(x = time, y = Cori.mean), color = "#006C67") +
    labs(x = NULL, y = expression(R[t]), color = NULL, 
         title = "Oracle estimate (Cori A)") 

  # Panel 2: CoriB (Cori pairs data)
  p2 <- ggplot() +
    common_layers +
    geom_ribbon(data = pairs_res, aes(x = t_start, ymin = `Quantile.0.025(R)`, 
                                      ymax = `Quantile.0.975(R)`), fill = "gray", 
                alpha = 0.3) +
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", 
              linetype = "dashed") +
    geom_line(data = pairs_res, aes(x = t_start, y = `Median(R)`), 
                                    color = "#006C67") +
    labs(x = NULL, y = expression(R[t]), color = NULL, 
         title = "Cori pairs data (Cori B)") 

  # Panel 3: CoriC (Parametric Cori)
  p3 <- ggplot() +
    common_layers +
    geom_ribbon(data = coriC, aes(x = time, ymin = lowerRt, ymax = upperRt), 
                fill = "gray", alpha = 0.3) +
    geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", 
              linetype = "dashed") +
    geom_line(data = coriC, aes(x = time, y = estimate), color = "#006C67") +
    labs(x = "Time (days)", y = expression(R[t]), color = NULL, 
         title = "Parametric Cori (Cori C)") 
  
  # Combine panels
  (p1 | p2 | p3) + 
    plot_layout(axes = "collect", guides = "collect") & 
    theme(legend.position = "bottom")
}
