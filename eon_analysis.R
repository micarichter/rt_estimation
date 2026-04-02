# Estimate Rt using data generated from Epidemics on Networks package in Python

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

n_samp <- 2000

eon_dsa <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/eon_dsa_deg10_r04a.csv")

# fix initial cases: etimes = 0
eon_dsa$etime[eon_dsa$itime == 0] <- 0

# everyone else with etime = NA is censored; they stayed in S forever
# these people should have etimes, itimes, rtimes of inf
eon_dsa$etime <- ifelse(is.na(eon_dsa$etime), Inf, eon_dsa$etime)
eon_dsa$itime <- ifelse(eon_dsa$etime == Inf, Inf, eon_dsa$itime)
eon_dsa$rtime <- ifelse(eon_dsa$etime == Inf, Inf, eon_dsa$rtime)

# add indicators
start <- 0
end <- 100
eon_dsa$estat <- ifelse(eon_dsa$etime < end, 1, 0)
eon_dsa$istat <- ifelse(eon_dsa$itime < end, 1, 0)
eon_dsa$rtsat <- ifelse(eon_dsa$rtime < end, 1, 0)

names(eon_dsa) <- c("X", "id", "eTime", "iTime", "rTime", "estat", "istat", "rstat")

# take a sample of 2000 (or read in previously selected sample)
eon_sample <- eon_dsa[sample(nrow(eon_dsa), n_samp), ]
#write.csv(eon_sample, "sampdeg10_r04a.csv")
#eon_sample <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/sampdeg10_r04a.csv")

# function to estimate Rt with network data
eon_est <- function(dat, begin, end, width, step, obs_end) {

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
                       "rt_var" = rep(NA, length(tstarts)))
  #inits <- c(0, 0, 0, 0, 0, 0)
  
  for (i in 1:length(tstarts)) {
    try({
      tstart <- tstarts[i]
      print(tstart + width)
      #print(inits)
      dat <- HSsubset(full_dat, tstart = tstart, tstop = tstart + width)
      #DSAest <- DSAmle(dat, init = inits, method = "L-BFGS-B")
      DSAest <- DSAmle(dat, method = "L-BFGS-B")
      pvec <- as.numeric(exp(DSAest$point$point))
      beta <- pvec[1]
      delta <- pvec[2]
      gamma <- pvec[3]
      xrhoE <- pvec[4]
      xrhoI <- pvec[5]
      xrhoR <- pvec[6]
      rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
      rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
      rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
      #inits <- c(log(beta), log(delta), log(gamma), log(rhoE), log(rhoI), log(rhoR))
      
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
      Rt_out[i, "estimate"] <- beta / gamma * S_est
 
      mlesamp <- DSApred_mlesamp(DSAest)
      cis <- DSApred_ci(mlesamp, times = seq(tstart, tstart + width, 0.1))
    
      Rt_out[i, c("upperRt", "lowerRt", "rt_var")] <- 
        last(cis$bounds[, c("upperRt", "lowerRt", "Rt_var")])
    })
  }
    Rt_out = Rt_out
}


system.time(res4 <- eon_est(dat = eon_sample, begin = 0, end = 50, width = 4, 
                           step = 1, obs_end = TRUE))
#plot(res$time, res$estimate, type = "l", ylim = c(0, 5))
#write.csv(res4, "dsanet_r04_w4_10k.csv")

# r02 res
# res <- read.csv("dsanet_r02_w4_10k.csv")
# res1 <- read.csv("dsanet_r04_w4_10k.csv")

# plot(res$time, res$estimate, type = "l", ylim = c(0, 3))
# grid()

# true Rt data, (R0 = 4)
true_sim <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/seir_deg10_r04a.csv")
#write.csv(true_sim, "seir_dat_r02.csv")
#plot(true_sim$time, true_sim$true_rt, type = "l", col = "red")
#res$true_rt <- true_sim$true_rt

# true Rt data (R0 = 2)
true_sim2 <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/seir_deg10_r02.csv")
#true_sim2 <- read.csv("seir_dat_r02.csv")
#plot(true_sim2$time, true_sim2$true_rt, type = "l")

# Cori estimates
cori_gt <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/cori_pairs_deg10_r04a.csv")
#pairs_sample <- cori_gt[sample(nrow(cori_gt), 50), ]
#pairs_sample <- cori_gt[c(1:1000), ]

# find point at which S = 0.9
c_cutoff <- min(true_sim$time[true_sim$S/500000 == 0.9])
pairs_sample <- cori_gt[cori_gt$etime_infectee < c_cutoff, ]
pairs_sample <- pairs_sample[sample(nrow(pairs_sample), 2000), ]

gt_df <- data.frame("EL" = floor(pairs_sample$etime_infector),
                    "ER" = ceiling(pairs_sample$etime_infector),
                    "SL" = floor(pairs_sample$etime_infectee),
                    "SR" = ceiling(pairs_sample$etime_infectee))
gt_df[] <- lapply(gt_df, as.integer)

cori_incidence <- read.csv("/Users/micaelarichter/Library/CloudStorage/OneDrive-TheOhioStateUniversity/python/incidence_deg10_r04a.csv")
names(cori_incidence) <- c("time", "incidence")

mcmc_control <- make_mcmc_control(burnin = 1000, thin = 10, seed = 13)

config <- make_config(incid = cori_incidence$incidence,
                      method = "si_from_data",
                      si_parametric_distr = "G",
                      mcmc_control = mcmc_control,
                      n1 = 500,
                      n2 = 50,
                      seed = 918)

system.time(cori_res <- estimate_R(incid = cori_incidence$incidence, method = "si_from_data",
                       si_data = gt_df, 
                       config = config))
#write.csv(cori_res$R, "1573c_r02_deg50.csv")
#write.csv(cori_res$R, "2k_r04_deg10a_cori.csv")
#write.csv(cori_res$R, "3k_r04_cori.csv")
#plot(cori_res)

#cori_res <- read.csv("2k_r04_deg10a_cori.csv")

# plotting function
eon_plots <- function(DSA, Cori, truth, R0, width, pop, ymax) {
 dsa_dat <- DSA
 cori_dat <- Cori$R
 true_dat <- truth
 
 annot <- data.frame(
   x = end - 40,
   y = ymax - c(0.15, 0.2, 0.25) * (ymax - 0),
   label = c(
     paste0("R0 = ", R0),
     paste0("DSA window = ", width, " days"),
     paste0("sample = ", pop, " human sensors")
   )
 )
 
 dsa_dat %>%
   ggplot() +
   theme_minimal() +
   
   # Cori CI
   geom_ribbon(data = cori_dat, aes(x = (t_start + t_end) / 2, ymin = `Quantile.0.025(R)`, 
                                    ymax = `Quantile.0.975(R)`), fill = "gray", alpha = 0.3) +
   
   # DSA CI
   geom_ribbon(aes(x = time, ymin = lowerRt, ymax = upperRt), fill = "#FFDBB5", 
                alpha = 0.4) +
   
   # true Rt
   geom_line(data = true_dat, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
   
   # Cori estimate
   #geom_smooth(data = cori_dat, aes(x = t_start, y = `Median(R)`, color = "Cori estimate"),
   #             se = FALSE) +
   geom_line(data = cori_dat, aes(x = (t_start + t_end) / 2, y = `Median(R)`, color = "Cori estimate")) +
   geom_vline(xintercept = c_cutoff, linetype = "dotted") +
   # DSA estimate 
   #geom_smooth(aes(x = time, y = estimate, color = "DSA estimate"), se = FALSE) +
   geom_line(aes(x = time, y = estimate, color = "DSA estimate")) +
   
   labs(x = "Time (days)", y = expression(R[t]), color = NULL) +
   scale_color_manual(
     values = c(
       "DSA estimate" = "#DA8210",
       "Cori estimate" = "#006C67"
     )
   ) +
   geom_text(data = annot, aes(x = x, y = y, label = label), hjust = 0) +
   coord_cartesian(ylim = c(0, ymax)) 
 }

eon_plots(DSA = res4, Cori = cori_res, truth = true_sim, width = 4, R0 = 4,
          pop = n_samp, ymax = 5)

# smoothed estimate using many windows
adaptive_smooth1 <- function(windows_vec = c(4, 6)) {
  
  # estimate Rt at different window sizes
  window_list <- lapply(windows_vec, function(x) {
    eon_est(dat = eon_sample, begin = 0, end = end, width = x, 
            step = 1, obs_end = TRUE) %>%
    dplyr::select(time, estimate, beta, gamma, rt_var) %>%
    dplyr::rename(!!paste0("rt", x) := estimate,
                  !!paste0("beta", x) := beta,
                  !!paste0("gamma", x) := gamma,
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
  
  precision <- rowSums(1/var_mat)
  rt_comb <- rowSums(rt_mat/var_mat) / precision
  var_comb <- 1/precision
  
  # CIs
  ntimes <- nrow(rt_trim)
  lower <- numeric(ntimes)
  upper <- numeric(ntimes)
  
  for(i in 1:ntimes){
    
    mu <- as.numeric(rt_trim[i, rt_cols])
    vars <- as.numeric(rt_trim[i, var_cols])
    Sigma <- diag(vars)
    
    samples <- mvrnorm(4000, mu = mu, Sigma = Sigma)
    
    weights <- 1/vars
    weights <- weights/sum(weights)
    rt_sim <- apply(samples, 1, function(x)  sample(x, size = 1, prob = weights))
    
    lower[i] <- quantile(rt_sim, 0.025)
    upper[i] <- quantile(rt_sim, 0.975)
  }
  
  res_df <- 
    data.frame(
    time = rt_trim$time,
    estimate = rt_comb,
    lower = lower,
    upper = upper
  )
  return(list(rt_trim, res_df))
}

windows_deg10_4a <- adaptive_smooth1()
write.csv(windows_deg100_2, "smooth_r04_deg100.csv")  


# adaptive window plot
adaptive_smooth_plot <- function(DSAsmooth, Cori, truth, R0, pop, ymax) {
  
  dsa_dat <- DSAsmooth
  #cori_dat <- Cori
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
    coord_cartesian(ylim = c(0, ymax)) 
}

adaptive_smooth_plot(DSA = windows_deg10_4a[[2]], Cori = cori_res, truth = true_sim, R0 = 4,
          pop = n_samp, ymax = 5)
