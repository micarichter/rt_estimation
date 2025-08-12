# Estimate R(t) using DSA and Cori methods
rt_est <- function(n, begin, end, width, step, R0, delta, gamma, I_init, 
                   R0_int = NULL, R0_end = NULL, int_time1 = NULL, 
                   int_time2 = NULL, days1 = NULL, days2 = NULL) {
  #epidemic <- SEIRepidemic(..., tmax = end + 1)
  #full_dat <- HSdat(n, epidemic, tstart = begin, tstop = end)
  delta_true <- delta
  gamma_true <- gamma
  
  if(!is.null(R0_int)) {
    R0 <- specify_arnaught(R0_vec = c(R0, R0_int, R0_end), 
                           change_start_vec = c(int_time1, int_time2), 
                           change_end_vec = c(int_time1 + days1, int_time2 + days2),
                                              NT = end)
    seir_dat <- simulate_seir_ode(arnaught = R0, t_E = 1/delta_true, t_I = 1/gamma_true, 
                                  N = n, S_init = n - I_init, E_init = 0, 
                                  I_init = I_init, n_t = end, n_steps_per_t = 10)
  } else {
    R0 <- rep(R0, end + 1)
      seir_dat <- simulate_seir_ode(arnaught = R0, t_E = 1/delta_true, t_I = 1/gamma_true, 
                                    N = n, S_init = n - I_init, E_init = 0, 
                                    I_init = I_init, n_t = end, n_steps_per_t = 10) 
  }
  seir_dat$incidence <- pmax(0, round(seir_dat$dS))
  seir_dat$incidence[is.na(seir_dat$incidence)] <- 0

  indiv_dat <- seir_dat %>%
    uncount(incidence) %>%
    select(eTime = time)

  # wiggle the eTimes
  indiv_dat$eTime <- indiv_dat$eTime - runif(nrow(indiv_dat), min = 0, max = 1)

  # add initially infected
  init_dat <- data.frame(eTime = rep(0, length(I_init)))
  indiv_dat <- rbind(indiv_dat, init_dat)
  indiv_dat$id <- seq(1, nrow(indiv_dat))
  indiv_dat$eTime <- as.numeric(indiv_dat$eTime)

  # augment data according to rates delta and gamma
  indiv_dat$iTime <- ifelse(indiv_dat$eTime == 0, 0,
                            indiv_dat$eTime + rexp(nrow(indiv_dat), delta_true))
  indiv_dat$rTime <- indiv_dat$iTime + rexp(nrow(indiv_dat), gamma_true)

  never_inf <- data.frame(id = seq(nrow(indiv_dat) + 1, n),
                          eTime = end + 1, iTime = end + 1, rTime = end + 1)
  full_dat <- rbind(indiv_dat, never_inf)
  full_dat <- full_dat[, c("id", "eTime", "iTime", "rTime")]

  # add indicators
  full_dat$estat <- ifelse(full_dat$eTime < end, 1, 0)
  full_dat$istat <- ifelse(full_dat$iTime < end, 1, 0)
  full_dat$rstat <- ifelse(full_dat$rTime < end, 1, 0)
  full_dat <- EPIdat(full_dat, begin, end)

  tstarts <- seq(begin, end - width, by = step)
  Rt_out <- data.frame("time" = rep(NA, length(tstarts)),
                       "estimate" = rep(NA, length(tstarts)),
                       "upperRt" = rep(NA, length(tstarts)),
                       "lowerRt" = rep(NA, length(tstarts)),
                       "beta" = rep(NA, length(tstarts)),
                       "delta" = rep(NA, length(tstarts)),
                       "gamma" = rep(NA, length(tstarts)),
                       "est_S" = rep(NA, length(tstarts)),
                       "true_S" = rep(NA, length(tstarts)),
                       "true_rt" = rep(NA, length(tstarts)))

  for (i in 1:length(tstarts)) {
      try({
        tstart <- tstarts[i]
        print(tstart + width)
        dat <- HSsubset(full_dat, tstart = tstart, tstop = tstart + width)
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

        S_est <- last(SEIRepidemic(beta = beta, delta = delta, gamma = gamma,
                              rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, tmin = tstart,
                              tmax = tstart + width, tstep = 0.01)[, "S"])

        S_true <- seir_dat$S[seir_dat$time == tstart + width] / n
        # true_rt <- R0[tstart] * seir_dat$S[tstart] / (seir_dat$S[tstart] +
        #                                         seir_dat$E[tstart] +
        #                                         seir_dat$I[tstart] +
        #                                         seir_dat$R[tstart])
        
        true_rt <- R0[tstart + width] * S_true * n / (S_true * n +
                                                      seir_dat$E[tstart + width] +
                                                      seir_dat$I[tstart + width] +
                                                      seir_dat$R[tstart + width])

        Rt_out[i, c("est_S", "true_S")] <- c(S_est, S_true)
        Rt_out[i, c("beta", "delta", "gamma")] <- c(beta, delta, gamma)
        Rt_out[i, "estimate"] <- beta / gamma * S_est
        Rt_out[i, "true_rt"] <- true_rt

        mlesamp <- DSApred_mlesamp(DSAest)
        cis <- DSApred_ci(mlesamp, times = seq(tstart, tstart + width, 0.1))

        Rt_out[i, c("time", "upperRt", "lowerRt")] <- last(cis[, c("time", "upperRt", "lowerRt")])
        })
  }
    parlist <- {
      list(
        t_E = 1 / delta_true,
        t_I = 1 / gamma_true
      )
    }
    parlist$true_mean_GI = (parlist$t_E + parlist$t_I)
    parlist$true_var_GI = parlist$t_E^2 + parlist$t_I^2
    assign("parlist", parlist, envir = .GlobalEnv)
    Cori_out <- get_cori(seir_dat, icol_name = 'incidence', window = 1)
    out <- list(
      Rt_out = Rt_out,
      Cori_out = Cori_out,
      parameters = list(delta = delta_true, gamma = gamma_true, width = width, end = end)
    )
}

# Plotting function for Rt est
rt_plot <- function(out, ymax) {  
  rt_out <- out$Rt_out
  cori_out <- out$Cori_out
  delta <- out$parameters$delta
  gamma <- out$parameters$gamma
  width <- out$parameters$width
  end <- out$parameters$end
  
  annot <- data.frame(
    x = end - 50,
    y = ymax - c(0.15, 0.2, 0.25) * (ymax - 0),
    label = c(
      paste0("Mean time in E = ", 1/delta, " days"),
      paste0("Mean time in I = ", 1/gamma, " days"),
      paste0("DSA window = ", width, " days")
    )
  )
  
  rt_out %>%
    ggplot() +
    theme_minimal() +
    geom_line(aes(x = time, y = true_rt), color = "black") +
    geom_ribbon(data = cori_out, aes(x = time, ymin = Cori.025, ymax = Cori.975), 
                fill = "gray", alpha = 0.3) +
    geom_line(data = cori_out, aes(x = time, y = Cori.mean), color = "black", 
              linetype = "dotted") +
    geom_ribbon(aes(x = time, ymin = lowerRt, ymax = upperRt), fill = "#53A548", 
                alpha = 0.4) +
    geom_line(aes(x = time, y = estimate, group = 1), color = "#53A548", linetype = "dashed") +
    labs(x = "Time (days)", y = expression(R[t])) +
    geom_text(data = annot, aes(x = x, y = y, label = label), hjust = 0)
}

