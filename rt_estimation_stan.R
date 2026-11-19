# MLE in stan to improve behavior & speed 
# script to accompany rt_window.stan 
# last edited 11.5.2026

library(cmdstanr)
set.seed(918)

# number of human sensors 
n_sens <- 2000

# read in true data
true_sim <- read.csv("eon_simulations/true_deg10_r02.csv")

# input and format data
eon_dsa <- read.csv("eon_simulations/dat_deg10_r02.csv")

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

  Rt_out <- data.frame(
    time = numeric(n_times),
    estimate = numeric(n_times),
    lowerRt = numeric(n_times),
    upperRt = numeric(n_times),
    beta = numeric(n_times),
    gamma = numeric(n_times),
    S_est = numeric(n_times),
    rt_var = numeric(n_times),
    prop_na_rt = numeric(n_times),
    prop_na_S = numeric(n_times)
  )
  
  Rt_draws <- vector("list", n_times)
  
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
      n_sens = nrow(dat_i),
      t0 = t0,
      ltimes = length(times),
      times = times,
      Etime = dat_i$Etime, Itime = dat_i$Itime, Rtime = dat_i$Rtime,
      Estat = dat_i$Estat, Istat = dat_i$Istat, Rstat = dat_i$Rstat,
      Erisk = dat_i$Erisk, Irisk = dat_i$Irisk, Rrisk = dat_i$Rrisk,
      Etime_idx = Etime_idx,
      use_empEIR = 1
    )
    
    # try with mcmc
    fit_draws <- tryCatch({
      mod$sample(
        data = stan_data,
        seed = 123,
        chains = 4,              # 4 chains is standard for robust HMC diagnostics
        parallel_chains = 4,     # Run in parallel if your machine has 4+ cores
        iter_warmup = 1000,      # Adjust these if you need more/fewer iterations
        iter_sampling = 1000,
        adapt_delta = 0.95,      # Helps prevent divergent transitions in early windows
        refresh = 0
      )
    }, error = function(e) {
      cat("HMC sampling failed for time window:", tstart, "\n")
      return(NULL)
    })
#Stan optimizer + Laplace approximation
#     fit_draws <- tryCatch({
#       # step a: Fast MLE Optimization
#       fit_mle <- mod$optimize(
#         data = stan_data,
#         seed = 123,
#         jacobian = TRUE,
#         algorithm = "lbfgs",
#         refresh = 0 # suppresses terminal output
#       )
#       # step b: Laplace approx
#       lap <- mod$laplace(
#         data = stan_data,
#         mode = fit_mle,
#         draws = 4000, # Number of draws to estimate variance/CIs
#         refresh = 0
#       )
#       dummy_check <- lap$draws("Rt_end")
#       lap # return laplace object if successful
#       #print("lap")
#     }, error = function(e) {
#       cat("going to pf \n")
# 
#     # step c: if laplace fails use pathfinder (variational inference)
#     tryCatch({
#       pf <- mod$pathfinder(
#         data = stan_data,
#         draws = 4000,
#         seed = 123,
#         refresh = 0
#       )
#       return(pf)
#     }, error = function(e) {
#     return(NULL)
#   })
# })
    if (!is.null(fit_draws)) {
      # laplace object supports the exact same $summary() calls as MCMC
      # extract raw draws first to check for nas
      raw_rt <- as.numeric(fit_draws$draws("Rt_end"))
      raw_S <- as.numeric(fit_draws$draws("S_end")) 
      
      # prop of draws thrown out
      prop_na_rt <- mean(!is.finite(raw_rt))
      prop_na_S <- mean(!is.finite(raw_S))
      
      sum_rt <- fit_draws$summary("Rt_end", .args = list(na.rm = TRUE))
      sum_beta <- fit_draws$summary("beta")
      sum_gamma <- fit_draws$summary("gamma")
      sum_S <- fit_draws$summary("S_end", .args = list(na.rm = TRUE))
      
      #if (is.na(sum_S$median)) browser()
      
      # calculate variance
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
        rt_var = rt_var,
        prop_na_rt = prop_na_rt,
        prop_na_S = prop_na_S
      )
      
      Rt_draws[[i]] <- raw_rt
    } else {
      Rt_out[i, ] <- c(tstart + width, rep(NA, 9))
      Rt_draws[[i]] <- rep(NA, 4000)
    }
  }
  names(Rt_draws) <- as.character(Rt_out$time)
  return(list(summary = Rt_out, draws = Rt_draws))
}

hmc8 <- eon_est_stan(eon_sample, 0, 150, 8, 1, EIRcumhaz = EIRcumhaz)

plot(res1$time, res1$estimate, type = "l", ylim = c(0, 2))
lines(true_sim$time, true_sim$true_rt, col = "red")
grid()

ggplot(data = res10) +
  geom_line(aes(x = time, y = estimate), color = "#DA8210") +
  geom_ribbon(aes(x = time, ymin = lowerRt, ymax = upperRt), fill = "#FFDBB5", alpha = 0.4) +
  geom_line(data = true_sim, aes(x = time, y = true_rt), color = "black", linetype = "dashed") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 150))

# combine different widths 
adaptive_smooth_stan <- function(eon_res_list, windows_vec = c(2, 4, 6, 8, 10, 12), CIs = TRUE) {
  
  windows_vec <- sort(windows_vec)
  
  # Select and rename columns for each window's output
  window_list <- lapply(seq_along(windows_vec), function(i) {
    x <- windows_vec[i]
    res <- eon_res_list[[i]]$summary |>
      dplyr::select(time, estimate, beta, gamma, rt_var) |>
      dplyr::rename(!!paste0("rt", x) := estimate,
                    !!paste0("beta", x) := beta,
                    !!paste0("gamma", x) := gamma,
                  #  !!paste0("R0", x) := R0,
                    !!paste0("var", x) := rt_var)
    return(res)
  })
  
  # housekeeping
  rt_all <- Reduce(function(x, y) {
    dplyr::left_join(x, y, by = "time", na_matches = "never")
  }, window_list)
  
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
  precision <- rowSums(1/var_mat, na.rm = TRUE)
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
      t_val <- as.character(rt_trim$time[i]) # look up draws
      
      mu <- as.numeric(rt_trim[i, rt_cols])
      vars <- as.numeric(rt_trim[i, var_cols])
      
      # keep only finite and not NA estimates 
      stable <- is.finite(mu) & is.finite(vars) & mu > 0 & vars > 0
      
      weights <- 1 / vars[stable]
      weights <- weights / sum(weights)
      
      # collect draws for only stable windows 
      stable_idx <- which(stable)
      all_draws <- vector("list", length(stable_idx))
      
      for(j in seq_along(stable_idx)) {
        w_idx <- stable_idx[j]
        all_draws[[j]] <- eon_res_list[[w_idx]]$draws[[t_val]]
      }
      
      # if there is only one stable window, take the draws directly
      if(length(stable_idx) == 1) {
        rt_sim <- all_draws[[1]]
      } else {
        # sample mixture component indices based on window weights  
        window_idx <- sample(seq_along(stable_idx), size = 4000, 
                                 replace = TRUE, prob = weights)
        rt_sim <- numeric(4000)
        # get specific number of draws needed from each window's posterior pool
        for (j in seq_along(stable_idx)) {
          mask <- window_idx == j
          n_j <- sum(mask)
          if (n_j > 0) {
            # sample from the specific window's stan draws
            rt_sim[mask] <- sample(all_draws[[j]], size = n_j, replace = TRUE)
          }
        }
      }
      # no exp needed since stan draws on the natural scale
      lower[i] <- quantile(rt_sim, 0.025, na.rm = TRUE)
      upper[i] <- quantile(rt_sim, 0.975, na.rm = TRUE)
    }
    res_df$lower <- lower
    res_df$upper <- upper
  }
  
  res_df <- subset(res_df, time >= max(windows_vec))
  
  return(list(rt_trim, rt_all, res_df))
}

window_list <- c(2, 4, 6, 8, 10, 12)
final_results <- adaptive_smooth_stan(
  eon_res_list = list(res2, res4, res6, res8, res10, res12), 
  windows_vec = window_list, 
  CIs = TRUE
)

adaptive_smooth_plot1(DSA = final_results[[3]], Cori = cori_ora,
                      truth = true_sim, ymax = 3)

