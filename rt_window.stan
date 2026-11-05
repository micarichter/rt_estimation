// stan version of Rt estimation 
// last edited 11.5.2026

// # Functions
functions {
  vector logSEIR(real t,                 // time
                vector y,                // state
                array[] real parms) {    // parameters
    real beta = parms[1];
    real delta = parms[2];
    real gamma = parms[3];
    vector[3] dydt;
    dydt[1] = -beta * exp(y[3]); // logS
    dydt[2] = beta * exp(y[1] + y[3] - y[2]) - delta; // logE
    dydt[3] = delta * exp(y[2] - y[3]) - gamma; // logI
    return dydt;
  }
}

// # Data Module
data {
  int<lower=0> N; // total number of human sensors 
  real t0; // start of epidemic
  array[N] real<lower=0> Etime; // times of exposure
  array[N] real<lower=0> Itime; // times of infection
  array[N] real<lower=0> Rtime; // times of recovery
  array[N] int<lower=0, upper=1> Estat; // t1 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Istat; // t2 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Rstat; // t3 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Erisk; 
  array[N] int<lower=0, upper=1> Irisk;
  array[N] int<lower=0, upper=1> Rrisk;
  
  int<lower=1> ltimes;
  array[ltimes] real times; 
  array[N] int<lower=1> Etime_idx;

  int<lower=0, upper=1> use_empEIR;
  real emp_Esurv;
  real emp_Isurv;
  real emp_Rsurv;
  real emp_Esurv_se;
  real emp_Isurv_se;
  real emp_Rsurv_se;
}

// Transformed Data
transformed data {
  real infected = 0;
  real time_in_E = 0;
  real recovered = 0;
  real time_in_I = 0;
  
  // Initial state indicator counts
  int<lower=0> nE_init = 0; // (Erisk == 0) & (Irisk == 1)
  int<lower=0> nI_init = 0; // (Irisk == 0) & (Rrisk == 1)
  int<lower=0> nR_init = 0; // (Rrisk == 0) & (Rstat == 1)
  
  for (i in 1:N) {
    // Initial state counts
    if (Erisk[i] == 0 && Irisk[i] == 1) {
      nE_init += 1;
    }
    if (Irisk[i] == 0 && Rrisk[i] == 1) {
      nI_init += 1;
    }
    if (Rrisk[i] == 0 && Rstat[i] == 1) {
      nR_init += 1;
    }
    // Likelihood transition counts
    if (Irisk[i] == 1 && Istat[i] == 1) {
      infected += 1.0;
    }
    time_in_E += (Itime[i] - Etime[i]);
    if (Rrisk[i] == 1 && Rstat[i] == 1){
      recovered += 1.0;
    }
    time_in_I += (Rtime[i] - Itime[i]);
  }
}

// # Parameter Module
parameters {
  real lnbeta;
  real lndelta;
  real lngamma;
  simplex[4] rho;
}

transformed parameters {
  real beta = exp(lnbeta);
  real delta = exp(lndelta);
  real gamma = exp(lngamma);
  
  real rhoS = rho[1];
  real rhoE = rho[2];
  real rhoI = rho[3];
  real rhoR = rho[4];
 
  vector[3] init; // to avoid log(0)
  init[1] = log(fmax(rhoS, 1e-12));
  init[2] = log(fmax(rhoE, 1e-12));
  init[3] = log(fmax(rhoI, 1e-12));
  
  array[3] real parms;  // parameters 
  parms[1] = beta;
  parms[2] = delta;
  parms[3] = gamma;
  
}

// # Model Defining Module
model {
  // Likelihood function
  // Solve ode to get logS and logI
  array[ltimes] vector[3] temp = ode_rk45(logSEIR, init, t0, times, parms);
  
  for (i in 1:N) {
    int idx = Etime_idx[i];
    // from E
    if (Erisk[i] == 1) {
      target += temp[idx, 1]; // logS(Etime[Erisk==1])
    }
    if (Erisk[i] == 1 && Estat[i] == 1) {
      target += lnbeta + temp[idx, 3]; // lnbeta + logI(Etime)
    }
  }
  
  if (nE_init > 0) {
    target += nE_init * log(rhoE); // initially exposed
  }
  // from I
  target += lndelta * (infected + 0.5) - exp(lndelta) * time_in_E;
  if (nI_init > 0) {
    target += nI_init * log(rhoI); // initially infected
  }
  // from R
  target += lngamma * (recovered + 0.5) - exp(lngamma) * time_in_I;
  if (nR_init > 0) {
    target += nR_init * log(rhoR); // initially recovered
  }
  // From empirical cumulative hazards 
  if (use_empEIR == 1) {
    target += normal_lpdf(rhoE | emp_Esurv, fmax(emp_Esurv_se, 1e-4));
    target += normal_lpdf(rhoI | emp_Isurv, fmax(emp_Isurv_se, 1e-4));
    target += normal_lpdf(rhoR | emp_Rsurv, fmax(emp_Rsurv_se, 1e-4));
  }
}

// # Generated quantities  
generated quantities {
  array[ltimes] vector[3] state;
  real Rt_end;
  real S_end;
  
  state = ode_rk45(logSEIR, init, t0, times, parms);
  
  S_end = exp(state[ltimes, 1]);
  Rt_end = exp(lnbeta - lngamma + state[ltimes, 1]);
}

