// stan version of Rt estimation 
// last edited 10.8.2026

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

// Debugging: show values of SEIR at each t
 // print("At time t=", t, ": S=", exp(y[1]), " E=", exp(y[2]), " I=", exp(y[3]), " R=", exp(y[4]));
//  print("dydt: ", dydt);


// # Data Module
data {
  int<lower=0> N; // total number of human sensors 
  int<lower=0> start; // start observation time
  int<lower=0> end; // end observation time
  int<lower=0> width;
  real<lower=0> t0; // start of epidemic
  array[N] real<lower=0> Etime; // times of exposure
  array[N] real<lower=0> Itime; // times of infection
  array[N] real<lower=0> Rtime; // times of recovery
  array[N] int<lower=0, upper=1> Estat; // t1 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Istat; // t2 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Rstat; // t3 event indicator (1 if observed, 0 if censored)
  array[N] int<lower=0, upper=1> Erisk; 
  array[N] int<lower=0, upper=1> Irisk;
  array[N] int<lower=0, upper=1> Rrisk;
  array[N] int<lower=1> Etime_idx;
  int<lower=1> Tmax; // final observation time
  array[Tmax - width + 1] real times;
  int<lower=1> ltimes;

  int<lower=0, upper=1> use_empEIR;
  array[ltimes] real emp_Esurv;
  array[ltimes] real emp_Isurv;
  array[ltimes] real emp_Rsurv;
  array[ltimes] real emp_Esurv_se;
  array[ltimes] real emp_Isurv_se;
  array[ltimes] real emp_Rsurv_se;
}

// Transformed Data
transformed data {
  array[0] real x_r;
  array[0] int x_i;
  
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
 
  vector[3] init;
  init[1] = log(rhoS);
  init[2] = log(rhoE);
  init[3] = log(rhoI);
}

// # Model Defining Module
model {
  array[3] real parms;  // parameters 
  parms[1] = beta;
  parms[2] = delta;
  parms[3] = gamma;
  
  // Debugging: show initial conditions 
  //print("Initial conditions: S=", exp(init[1]), "E=", exp(init[2]), "I=", exp(init[3]), "R=", exp(init[4]));
  
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
    target += normal_lpdf(rhoE | emp_Esurv, emp_Esurv_se);
    target += normal_lpdf(rhoI | emp_Isurv, emp_Isurv_se);
    target += normal_lpdf(rhoR | emp_Rsurv, emp_Rsurv_se);
  }
}
