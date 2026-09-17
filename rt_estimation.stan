// stan version of Rt estimation 
// last edited 9.17.2026

// # Functions
functions {
  array[] real logSEIR(real t, array[] real y, array[] real parms, 
  array[] real rdata, array[] int idata) {
    real beta = parms[1];
    real delta = parms[2];
    real gamma = parms[3];
    
    array[3] real dydt;
    dydt[1] = -beta * exp(y[3]); // logS
    dydt[2] = beta * exp(y[1] + y[3] - y[2]) - delta; // logE
    dydt[3] = delta * exp(y[2] - y[3]) - gamma; // logI
    return dydt;
  }
}

// Debugging: show values of SEIR at each t
 // print("At time t=", t, ": S=", y[1], " E=", y[2], " I=", y[3], " R=", y[4]);
//  print("dydt: ", dydt);


// # Data Module
data {
  int<lower=0> k; // number of infected individuals
  int<lower=0> r; // number of right censored individuals
  int<lower=0> l; // number of left censored individuals
  int<lower=0> start; // start observation time
  int<lower=0> end; // end observation time
  real<lower=0> t0; // start of epidemic
  array[k+1] real<lower=0> Etime; // times of exposure
  array[k+l] real<lower=0> Itime; // times of infection
  array[k+l] real<lower=0> Rtime; // times of recovery
  real<lower=0> Tmax; // final observation time
  array[k+1] int<lower=0, upper=1> Estat; // t1 event indicator (1 if observed, 0 if censored)
  array[k+l] int<lower=0, upper=1> Istat; // t2 event indicator (1 if observed, 0 if censored)
  array[k+l] int<lower=0, upper=1> Rstat; // t3 event indicator (1 if observed, 0 if censored)
  array[k+l] int<lower=0, upper=1> Erisk; 
  array[k+l] int<lower=0, upper=1> Irisk;
  array[k+l] int<lower=0, upper=1> Rrisk;
  
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
  array[0] real x_r;
  array[0] int x_i;
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
  real gamma = exp(gamma);
  
  real rhoS = rho[1];
  real rhoE = rho[2];
  real rhoI = rho[3];
  real rhoR = rho[4];
 
  array[3] real init;
  init[1] = log(rhoS);
  init[2] = log(rhoE);
  init[3] = log(rhoI);
}

// # Model Defining Module
model {
  array[k+1,3] real temp; // SEIR time series
  array[3] real parms;  // parameters 
  array[3] real init;   // initial conditions 
  
  parms[1] = beta;
  parms[2] = delta;
  parms[3] = gamma;
  
  // Debugging: show initial conditions 
  print("Initial conditions: S=", exp(init[1]), "E=", exp(init[2]), "I=", exp(init[3]), "R=", exp(init[4]));
  
  // Likelihood function
  temp = ode_rk45(logSEIR, init, t0, t1, parms, x_r, x_i, 1.0E-4, 1.0E-4, 1.06); 
  
  
  for (i in 1:size(t1) - 1) {
    real log_lik =
    event1[i] * log(beta * temp[i,1]  * temp[i,3]) - delta * (fmin(t2[i], end) - fmax(t1[i], start)) +
    event2[i] * log(delta) - gamma * (fmin(t3[i], end) - fmax(t2[i], start)) +
    event3[i] * log(gamma);
  
  target+= log_lik;
    
  // Debugging: show i, S, I, and log likeilhood
  print("Iteration: ", i);
  print("event1: ", event1[i], " temp[i,1]: ", temp[i,1], " temp[i,3]: ", temp[i,3]);
  print("log_lik: ", log_lik);
  }
  
  target += r * log(temp[k+1, 1]);
  
  target += l * log(1 - temp[start, 1]);
  
  if (use_empEIR == 1) {
    target += normal_lpdf(rhoE | emp_Esurv, emp_Esurv_se);
    target += normal_lpdf(rhoI | emp_Isurv, emp_Isurv_se);
    target += normal_lpdf(rhoR | emp_Rsurv, emp_Rsurv_se);
  }
  
  // Priors
  //target += normal_lpdf(sigma | m_sigma / T_sigma, sqrt(m_sigma) / T_sigma);
  //target += normal_lpdf(gamma | m_gamma / T_gamma, sqrt(m_gamma) / T_gamma);
  
  // target += -l * log(1-temp[k+1, 1]); 
}
