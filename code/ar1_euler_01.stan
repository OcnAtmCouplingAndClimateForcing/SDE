// ar(1) constinuous stochastic differential equation
// code modified from discussion here,
// https://discourse.mc-stan.org/t/ito-process-as-numerical-solution-of-stochastic-differential-equation/9192/25
data {
  int N;                        // nb. of data points
  int M;                        // nb. of auxiliary data between each pair of data points
  vector[N] obs_y;              // observed data y
  vector[M*N] obs_x;              // observed data x, aka forcing variable
}
transformed data {
  real dt = 1/(M+0.0);            // fixed time step
}
parameters {
  real<lower=0.0> sigma;      // data simulated from given sigma
  //real<lower=0.0> obs_sigma;  // obs error
  //real<lower=0.0> theta;        // theta = 1/tau = ocean 
  real x00;                   // initial condition
  //real sn[M*N];               // Wiener process' std normal steps
}
transformed parameters {
  real sigma_sqrt_dt;
  real x0;
  real tau;
  real theta;
  vector[M] pred_interval; // dimensioned by M, where M represents end of time step
  vector[N] pred_y;
  sigma_sqrt_dt = sigma * sqrt(dt);
  theta = 0.2;
  tau = 1/theta;
  x0 = x00;
  for (i in 1:N) {
    // nested loop equivalent to SDE integrator
    // that outputs pred y, given white noise sn, time step dt, and parameter controlling variability sigma
    for (j in 1:M) {
      // wi+1 = wi + µwi∆ti + σwi∆W
      if(j==1) {
        pred_interval[j] = x0 + sigma_sqrt_dt * obs_x[(i - 1) * M + j];
      } else {
        pred_interval[j] = pred_interval[j - 1] - theta * pred_interval[j - 1]*dt + sigma_sqrt_dt * obs_x[(i - 1) * M + j];
      }
    }
    pred_y[i] = pred_interval[M]; // prediction stored for discrete time step
    x0 = pred_interval[M]; // reset initial condition for next time step
  }
}
model {
  // priors
  //x00 ~ normal(0,3); // initial state
  #theta ~ normal(0.5,1); // ocean memory
  //sn ~ normal(0, 1); // white noise
  //sigma ~ student_t(3,0,3); // variability of sn
  //obs_sigma ~ student_t(3,0,3); // obs error sigma
  
  // observation model
  obs_y ~ normal(pred_y, 0.01);
}
