// ar(1) constinuous stochastic differential equation
// code modified from discussion here,
// https://discourse.mc-stan.org/t/ito-process-as-numerical-solution-of-stochastic-differential-equation/9192/25
data {
  int N;                        // nb. of data points
  int M;                        // nb. of auxiliary data between each pair of data points
  vector[N] obs_y;              // observed data y
  vector[N] obs_x;              // observed data x, not used
}
transformed data {
  real dt = 1/(M+0.0);            // fixed time step
}
parameters {
  //real<lower=0.0> sigma;      // data simulated from given sigma
  real<lower=0.0> obs_sigma;  // obs error
  real<lower=0.0,upper=1.0> gamma;        // theta = 1/tau = ocean
  real x00;                   // initial condition
  //real sn[M*N];               // Wiener process' std normal steps
}
transformed parameters {
  //real sigma_sqrt_dt;
  real x0;
  real theta;
  vector[M] pred_interval; // dimensioned by M, where M represents end of time step
  vector[N] pred_y;
  real sigma;
  real sigma_sqrt_dt;
  sigma = 1;
  sigma_sqrt_dt = sigma * sqrt(dt);
  x0 = x00;

  for (i in 1:N) {
    // nested loop equivalent to SDE integrator
    // that outputs pred y, given white noise sn, time step dt, and parameter controlling variability sigma
    pred_interval[1] = (1-gamma)*x0 + sigma_sqrt_dt * obs_x[(i - 1) * M + 1];
    pred_y[i] = pred_interval[1]; // prediction stored for discrete time step
    x0 = pred_interval[1]; // reset initial condition for next time step
  }
}
model {
  // priors
  //real sig;
  x00 ~ normal(obs_y[1],2); // initial state
  gamma ~ normal(0.5,1); // 1/theta
  //sn ~ normal(0, 1); // white noise with sd = 1
  obs_sigma ~ normal(0,0.1); // variability of sn
  //sigma ~ normal(0,1); // variability of sn
  //sig = sqrt(sigma^2*(1-exp(-2*theta*(1.0/M)))/(2*theta));
  obs_y ~ normal(pred_y, obs_sigma);
}
