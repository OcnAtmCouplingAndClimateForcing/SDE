// ar(1) constinuous stochastic differential equation
// code modified from discussion here,
// https://discourse.mc-stan.org/t/ito-process-as-numerical-solution-of-stochastic-differential-equation/9192/25
data {
  int N;                        // nb. of data points
  //int M;                        // nb. of auxiliary data between each pair of data points
  vector[N] obs_y;              // observed data y
  //vector[N] obs_x;              // observed data x, not used
}
transformed data {
  //real dt = 1/(M+0.0);            // fixed time step
}
parameters {
  real<lower=0.0> sigma;      // data simulated from given sigma
  //real<lower=0.0> obs_sigma;  // obs error
  real<lower=0.0> tau;        // theta = 1/tau = ocean 
  //real x00;                   // initial condition
  //real sn[M*N];               // Wiener process' std normal steps
}
transformed parameters {
  //real sigma_sqrt_dt;
  //real x0;
  real theta;
  real sigma_analytic;
  vector[N] pred_y;
  //vector[M] pred_interval; // dimensioned by M, where M represents end of time step
  theta = 1/tau;
  sigma_analytic = sqrt((sigma^2)*(1-exp(-2*theta))/(2*theta));
  pred_y[1]=0;
  for(i in 2:N) {
    pred_y[i] = obs_y[i-1]*exp(-theta);
  }
}
model {
  // priors
  tau ~ normal(3,3); // 1/theta
  sigma ~ normal(0,1); // variability of sn
  for(i in 2:N) {
    obs_y[i] ~ normal(pred_y[i], sigma_analytic); 
  }

}
