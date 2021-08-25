data {
  int N;                        // nb. of data points
  real obs_y[N];              // observed data y
  real obs_x[N];              // observed data x, not used
  real obs_sigma; // parameter for observation error sigma
}
parameters {
  //real<lower=0.0> obs_sigma;  // obs error
  real<lower=0.0,upper=1.0> gamma;        // theta = 1/tau = ocean
}
transformed parameters {
  real pred_y[N]; // array rather than vector so we can calculate sd
  real sd_y;
  pred_y[1] = 0;
  for (t in 2:(N)) {
    pred_y[t] = (1-gamma)*pred_y[t-1] + obs_x[t-1]; // prediction stored for discrete time step
  }

  // correct for mid-point -- this is just to make sure results match R version
  // for(t in 1:(N-1)) {
  //   pred_y[t] = (pred_y[t] + pred_y[t+1])/2.0;
  // }
  // // should be normalized already, but force it
  // sd_y = sd(pred_y[1:(N-1)]);
  // for(t in 1:N) {
  //   pred_y[t] = pred_y[t]/sd_y;
  // }
}
model {
  // priors
  gamma ~ normal(1,0.1); // 1/theta
  //obs_sigma ~ student_t(3,0,0.1); // variability of sn
  for(t in 1:N) {
    obs_y[t] ~ normal(pred_y[t], obs_sigma);
  }
}
