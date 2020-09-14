data {
  int<lower=1> N; // number of data points
  real y[N]; // data, SST in this case
  real ts[N]; // time steps to evaluate at
  real x[N]; // SLP data for this example
  real a;
}
parameters {
  //real<lower=0> a;
  real<lower=0> sigma_resid;
}
transformed parameters {
  vector[N] pred_y;
  vector[N] sigma;
  vector[N] b; // derived, as theta/a or x[t]/a
  real cnst;
  cnst = sigma_resid*sigma_resid/(2*a);
  pred_y[1] = 0;
  for(t in 2:N) {
    // dy(t)/dt = a*(b-y(t)) = ab - ay(t) = a*x[t]/a - ay(t) = x(t) - ay(t)
    b[t] = x[t] / a;
    // y[t-1] equivalent to r0
    pred_y[t] = y[t-1]*exp(-a*ts[t]) + b[t]*(1-exp(-a*ts[t]));
    sigma[t] = sqrt(cnst*(1 - exp(-2*a*ts[t])));
  }
}
model {
  sigma_resid ~ student_t(3,0,2);
  a ~ student_t(3,0,2);
  for (t in 2:N) {
    y[t] ~ normal(pred_y[t], sigma[t]);
  }
}
