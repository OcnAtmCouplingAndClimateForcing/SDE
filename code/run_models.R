library(rstan)
library(dplyr)
library(ggplot2)

d = read.csv("example slp sst data.csv", 
  stringsAsFactors = FALSE)

# filter out NAs
d = dplyr::filter(d, !is.na(sst.anom))
# scale?
d$sst.anom = scale(d$sst.anom)
d$slp.anom.lag2 = scale(d$slp.anom.lag2)

# fit first era 
d = dplyr::filter(d, d$era == unique(d$era)[1])

# Fit continuous model
data_list = list(
  N = nrow(d),
  y = c(d$sst.anom),
  x = c(d$slp.anom.lag2),
  ts = rep(1, nrow(d)), # monthly time step
  a = 1/3 # monthly time step
)

fit_1 = stan("continuous_toy.stan", 
  data = data_list, pars = c("pred_y","sigma_resid"),
  iter=3000, chains=4)
pars = rstan::extract(fit_1)

d$pred = apply(pars$pred_y,2,mean)
d$low = apply(pars$pred_y,2,quantile,0.025)
d$hi = apply(pars$pred_y,2,quantile,0.975)


