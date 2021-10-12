# compare the ability to reconstruct variability in PDO and NPGO
# using AR(1) models forced by SLP variability

library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)

# load data
d = read.csv("data/slp_sst_PCs_1948-2021.csv",
             stringsAsFactors = FALSE)

# attempt to reconstruct SST PC1
d$date = lubridate::parse_date_time(x = paste(d$year,d$month,"01"),orders="ymd",tz="Pacific")

# make a bespoke data set with lagged slp
sst <- filter(d, variable == "sst")
slp <- filter(d, variable == "slp")

dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc1[2:nrow(sst)],
                  slp = slp$pc1[1:nrow(slp)-1])

dat = dat[1:100,]
# Bayesian version - throws error
data_list = list(
  N = nrow(dat),
  obs_y = c(dat$sst),
  obs_x = c(dat$slp)
)

model <- stan_model("code/ar1_forcing_ss_simple.stan")
fit = optimizing(model, data = data_list)

fit = sampling(model,
           data = data_list,
           iter=3000,
           chains=1,
           control=list(adapt_delta=0.99, max_treedepth=20))
pars = rstan::extract(fit)
