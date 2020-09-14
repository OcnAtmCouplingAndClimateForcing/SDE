library(rstan)
library(dplyr)
library(ggplot2)

d = read.csv("data/example slp sst data.csv", 
  stringsAsFactors = FALSE)

# filter out NAs
d = dplyr::filter(d, !is.na(sst.anom))
# scale?
d$sst.anom = scale(d$sst.anom)
d$slp.anom.lag2 = scale(d$slp.anom.lag2)

# d = d[1:30,] # subset first 30 observations, just for speed / illustration

d <- d %>%
  filter(year <= 1988) # restricting to 'era1'

d$t = seq(1,nrow(d))
# Fit continuous model
data_list = list(
  N = nrow(d),
  M = 10, # resolution, should be higher than 10 (more like 100 or 1000)
  obs_y = c(d$sst.anom)
)

fit_1 = stan("scripts/ar1_euler_00.stan", 
  data = data_list, pars = c("sigma","obs_sigma","pred_y","tau","x00"),
  iter=3000, chains=3, control=list(adapt_delta=0.99, max_treedepth=20))
pars = rstan::extract(fit_1)

d$pred = apply(pars$pred_y,2,mean)
d$low = apply(pars$pred_y,2,quantile,0.025)
d$hi = apply(pars$pred_y,2,quantile,0.975)

ggplot(d, aes(t, pred)) + geom_ribbon(aes(ymin=low,ymax=hi),fill="blue",col=NA,alpha=0.1) + 
  geom_line(col="dark blue",linetype = "dashed") + 
  geom_point(aes(t,sst.anom),col="blue",size=1)

ggplot(d, aes(pred, sst.anom)) + geom_point(col="red",size=1) # fit is apparently perfect!

cor(d$pred, d$sst.anom)

ggplot(d, aes(t, sst.anom - pred)) + geom_point(col="red",size=1) + geom_smooth()

# try with 1989-2013!
d2 = read.csv("data/example slp sst data.csv", 
             stringsAsFactors = FALSE)

# filter out NAs
d2 = dplyr::filter(d2, !is.na(sst.anom))
# scale?
d2$sst.anom = scale(d2$sst.anom)
d2$slp.anom.lag2 = scale(d2$slp.anom.lag2)

# d2 = d2[1:30,] # subset first 30 observations, just for speed / illustration

d2 <- d2 %>%
  filter(year %in% 1989:2013) # restricting to 'era2'

d2$t = seq(1,nrow(d2))
# Fit continuous model
data_list = list(
  N = nrow(d2),
  M = 10, # resolution, should be higher than 10 (more like 100 or 1000)
  obs_y = c(d2$sst.anom)
)

fit_2 = stan("scripts/ar1_euler_00.stan", 
             data = data_list, pars = c("sigma","obs_sigma","pred_y","tau","x00"),
             iter=3000, chains=3, control=list(adapt_delta=0.99, max_treedepth=20))
pars2 = rstan::extract(fit_2)

d2$pred = apply(pars2$pred_y,2,mean)
d2$low = apply(pars2$pred_y,2,quantile,0.025)
d2$hi = apply(pars2$pred_y,2,quantile,0.975)

ggplot(d2, aes(t, pred)) + geom_ribbon(aes(ymin=low,ymax=hi),fill="blue",col=NA,alpha=0.1) + 
  geom_line(col="dark blue",linetype = "dashed") + 
  geom_point(aes(t,sst.anom),col="blue",size=1)

ggplot(d2, aes(pred, sst.anom)) + geom_point(col="red",size=1) # fit is apparently perfect!

cor(d2$pred, d2$sst.anom)

ggplot(d2, aes(t, sst.anom - pred)) + geom_point(col="red",size=1) + geom_smooth()

