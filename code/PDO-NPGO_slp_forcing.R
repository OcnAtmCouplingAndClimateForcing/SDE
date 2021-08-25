# compare the ability to reconstruct variability in PDO and NPGO
# using AR(1) models forced by SLP variability

library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)

# plot settings
theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load data
d = read.csv("data/slp_sst_PCs_1948-2021.csv",
             stringsAsFactors = FALSE)

# these data are PC1/PC2 of NE Pacific SLPa/SSTa for 1/1948 - 5/2021
# each PC time series has already been scaled

# plot to check
plot.dat <- d %>%
  pivot_longer(cols = c(-m, -month, -year, -variable)) %>%
  mutate(decimal.year = year + (m-0.5)/12)

ggplot(plot.dat, aes(decimal.year, value)) +
  geom_line() +
  facet_grid(name~variable)

# looks as we would expect - white noise for the SLP PCs, red noise for the SST PCs (especially PC1)

# find the decorrelation scale for each SST PC to estimate theta

print(acf(filter(d, variable == "sst")$pc1)) # above 0.5 for lags 1-7
print(acf(filter(d, variable == "sst")$pc2)) # above 0.5 for lags 1-3

# and find peak cross correlation for PC1s / PC2s
print(ccf(filter(d, variable == "slp")$pc1, filter(d, variable == "sst")$pc1)) # peak at lag1
print(ccf(filter(d, variable == "slp")$pc2, filter(d, variable == "sst")$pc2)) # also peak at lag1; much weaker


# attempt to reconstruct SST PC1
d$date = lubridate::parse_date_time(x = paste(d$year,d$month,"01"),orders="ymd",tz="Pacific")

# make a bespoke data set with lagged slp
sst <- filter(d, variable == "sst")
slp <- filter(d, variable == "slp")

dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc1[2:nrow(sst)],
                  slp = slp$pc1[1:nrow(slp)-1])

ar_ls = function(time,forcing,gamma) {
  #S(t+1) = (1-GAMMA*DT)*S(t) + F(t)*DT
  forcing = c(forcing - mean(forcing))
  T=length(forcing)
  sig = 0

  for(t in 1:(T-1)) {
    #sig[t+1] = -theta*sig[t] + forcing[t]
    sig[t+1] = (1-gamma)*sig[t] + forcing[t]
  }

  # next estimates are linearly de-trended
#  s.sig = sig
  sig = sig - lm(sig ~ time)$fitted.values
  # interpolate output on the original time grid
  s.sig=(sig[-1]+sig[-T])/2 # midpoint
  # final step is normalize
  s.sig=s.sig/sd(s.sig)
  return(s.sig)
}

calc_ss = function(theta) {
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, theta)
  ss = -sum((pred_ts - dat$sst[-1])^2) # return -SS for optim
}

o = optimize(f=calc_ss, interval = c(0,1))

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, gamma = o$minimum)

pred.pdo = data.frame(t = dat$date,
                      sst.pc1 = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor(pred.pdo$sst.pc1, pred.pdo$integrated.slp)

ggplot(pred.pdo, aes(integrated.slp, sst.pc1)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = F)

# this is getting into the weeds, but doesn't seem quite like a linear relationship
summary(mgcv::gam(sst.pc1 ~ s(integrated.slp, k = 4), data = pred.pdo))

pred.pdo <- pred.pdo %>%
  pivot_longer(cols = -t)

pdo.reconstruct <- ggplot(pred.pdo, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP PC1", "SST PC1")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("PDO (r = 0.37)") +
  ylab("Anomaly")

# compare with NPGO
dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc2[2:nrow(sst)],
                  slp = slp$pc2[1:nrow(slp)-1])

calc_ss = function(theta) {
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, theta)
  ss = -sum((pred_ts - dat$sst[-1])^2) # return -SS for optim
}

o = optimize(f=calc_ss, interval = c(0,1))

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, gamma = o$minimum)

pred.npgo = data.frame(t = dat$date,
                      sst.pc2 = dat$sst,
                      integrated.slp = c(0, -as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor(pred.npgo$sst.pc2, pred.npgo$integrated.slp)

ggplot(pred.npgo, aes(integrated.slp, sst.pc2)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = F)

pred.npgo <- pred.npgo %>%
  pivot_longer(cols = -t)

npgo.reconstruct <-ggplot(pred.npgo, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP PC2", "SST PC2")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("NPGO (r = 0.08)") +
  ylab("Anomaly")

png("./figs/PDO_NPGO_reconsturctions.png", width=6, height=8, units='in', res=300)
ggarrange(pdo.reconstruct, npgo.reconstruct, ncol=1)
dev.off()


# Bayesian version. M = 1 here, because the forcing variable obs_x
# is already discretized and there's no missing values
data_list = list(
  M = 1,
  N = nrow(dat),
  obs_y = c(dat$sst),
  obs_x = c(dat$slp)
)

# first Bayesian model is estimating both gamma, obs_sigma, and sigma,
# which scales the forcing.
library(bayesplot)
fit = stan("code/ar1_forcing_ss.stan",
           data = data_list,
           pars = c("sigma","pred_y","gamma","obs_sigma"),
           iter=5000,
           chains=3)
pars = rstan::extract(fit)

mcmc_areas(fit, pars = c("sigma","gamma","obs_sigma"))

## Second, we can fix the process sd, "sigma" at 1 and re-fit
fit2 = stan("code/ar1_forcing_ss_fixsig.stan",
           data = data_list,
           pars = c("pred_y","gamma","obs_sigma"),
           iter=5000,
           chains=3)
pars = rstan::extract(fit2)
mcmc_areas(fit2, pars = c("gamma","obs_sigma"))

## Third, if we knew the observation error standard deviation of these SST
# measurements, we could pass that in as a known quantity and only estimate
# 'gamma' -- this is most directly similar to the SS approach in the optim
# function above. Note: the mean estimates of the trajectory won't change if
# you vary obs_sigma, but the precision of the estimates will
data_list$obs_sigma = 0.1
fit3 = stan("code/ar1_forcing_ss_fixbothsig.stan",
            data = data_list,
            pars = c("pred_y","gamma"),
            iter=5000,
            chains=3)
pars = rstan::extract(fit3)
mcmc_areas(fit3, pars = c("gamma"))

## Now we can also compare estimates across the three models
library(broom.mixed)
y1 = tidy(fit, pars=c("pred_y"))
y1$time = seq(1,nrow(y1))
y1$model = "ss"
y2 = tidy(fit2, pars=c("pred_y"))
y2$time = seq(1,nrow(y2))
y2$model = "fixsig"
y3 = tidy(fit3, pars=c("pred_y"))
y3$time = seq(1,nrow(y3))
y3$model = "fixboth"
png("figs/bayes_estimates.png", width=6, height=8, units='in', res=300)
rbind(y1,y2,y3) %>%
  ggplot(aes(time, estimate,fill=model,col=model)) +
  geom_ribbon(aes(ymin=estimate-std.error,ymax=estimate+std.error),alpha=0.1) +
  geom_line(alpha=0.1)
dev.off()
# look at the std error rather than the estimate. this shows that
png("figs/bayes_sds_of_estimates.png", width=6, height=8, units='in', res=300)
rbind(y1,y2,y3) %>%
  ggplot(aes(time, log(std.error),fill=model,col=model)) +
  #geom_ribbon(aes(ymin=estimate-std.error,ymax=estimate+std.error),alpha=0.1) +
  geom_line(alpha=0.5) +
  ylab("Ln(std.error of estimate")
dev.off()
