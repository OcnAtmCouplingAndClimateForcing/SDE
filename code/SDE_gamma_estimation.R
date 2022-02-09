# use slp-sst relationships affecting Bering Sea ecosystem
# as a test case to evaluate the esimation of gamm in SDE script

library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)
library(FactoMineR)

# plot settings
theme_set(theme_bw())

# uncomment this section to fit EOF/PC

# # load slp and cell weights for calculating SLP PC1
# slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)
# 
# # weight by the size of cells
# weights <- read.csv("./data/north.pacific.slp.weights.csv", row.names = 1)
# 
# pca <- svd.triplet(cov(slp), col.w=weights[,1]) #weighting the columns
# 
# pc1 <- as.matrix(slp) %*% pca$U[,1]
# 
# # save so FactoMineR isn't needed
# pc.save <- data.frame(pc1 = pc1)
# 
# write.csv(pc.save, "./data/N._Pac_slp_pc1.csv", row.names = F)

pc1 <- read.csv("./data/N._Pac_slp_pc1.csv")

# and scale!
pc1 <- as.vector(scale(pc1))

# load mean EBS sst anomaly
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

# examine cross-correlations
ccf(sst[,2], pc1, lag.max = 10)
print(ccf(sst[,2], pc1, lag.max = 10))
# lags 1-4 are broadly similar


# fit AR model to the entire time series

# make a data frame to hold the sst and slp time series
yr <- as.numeric(as.character(chron::years(sst$date)))

# and fix incorrect years!
fix <- yr > 2030
yr[fix] <- yr[fix] - 100

m <- as.numeric(months(sst$date))

dat <- data.frame(date = lubridate::parse_date_time(x = paste(yr, m, "01"), orders="ymd", tz="America/Anchorage"),
                  sst = sst[,2],
                  slp = c(NA, pc1[1:767])) # lagging slp by one month

# and drop NAs
dat <- na.omit(dat)

# ar_ls calculates the process deviations after
# accounting for forcing variables and autocorrelation,
# (1-gamma)

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
  #s.sig = sig
  sig = sig - lm(sig ~ time)$fitted.values
  # interpolate output on the original time grid
  s.sig=(sig[-1]+sig[-T])/2 # midpoint
  # final step is normalize
  s.sig=s.sig/sd(s.sig)
  return(s.sig)
}

calc_ss = function(theta) {
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, theta)
  ss = sum((pred_ts - dat$sst)^2) # return SS for optim (minimizes by default)
}

## loop through a range of maximum values of gamma (0.01 - 1)
interval.max <- seq(0.01, 1, by = 0.01)

output <- data.frame()

for(i in 1:length(interval.max)){
# i <- 1  

# optimize by default is minimizing (with maximum = FALSE)
o = optimize(f=calc_ss, interval = c(0,interval.max[i]), maximum=FALSE)

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
                gamma = o$minimum)

pred.sst = data.frame(t = dat$date,
                      sst = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP

output <- rbind(output,
                data.frame(interval.max = interval.max[i],
                           tau = 1/o$minimum,
                           sst.correlation = cor(pred.sst$sst, pred.sst$integrated.slp),
                           sum.squares = sum(pred_ts - dat$sst[2:767])^2)) # is this correct??
# dropping 1st value of dat$sst to make the same length as pred_ts

# still getting the following warning after the loop runs??

# In pred_ts - dat$sst :
#   longer object length is not a multiple of shorter object length

}

head(output, n = 30)

# funny break - tau = 5 to tau = 22027!

# plot
output <- output %>%
  pivot_longer(cols = c(-interval.max, -tau))

ggplot(output, aes(tau, value)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scale = "free_y")

# screwy high tau values - leave these out of plot
ggplot(filter(output, tau < 6), aes(tau, value)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scale = "free_y")

# different answers? 
# tau ~ 1.7 minimizes SS, tau = 5 maximizes correlation

## try to just loop through fixed values of tau (1 to 24)

tau <- 1:24

output.2 <- data.frame()

for(i in 1:length(tau)){
  # i <- 1  
  
  # optimize by default is minimizing (with maximum = FALSE)
  # o = optimize(f=calc_ss, interval = c(0,1), maximum=FALSE)
  
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
                  gamma = 1/tau[i])
  
  pred.sst = data.frame(t = dat$date,
                        sst = dat$sst,
                        integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP
  
  output.2 <- rbind(output.2,
                  data.frame(tau = tau[i],
                             sst.correlation = cor(pred.sst$sst, pred.sst$integrated.slp),
                             sum.squares = sum(pred_ts - dat$sst[2:767])^2)) # is this correct??

  # no longer getting errors re. length of pred_ts and dat$sst!
  

}

output.2

# plot
output.2 <- output.2 %>%
  pivot_longer(cols = -tau)

ggplot(output.2, aes(tau, value)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scale = "free_y")

# again, a discrepancy between correlation and sum of squares
