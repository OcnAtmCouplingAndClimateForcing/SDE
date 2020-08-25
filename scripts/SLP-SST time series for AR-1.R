library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)

# # script for calculating GOA sst anomalies wrt 1951-1980
# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1900-01-01):1:(2019-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "data/temp.nc")
# # load and process SST data
# nc <- nc_open("data/temp.nc")

# downloading sst data in browser!
nc <- nc_open("data/nceiErsstv5_015b_edae_d976.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# these are nominally the first day of each month!

# extract study area
# 54-62 deg. N, 200-226 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

month.sst <- rowMeans(SST, na.rm = T)

# get anomalies
yr <- as.numeric(as.character(years(d)))
m <- months(d)
julian <- lubridate::yday(dates(as.character(d)))

# limit to 1950-2012
m <- m[yr %in% 1950:2012]
SST <- SST[yr %in% 1950:2012,]
julian <- julian[yr %in% 1950:2012]
yr <- yr[yr %in% 1950:2012]

f <- function(x) tapply(x, m, mean)  
mu <- apply(SST, 2, f)	# Compute monthly means for each cell
mu <- mu[rep(1:12, length(m)/12),]  # Replicate means matrix for each year at each cell

anom <- SST - mu

# and get area mean
sst.anom <- rowMeans(anom, na.rm=T)

# set up a data frame
sst.dat <- data.frame(year=yr,
                  month=m,
                  day.of.year=julian,
                  sst.anom=sst.anom)

# now slp
# downloading global NCEP/NCAR slp

# URL <- 
#   "https://upwell.pfeg.noaa.gov/erddap/griddap/noaa_esrl_118e_d5aa_117b.nc?slp[(1950-01-01):1:(2012-12-01)][(90.0):1:(20)][(0.0):1:(357.5)]"

URL <- ("https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.nc?slp[(1949-01-01):1:(2012-12-01)][(57.5):1:(47.5)][(190):1:(207.5)]")

# download.file(URL, "data/NCEP.NCAR.slp.nc")
# 
# # and test
# slp <- nc_open("data/NCEP.NCAR.slp.nc")

# there's some bug that doesn't allow nc_open to work with files that have been downloaded this way. I downloaded directly from browser
# and it works...
# I think this is a PC issue and the download.file approach will work on a Mac

slp <- nc_open("data/esrlNcepRe_2024_c172_9b57.nc")

slp

# extract dates

# seconds since 1-1-1970
raw <- ncvar_get(slp, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

x <- ncvar_get(slp, "longitude")
y <- ncvar_get(slp, "latitude")
slp <- ncvar_get(slp, "slp", verbose = F)
dim(slp) # 8 long, 5 lat, 23,346 days

# need to reverse latitude for plotting!
y <- rev(y)
slp <- slp[,5:1,]

# Change data into a matrix with months / cells for rows / columns
slp <- aperm(slp, 3:1)  
slp <- matrix(slp, nrow=dim(slp)[1], ncol=prod(dim(slp)[2:3]))  

z <- colMeans(slp, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", ylim=c(30,66), xlim=c(180, 230))

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)

# get anomalies
# I'll figure those for each day...or maybe for three-day rolling means??
julian <- lubridate::yday(dates(as.character(d)))


f <- function(x) tapply(x, julian, mean)  
mu <- apply(slp, 2, f)	# Compute daily means for each cell


# and! quick plot to check
plot(1:366,rowMeans(mu), type="l") # too noisy!

count <- function(x) sum(!is.na(x))

mu.mean <- rowMeans(mu)
mu.se <- apply(mu, 1, sd)/sqrt(apply(mu, 1, count)) 

plot.dat <- data.frame(day=1:366,
                       mean=mu.mean,
                       lowCI=mu.mean-1.96*mu.se,
                       highCI=mu.mean+1.96*mu.se)

ggplot(plot.dat, aes(day, mean)) +
  theme_bw() +
  geom_line() +
  geom_ribbon(aes(ymin=highCI, ymax=lowCI), alpha=0.2)

# that seems so weird! small CIs... is that really a meaningful pattern??

# and plot each year

year <- years(d)

plot.dat <- data.frame(day=julian,
                       year=year,
                       slp=rowMeans(slp))


sum.dat <- plot.dat %>%
  filter(year %in% c(1949:1950))

ggplot(sum.dat, aes(day, slp, color=year)) +
  theme_bw() +
  geom_line(color="light gray")

# ok! that's a lot of noise among years so the CIs I plotted appear to be too small...

# plot.dat <- plot.dat %>%
#   pivot_longer(cols=)

#########################################
# will smooth a bit...

mu3 <- rbind(mu[366,], mu, mu[1,])

# plot
plot(1:368, zoo::rollmean(rowMeans(mu3), 3, fill=NA), type="l") # still too noisy!

# try a 7-day rolling mean
mu7 <- rbind(mu[364:366,], mu, mu[1:3,])
plot(1:372, zoo::rollmean(rowMeans(mu7), 7, fill=NA), type="l") # still noisy...could fit a regr model instead of smoothing?

gam.dat <- data.frame(day=1:366,
                      mu=rowMeans(mu))

mod <- mgcv::gam(mu ~ s(day, bs='cc', k=9), data=gam.dat)

summary(mod)
plot(mod, rug=F, resid=T, pch=1)

# get predicted value
pred.mu <- predict(mod)

# and a data frame to hand off to Eric
daily.slp <- data.frame(year=as.numeric(as.character(year)), 
                        month=months(d),
                        day.of.year=julian, 
                        predicted.slp=pred.mu[match(julian, names(pred.mu))], 
                        observed.slp=rowMeans(slp))

daily.slp$anomaly.slp <- daily.slp$observed.slp-daily.slp$predicted.slp



# and combine sst-slp

daily.data <- left_join(daily.slp, sst.dat)
write.csv(daily.data, "data/daily.slp.monthly.sst.csv")
