library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)
library(zoo)

# script for calculating GOA sst anomalies wrt 1951-1980
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1900-01-01):1:(2019-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "data/temp.nc")
# load and process SST data
nc <- nc_open("data/temp.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

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

# limit to 1950-2012
m <- m[yr %in% 1950:2012]
SST <- SST[yr %in% 1950:2012,]
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
                  sst.anom=sst.anom)

# now slp
# downloading global NCEP/NCAR slp

# URL <- 
#   "https://upwell.pfeg.noaa.gov/erddap/griddap/noaa_esrl_118e_d5aa_117b.nc?slp[(1950-01-01):1:(2012-12-01)][(90.0):1:(20)][(0.0):1:(357.5)]"

URL <- ("https://upwell.pfeg.noaa.gov/erddap/griddap/noaa_esrl_118e_d5aa_117b.nc?slp[(1949-01-01):1:(2012-12-01)][(57.5):1:(47.5)][(190):1:(207.5)]")

# download.file(URL, "data/NCEP.NCAR.slp.nc")

# and test
slp <- nc_open("data/NCEP.NCAR.slp.nc")
slp

# extract dates

# seconds since 1-1-1970
raw <- ncvar_get(slp, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

x <- ncvar_get(slp, "longitude")
y <- ncvar_get(slp, "latitude")
slp <- ncvar_get(slp, "slp", verbose = F)
dim(slp) # 8 long, 5 lat, 768 months

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
yr <- as.numeric(as.character(years(d)))
m <- months(d)


f <- function(x) tapply(x, m, mean)  
mu <- apply(slp, 2, f)	# Compute monthly means for each cell
mu <- mu[rep(1:12, length(m)/12),]  # Replicate means matrix for each year at each cell

anom <- slp - mu

# and get area mean
slp.anom <- rowMeans(anom, na.rm=T)

# set up a data frame
slp.dat <- data.frame(year=yr,
                      month=m,
                      slp.anom=slp.anom)

dat <- left_join(slp.dat, sst.dat)

################
# now some exploratory analyses
temp.dat <- dat %>%
  filter(year <= 1988)

l <- nrow(temp.dat)

cor0 <- cor(temp.dat$slp.anom, temp.dat$sst.anom, use="p")

c <- 1
cor1 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

c <- 2
cor2 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

c <- 3
cor3 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

c <- 4
cor4 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

c <- 5
cor5 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

c <- 6
cor6 <- cor(temp.dat$slp.anom[1:(l-c)], temp.dat$sst.anom[(1+c):l], use="p")

cor <- data.frame(lag=0:6,
                       correlation=c(cor0, cor1, cor2, cor3, cor4, cor5, cor6))

ggplot(cor, aes(lag, correlation)) +
  theme_bw() + 
  geom_bar(stat="identity", position="dodge")+
  ggtitle("GOA SSTa-Aleutian Low SLPa: cross-cor 1950-1988")

ggsave("figs/slp-sst correlation diff lags.png", height=4, width=5)

# looks like lags 1:3 host the strongest correlations

dat$slp.anom.3 <- rollmean(dat$slp.anom, 3, align = "right", fill=NA)
dat$slp.anom.3.lag1 <- c(NA, dat$slp.anom.3[1:(nrow(dat)-1)])

dat$era <- ifelse(dat$year<=1988, "1950-1988", "1989-2012")

ggplot(dat, aes(slp.anom.3.lag1, sst.anom, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method = "lm", se=F)

# save a version for testing out time-dependent AR(1) models
dat$slp.anom.lag2 <- c(NA, NA, slp.anom[1:(nrow(dat)-2)])

save.dat <- dat %>%
  select(year, month, slp.anom.lag2, sst.anom, era)

write.csv(save.dat, "data/example slp sst data.csv", row.names = F)
