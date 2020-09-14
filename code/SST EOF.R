# download SST and fit EOFs for different eras

library(ncdf4)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
library(FactoMineR)

year <- 2008
month <- "12"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1949-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(120):1:(250)]", sep="")

download.file(URL, "data/North.Pacific.ersst")

nc <- nc_open("data/North.Pacific.ersst")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

sst <- ncvar_get(nc, "sst") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
sst <- aperm(sst, 3:1)  # First, reverse order of dimensions ("transpose" array)

sst <- matrix(sst, nrow=dim(sst)[1], ncol=prod(dim(sst)[2:3]))  # Change to matrix

dimnames(sst) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# check
sst.mean <- colMeans(sst)
z <- t(matrix(sst.mean,length(y)))  
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)

# get anomalies for EOF
# identify land and remove
land <- is.na(colMeans(sst))  # Logical vector that's true over land!
X <- sst[,!land] 

m <- months(d)  
f <- function(x) tapply(x, m, mean)  
mu <- apply(X, 2, f)	# Compute monthly means for each cell
mu <- mu[rep(1:12, length(d)/12),]  # Replicate means matrix for each year at each cell

X.anom <- X - mu   # Compute matrix of anomalies!

# now EOF by era (try 20-year periods)
yr <- years(d)

# get a vector of weights (square root of the cosine of latitude)
lat.weights <- lat[!land]
weight <- sqrt(cos(lat.weights*pi/180))

EOF.1 <- svd.triplet(cov(X.anom[yr <= 1968,]), col.w=weight)
EOF.2 <- svd.triplet(cov(X.anom[yr %in% 1969:1988,]), col.w=weight)
EOF.3 <- svd.triplet(cov(X.anom[yr %in% 1989:2008,]), col.w=weight)

# EOF.1 <- svd(cor(X.anom[yr <= 1988,]))
# EOF.2 <- svd(cor(X.anom[yr > 1988,]))

# plot loadings

png("figs/SST EOF by era.png", 6, 7.5, units='in', res=300)

par(mfcol=c(3,2), mar=c(0.5, 0.5, 1.5, 0.5))
# start with EOF1

zlim <- range(EOF.1$U[,1], EOF.2$U[,1], EOF.3$U[,1])

z <- rep(NA, ncol(sst))
z[!land] <- EOF.1$U[,1]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF1 1949:1968")

z <- rep(NA, ncol(sst))
z[!land] <- EOF.2$U[,1]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF1 1969:1988")

z <- rep(NA, ncol(sst))
z[!land] <- EOF.3$U[,1]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF1 1989:2008")


# now EOF2
zlim <- range(EOF.1$U[,2], EOF.2$U[,2], EOF.3$U[,2])

z <- rep(NA, ncol(sst))
z[!land] <- EOF.1$U[,2]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF2 1949:1968")

z <- rep(NA, ncol(sst))
z[!land] <- EOF.2$U[,2]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF2 1969:1988")

z <- rep(NA, ncol(sst))
z[!land] <- EOF.3$U[,2]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64), zlim=c(-zlim[2], zlim[2]), 
      ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
mtext("EOF2 1989:2008")

dev.off()