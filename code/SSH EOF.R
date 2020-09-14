# download SODA SSH and fit EOFs for different eras

library(ncdf4)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
library(FactoMineR)

URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_f29d_0f00_8b34.nc?ssh[(1948-01-15):1:(2010-12-15T00:00:00Z)][(20.25):1:(66.25)][(150.25):1:(250.25)]"
download.file(URL, "data/SODA.SSH")

nc <- nc_open("data/SODA.SSH")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "ssh") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
ssh <- aperm(ssh, 3:1)  # First, reverse order of dimensions ("transpose" array)

ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  # Change to matrix

dimnames(ssh) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# drop area N of 60N
drop <- lat > 60
ssh[,drop] <- NA

# check
SSH.mean <- colMeans(ssh)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)

# get anomalies for EOF
# identify land and remove
land <- is.na(colMeans(ssh))  # Logical vector that's true over land!
X <- ssh[,!land] 

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

EOF.1 <- svd.triplet(cov(X.anom[yr <= 1988,]), col.w=weight)
EOF.2 <- svd.triplet(cov(X.anom[yr > 1988,]), col.w=weight)

EOF.1 <- svd(cor(X.anom[yr <= 1988,]))
EOF.2 <- svd(cor(X.anom[yr > 1988,]))

# plot loadings
z <- rep(NA, ncol(ssh))
z[!land] <- EOF.1$u[,2]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)

z <- rep(NA, ncol(ssh))
z[!land] <- EOF.2$u[,2]
z <- t(matrix(z,length(y)))  
image(x,y,z, col=oceColorsPalette(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
