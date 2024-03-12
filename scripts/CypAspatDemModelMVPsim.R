########################################################################################
## Aspatial demographic projection model for Cyprus based on Hadley CM3 NPP hindcasts
## Corey Bradshaw
## March 2024
########################################################################################

# R libraries required
library(plotrix)
library(boot)
library(tcltk)
library(MASS)
library(sp)
library(terra)
library(rgdal)
library(raster)
library(oceanmap)
library(OceanView)
library(abind)
library(pracma)
library(rgl)
library(spatstat)
library(spatialEco)
library(SpatialPack)

## source
source("~/scripts/SourceFunctions/matrixOperators.r") # matrix operators

## custom functions
# stochastic beta sampler (single sample)
stoch.beta.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# stochastic beta sampler (n samples)
stoch.n.beta.func <- function(n, mu, var) {
  Sx <- rbeta(n, (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# stochastic beta sampler for a survival vector (replaces above function for faster processing)
stoch.surv.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# dynamics model function
Nproj.func <- function(Nt, rm, K) {
  Nt1 <- round(Nt * exp(rm*(1-(Nt/K))), 0)
  return(Nt1)
}

# rescale a range
rscale <- function (x, nx1, nx2, minx, maxx) {
  nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
  return(nx)
}

# matrix rotation
rot.mat <- function(x) t(apply(x, 2, rev))

# matrix poisson resampler
rpois.fun <- function(x,y,M) {
  rpois(1,M[x,y])
}
rpois.vec.fun <- Vectorize(rpois.fun,vectorize.args = c('x','y'))

## list coordinates to xyz
coordlist2xyz <- function (list) {
  rl <- length(list[[1]]); cl <- length(list[[2]])
  coords <- c(NA,NA)
  for (r in 1:rl) {
    for (c in 1:cl) {
      coords <- rbind(coords, c(list[[1]][r],list[[2]][c]))
    }
  }
  coords <- coords[-1,]
  return(coordxyz=coords)
}


## net primary production (NPP) Hadley Centre Climate Model version 3 (HadCM3)
nppH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_NPP(absolutevalues).csv", header=T, sep=",") # 0.5 deg lat resolution
not.naH <- which(is.na(nppH[,3:dim(nppH)[2]]) == F, arr.ind=T)
upper.rowH <- as.numeric(not.naH[1,1])
lower.rowH <- as.numeric(not.naH[dim(not.naH)[1],1])
min.latH <- min(nppH[not.naH[,1], 1])  
max.latH <- max(nppH[not.naH[,1], 1])
min.lonH <- min(nppH[not.naH[,1], 2])
max.lonH <- max(nppH[not.naH[,1], 2])

as.numeric(attr(table(nppH$Lat), "names")) # lats
as.numeric(attr(table(nppH$Lon), "names")) # lons

# Cyprus region
cypr.subH <- rep(0, dim(nppH)[1])
for (n in 1:dim(nppH)[1]) {
  cypr.subH[n] <- ifelse(nppH[n,1] >= min.latH & nppH[n,1] <= max.latH & nppH[n,2] >= min.lonH & nppH[n,2] <= max.lonH, 1, 0)
}  
cypr.keepH <- which(cypr.subH == 1)
nppH.cypr <- nppH[cypr.keepH,]

sub.entryH <- which(colnames(nppH.cypr) == paste("X",14000,sep=""))
nppH.cypr.entry <- nppH.cypr[,c(1,2,sub.entryH)]

coordinates(nppH.cypr.entry) = ~ Lon + Lat
proj4string(nppH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(nppH.cypr.entry) = TRUE
nppH.entry = raster(nppH.cypr.entry)
image(nppH.entry, col=rev(grey(1:100/100)))

# transform to array
lzH <- dim(nppH.cypr)[2] - 2
nppH.array <- array(data=NA, dim=c(dim(raster2matrix(nppH.entry)),lzH))
for (k in 3:(lzH+2)) {
  nppH.cypr.k <- nppH.cypr[,c(1,2,k)] 
  coordinates(nppH.cypr.k) = ~ Lon + Lat
  proj4string(nppH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(nppH.cypr.k) = TRUE
  nppH.k = raster(nppH.cypr.k)
  nppH.array[,,k-2] <- raster2matrix(nppH.k)
}
image((nppH.array[,,5]), col=rev(grey(1:100/100)))
dim(nppH.array)

# only Cyprus proper
max.lat <- 36
min.lat <- 34.25
min.lon <- 32.25
max.lon <- 34.75

cyp.subH <- rep(0, dim(nppH)[1])
for (n in 1:dim(nppH)[1]) {
  cyp.subH[n] <- ifelse(nppH[n,1] >= min.lat & nppH[n,1] <= max.lat & nppH[n,2] >= min.lon & nppH[n,2] <= max.lon, 1, 0)
}  
cyp.keepH <- which(cyp.subH == 1)
nppH.cyp <- nppH[cyp.keepH,]

sub.entryH <- which(colnames(nppH.cyp) == paste("X",14000,sep=""))
nppH.cyp.entry <- nppH.cyp[,c(1,2,sub.entryH)]

coordinates(nppH.cyp.entry) = ~ Lon + Lat
proj4string(nppH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(nppH.cyp.entry) = TRUE
nppHcyp.entry = raster(nppH.cyp.entry)
image(nppHcyp.entry, col=rev(grey(1:100/100)))
image(nppH.entry, col=rev(grey(1:100/100)))

plot(nppH.entry)

# transform to array
lzH <- dim(nppH.cyp)[2] - 2
nppHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(nppHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  nppH.cyp.k <- nppH.cyp[,c(1,2,k)] 
  coordinates(nppH.cyp.k) = ~ Lon + Lat
  proj4string(nppH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(nppH.cyp.k) = TRUE
  nppHcyp.k = raster(nppH.cyp.k)
  nppHcyp.array[,,k-2] <- raster2matrix(nppHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((nppHcyp.array[,,ka.show]), col=rev(grey(1:100/100)))
image((nppH.array[,,ka.show]), col=rev(grey(1:100/100)))
par(mfrow=c(1,1))
dim(nppHcyp.array)

## NPP temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
nppHcyp.array.20pres <- nppHcyp.array[,,1:length(t1000Hvec)]
dim(nppHcyp.array.20pres)
image((nppHcyp.array.20pres[,,4]), col=rev(grey(1:100/100)))

cyp.nppH.mn <- cyp.nppH.lo <- cyp.nppH.up <- rep(NA,dim(nppHcyp.array.20pres)[3])
for (t in 1:dim(nppHcyp.array.20pres)[3]) {
  cyp.nppH.mn[t] <- mean(nppHcyp.array.20pres[,,t],na.rm=T)
  cyp.nppH.lo[t] <- quantile(nppHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.nppH.up[t] <- quantile(nppHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.nppH.mn, type="l", xlab="ka", ylab="NPP", ylim=c(min(cyp.nppH.lo), max(cyp.nppH.up)))
lines(t1000Hvec, cyp.nppH.lo, lty=2, col="red")
lines(t1000Hvec, cyp.nppH.up, lty=2, col="red")


## estimated sea level and change in area of Cyprus
## Global ESL reconstruction - Lambeck et al. (2014) https://doi.org/10.1073/pnas.1411762111
globESL <- read.table("~/data/seaLevel/lambeckESL.csv", sep=",", header=T)
head(globESL)
tail(globESL)

age1yr.vec <- seq(0,20,.001)

siter <- 1000
esl.mat <- matrix(data=NA, nrow=siter, ncol=length(age1yr.vec))
for (s in 1:siter) {
  esl.it <- rnorm(dim(globESL)[1], mean=globESL$esl, sd=globESL$eslsd)
  eslT.it <- data.frame(globESL$ka, esl.it)
  colnames(eslT.it) <- c("age", "esl")
  esl.mat[s, ] <- (approx(x=eslT.it$age, y=eslT.it$esl, xout=age1yr.vec))$y
} # end s
esl.mn <- apply(esl.mat, MARGIN=2, mean, na.rm=T)
esl.sd <- apply(esl.mat, MARGIN=2, sd, na.rm=T)
esl1yr.dat <- data.frame(1000*age1yr.vec, esl.mn, esl.sd)
colnames(esl1yr.dat)[1] <- "age"
head(esl1yr.dat)
tail(esl1yr.dat)

plot(esl1yr.dat$age, esl1yr.dat$esl.mn, type="l", xlab="years before present", ylab="esl (m)")
lines(esl1yr.dat$age, esl1yr.dat$esl.mn+esl1yr.dat$esl.sd, lty=2, col="red")
lines(esl1yr.dat$age, esl1yr.dat$esl.mn-esl1yr.dat$esl.sd, lty=2, col="red")

# truncate esl-age at 20 ka
esl1yr20ka.dat <- esl1yr.dat[esl1yr.dat$age <= 20000,]
tail(esl1yr20ka.dat)

write.table(esl1yr20ka.dat, "esl1yr20ka.csv", sep=",", row.names = F)


## GEBCO 2022 sea level vs. area of Cyprus
cypAesl <- read.table("~/data/seaLevel/cypareaSL.csv", sep=",", header=T)
head(cypAesl)
plot(cypAesl$esl, cypAesl$cyp.area, type="l", xlab="esl (m)", ylab="area of Cyprus (km2)")


aiter <- 1000
area.mat <- matrix(data=NA, nrow=aiter, ncol=dim(esl1yr20ka.dat)[1])
for (a in 1:aiter) {
  esl.it <- rnorm(dim(esl1yr20ka.dat)[1], esl1yr20ka.dat$esl.mn, esl1yr20ka.dat$esl.sd)
  
  area.it <- rep(NA, length(esl.it))
  for (c in 1:length(esl.it)) {
    area.it[c] <- cypAesl$cyp.area[which.min(abs(esl.it[c] - cypAesl$esl))]
  } # end c
  area.mat[a,] <- area.it
  
} # end a
areaC.mn <- apply(area.mat, MARGIN=2, mean, na.rm=T)
areaC.lo <- apply(area.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
areaC.up <- apply(area.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
areaC.sd <- apply(area.mat, MARGIN=2, sd, na.rm=T)
areaT.dat <- data.frame(esl1yr20ka.dat$age, areaC.mn, areaC.up, areaC.lo, areaC.sd)
colnames(areaT.dat)[1] <- "age"
head(areaT.dat)

plot(areaT.dat$age, areaT.dat$areaC.mn, type="l", xlab="age", ylab="area of Cyprus (km2)")
lines(areaT.dat$age, areaT.dat$areaC.mn+areaT.dat$areaC.sd, lty=2, col="red")
lines(areaT.dat$age, areaT.dat$areaC.mn-areaT.dat$areaC.sd, lty=2, col="red")

write.table(areaT.dat, "areaT.csv", sep=",", row.names = F)

par(mfrow=c(1,3))
plot(esl1yr20ka.dat$age, esl1yr20ka.dat$esl.mn, type="l", xlab="years before present", ylab="esl (m)")
lines(esl1yr20ka.dat$age, esl1yr20ka.dat$esl.mn+esl1yr20ka.dat$esl.sd, lty=2, col="red")
lines(esl1yr20ka.dat$age, esl1yr20ka.dat$esl.mn-esl1yr20ka.dat$esl.sd, lty=2, col="red")

plot(cypAesl$esl, cypAesl$cyp.area, type="l", xlab="esl (m)", ylab="area of Cyprus (km2)")

plot(areaT.dat$age, areaT.dat$areaC.mn, type="l", xlab="age", ylab="area of Cyprus (km2)")
lines(areaT.dat$age, areaT.dat$areaC.mn+areaT.dat$areaC.sd, lty=2, col="red")
lines(areaT.dat$age, areaT.dat$areaC.mn-areaT.dat$areaC.sd, lty=2, col="red")
par(mfrow=c(1,1))

## relative density, carrying capacity
# npp to K
hum.dens.med <- 6.022271e-02
hum.dens.lq <- 3.213640e-02
hum.dens.uq <- 1.439484e-01
hum.dens.max <- 1.152206e+00
hum.dens.min <- 1.751882e-02

cyp.area1ka.mn <- (approx(x=areaT.dat$age, y=areaT.dat$areaC.mn, xout=t1000Hvec*1000))$y
cyp.area1ka.lo <- (approx(x=areaT.dat$age, y=areaT.dat$areaC.lo, xout=t1000Hvec*1000))$y
cyp.area1ka.up <- (approx(x=areaT.dat$age, y=areaT.dat$areaC.up, xout=t1000Hvec*1000))$y
#cyp.area <- 9251  # area of Cyprus km2

# modify underlying K magnitude by modifying NPP across the board
K.mean <- rscale(cyp.nppH.mn, round(hum.dens.min*cyp.area1ka.mn, 0), round(hum.dens.max*cyp.area1ka.mn, 0), min(cyp.nppH.lo, na.rm=T), max(cyp.nppH.up, na.rm=T))
K.up <- rscale(cyp.nppH.up, round(hum.dens.min*cyp.area1ka.up, 0), round(hum.dens.max*cyp.area1ka.up, 0), min(cyp.nppH.lo, na.rm=T), max(cyp.nppH.up, na.rm=T))
K.lo <- rscale(cyp.nppH.lo, round(hum.dens.min*cyp.area1ka.lo, 0), round(hum.dens.max*cyp.area1ka.lo, 0), min(cyp.nppH.lo, na.rm=T), max(cyp.nppH.up, na.rm=T))

plot(t1000Hvec, K.mean, type="l", xlab="ka", ylab="K", ylim=c(min(K.lo), max(K.up)))
lines(t1000Hvec, K.lo, lty=2, col="red")
lines(t1000Hvec, K.up, lty=2, col="red")

## apply resampling procedure to produce decadal values within K limits
intyrs.vec <- seq(min(t1000Hvec*1000), max(t1000Hvec*1000), 1)

Kiter <- 10000
K.mat <- matrix(data=NA, nrow=Kiter, ncol=length(intyrs.vec))
for (i in 1:Kiter) {
  K.it <- runif(length(t1000Hvec), min=K.lo, max=K.up)
  K.expand <- approx(t1000Hvec*1000, K.it, xout=intyrs.vec)
  K.mat[i, ] <- K.expand$y
}
Ksmooth.md <- round(apply(K.mat, MARGIN=2, median, na.rm=T), 0)
Ksmooth.lo <- round(apply(K.mat, MARGIN=2, quantile, probs=0.005, na.rm=T), 0)
Ksmooth.up <- round(apply(K.mat, MARGIN=2, quantile, probs=0.995, na.rm=T), 0)
plot(intyrs.vec, Ksmooth.up, type="l", lwd=2, lty=1, col="red", xlab="years ago", ylab="K", ylim=c(min(Ksmooth.lo), max(Ksmooth.up)))
lines(intyrs.vec, Ksmooth.lo, lty=1, lwd=2, col="red")
Ksmooth.sd <- round(apply(K.mat, MARGIN=2, sd, na.rm=T), 0)
K.dat <- data.frame(intyrs.vec, Ksmooth.md, Ksmooth.up, Ksmooth.lo, Ksmooth.sd)
head(K.dat)




#######################
## demographic model ##
#######################

# Siler hazard h(x) (Gurven et al. 2007)
# average hunter-gatherer
a1 <- 0.422 # initial infant mortality rate (also known as αt)
b1 <- 1.131 # rate of mortality decline (also known as bt)
a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.47e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.086 # rate of mortality increase
longev <- 80
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
plot(x,l.x,type="l")

l.inf <- exp(-a1/b1) # survival at infinite time
T.m <- 1/b1 # time constant at which maturity is approached
h.m <- a2 # hazard for mature animals
l.m <- exp(-a2*x) # survival
h.s <- a3*exp(b3*x) # hazard for senescence
l.s <- exp((a3/b3)*(1 - exp(b3*x))) # survival for senescence
f.x <- a3*exp(b3*x)*exp((a3/b3)/(1-exp(b3*x))) # probability density function
(log(a3) - log(a1)) / a3
T.s <- (1/b3) # modal survival time

# average forager-horticuluralist
a1.fh <- 0.418; b1.fh <- 1.657; a2.fh <- 0.012; a3.fh <- 3.65e-04; b3.fh <- 0.074
h.x.fh <- a1.fh * exp(-b1.fh*x) + a2.fh + a3.fh * exp(b3.fh * x)
plot(x,h.x.fh,pch=19,type="l")
plot(x,log(h.x.fh),pch=19,type="l")
l.x.fh <- exp((-a1.fh/b1.fh) * (1 - exp(-b1.fh*x))) * exp(-a2.fh * x) * exp(a3.fh/b3.fh * (1 - exp(b3.fh * x)))
plot(x,l.x.fh,type="l")

# NT Aboriginal
a1.nta <- 0.242; b1.nta <- 1.031; a2.nta <- 0.000; a3.nta <- 7.13e-04; b3.nta <- 0.063
h.nta <- a1.nta * exp(-b1.nta*x) + a2.nta + a3.nta * exp(b3.nta * x)
plot(x,h.nta,pch=19,type="l")
plot(x,log(h.nta),pch=19,type="l")
l.x.nta <- exp((-a1.nta/b1.nta) * (1 - exp(-b1.nta*x))) * exp(-a2.nta * x) * exp(a3.nta/b3.nta * (1 - exp(b3.nta * x)))
plot(x,l.x.nta,type="l")

## survival
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
dx
qx <- dx/lx[1:(length(lx)-1)]
qx
Sx <- 1 - qx
Sx
sx <- lx[2:len.lx]/lx[1:(len.lx-1)]
mx <- 1 - sx
Lx <- (lx[1:(len.lx-1)] + lx[2:len.lx])/2
ex <- rev(cumsum(rev(Lx)))/lx[-len.lx]
ex
ex.avg <- ex + x[-len.lx]
ex.avg

# average forager-horticulturalist
lx.fh <- round(init.pop*l.x.fh,0)
len.lx.fh <- length(lx.fh)
dx.fh <- lx.fh[1:(len.lx.fh-1)]-lx.fh[2:len.lx.fh]
qx.fh <- dx.fh/lx.fh[1:(length(lx.fh)-1)]
Sx.fh <- 1 - qx.fh
sx.fh <- lx.fh[2:len.lx.fh]/lx.fh[1:(len.lx.fh-1)]
mx.fh <- 1 - sx.fh
Lx.fh <- (lx.fh[1:(len.lx.fh-1)] + lx.fh[2:len.lx.fh])/2
ex.fh <- rev(cumsum(rev(Lx.fh)))/lx[-len.lx.fh]
ex.avg.fh <- ex.fh + x[-len.lx.fh]

# average NT aboriginal
lx.nta <- round(init.pop*l.x.nta,0)
len.lx.nta <- length(lx.nta)
dx.nta <- lx.nta[1:(len.lx.nta-1)]-lx.nta[2:len.lx.nta]
qx.nta <- dx.nta/lx.nta[1:(length(lx.nta)-1)]
Sx.nta <- 1 - qx.nta
sx.nta <- lx.nta[2:len.lx.nta]/lx.nta[1:(len.lx.nta-1)]
mx.nta <- 1 - sx.nta
Lx.nta <- (lx.nta[1:(len.lx.nta-1)] + lx.nta[2:len.lx.nta])/2
ex.nta <- rev(cumsum(rev(Lx.nta)))/lx[-len.lx.nta]
ex.avg.nta <- ex.nta + x[-len.lx.nta]

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

par(mfrow=c(2,1))
plot(x[-1], Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
plot(x[-1], ex, pch=19, type="l", xlab="age (years)", ylab="life expectancy")
par(mfrow=c(1,1))
surv.out <- data.frame(x[-1],Sx,ex)
colnames(surv.out)[1] <- "x"

# plot average hunter-gatherer, average forager-horticulturalist, and NT aboriginal together
plot(x[-1], Sx, pch=19, type="l", xlab="age (years)", ylab="Sx") # average hunter-gatherer
lines(x[-1], Sx.fh, lty=2, lwd=3, col="red") # average forager-horticulturalist
lines(x[-1], Sx.fh, lty=3, lwd=3, col="green") # NT aboriginal

# fertility (Walker et al. 2006)
primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
prim.mean <- round(mean(primiparity.walker),0)
prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
print(c(prim.lo, prim.mean, prim.hi))
dat.world13 <- read.table("~/data/demography/world2013lifetable.csv", header=T, sep=",")
fert.world13 <- dat.world13$m.f
fert.trunc <- fert.world13[1:(longev+1)]
pfert.trunc <- fert.trunc/sum(fert.trunc)
fert.bentley <- 4.69/2 # Bentley 1985 for !Kung
fert.vec <- fert.bentley * pfert.trunc
plot(x,fert.vec, type="l", xlab="age (years)", ylab="fertility")

fert.out <- data.frame(x,fert.vec)
colnames(fert.out)[2] <- "fert"

## construct matrix
stages <- len.lx
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- x
rownames(popmat) <- x

## populate matrix
popmat[1,] <- fert.vec
diag(popmat[2:stages,]) <- Sx
popmat[stages,stages] <- 0 # Sx[stages-1]
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
plot(x,stable.stage.dist(popmat),type="l")
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)

## different Siler model parameters' influence on base matrix
popmat.fh <- popmat.orig # forager-horticulturalist
diag(popmat.fh[2:stages,]) <- Sx.fh
popmat.fh[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.fh)

popmat.nta <- popmat.orig # NT Aboriginal
diag(popmat.nta[2:stages,]) <- Sx.nta
popmat.nta[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.nta)

## initial population vector
pop.found <- 2750 # founding population size
init.vec <- stable.stage.dist(popmat) * pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 1 # 
#************************
yr.end <- 15000 # set projection end date
#************************
t <- (yr.end - yr.now)

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig
yr.vec <- seq(yr.now,yr.end)

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i] 
}

## year projection vector
yrs <- seq(yr.end,yr.now,-1)

# plot
#plot(yrs,log10(as.vector(colSums(n.mat))),type="l",xlab="year",ylab="log10 N",xlim=c(yr.now,yr.end))


## density feedback function on survival
popmat <- popmat.orig
Sx.mod <- 0.99615*Sx
diag(popmat[2:stages,]) <- Sx.mod
popmat[stages,stages] <- 0 # Sx[stages-1]
max.r(popmat) # rate of population change, 1-yr
surv.min <- 0.98

# initial colonisation window
col.old <- 14257
col.yng <- 13613

# generations to project
gen.proj <- 150

# iterations
iter <- 10000

# set quasi-extinction threshold
min.thresh <- 50

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- round(runif(iter,col.yng,col.old), 0) 
#************************
yr.en <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
#************************
t <- (yr.st - yr.en)

## initial population vector
popmat <- popmat.orig
init.vec <- stable.stage.dist(popmat) * pop.found/2 # stable stage distribution x founding population size (female only)

## set population storage matrices
f.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
f.mat[,1] <- init.vec

## iterate projection
itdiv <- iter/100

# set storage matrices & vectors
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))

# progress counter
pb <- txtProgressBar(min=1, max=iter, style=3)

for (e in 1:iter) { # e iterations loop
  
  yr.st.sub <- which(K.dat$intyrs.vec == yr.st[e])
  yr.en.sub <- which(K.dat$intyrs.vec == yr.en[e])
  yr.vec.run <- yr.st[e]:yr.en[e]
  K.run <- rnorm(length(yr.vec.run), mean=K.dat$Ksmooth.md[yr.st.sub:yr.en.sub], sd=K.dat$Ksmooth.sd[yr.st.sub:yr.en.sub]) # this iteration's realised K series
  #plot(yr.vec.run, K.run, type="l", xlab="year",ylab="K (N)")
  # Krun.dat <- data.frame(yr.vec.run, K.run) # write example K run
  # setwd("/Users/brad0317/Documents/Papers/Palaeo/Cyprus/out")
  # write.table(Krun.dat, "Krun.csv", sep=",", row.names = F)
  
  ## reset popmat to original values
  popmat <- popmat.orig
  
  ## set up projection loop
  for (i in 1:t[e]) { # i projection loop
    
    ## reconstruct popmat with stochastic elements
    popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
    diag(popmat[2:(stages), ]) <- ((ifelse(sum(f.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
    
    # survival (+ catastrophic mortality at 50%)
    if ((rbinom(1, 1, cat.pr)) == 1) {
      diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
    popmat[stages,stages] <- 0
    
    ## project over interval
    f.mat[,i+1] <- popmat %*% f.mat[,i] # females
    
  } # end i projection-length loop
  
  n.sums.mat[e,] <- 2 * as.vector(colSums(f.mat)) # total population (F + M)
  
  setTxtProgressBar(pb, e)
} # end e iterations loop

# N confidence limits
n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

n.min <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # minimum over all projected years
q.ext.vec <- ifelse(n.min < (min.thresh), 1, 0)
pr.ext <- sum(q.ext.vec)/iter
pr.ext

# plot
plot(yr.vec.run, (n.mn), type="l", xlab="ka", ylab="N", ylim=c((min(n.lo)),(max(n.up))))
lines(yr.vec.run, (n.lo), lty=2, col="red")
lines(yr.vec.run, (n.up), lty=2, col="red")

Nproj.out <- data.frame(yr.vec.run, n.mn, n.up, n.lo)
write.table(Nproj.out, "Nproj.csv", row.names = F, sep=",")


##################
## estimate MVP ##
##################
iter <- 10000

yr.st <- round(runif(iter,col.yng,col.old), 0) 
#************************
yr.en <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
#************************
t <- (yr.st - yr.en)

found.vec <- seq(500, 4000, 250)
pr.ext.vec <- rep(NA, length(found.vec))

for (f in 1:length(found.vec)) {

  ## initial population vector
  popmat <- popmat.orig
  init.vec <- stable.stage.dist(popmat) * found.vec[f]/2 # stable stage distribution x founding population size (female only)
  
  ## set population storage matrices
  f.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
  f.mat[,1] <- init.vec
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))
  
  # progress counter
  pb <- txtProgressBar(min=1, max=iter, style=3)
  
  for (e in 1:iter) { # e iterations loop
    
    yr.st.sub <- which(K.dat$intyrs.vec == yr.st[e])
    yr.en.sub <- which(K.dat$intyrs.vec == yr.en[e])
    yr.vec.run <- yr.st[e]:yr.en[e]
    K.run <- rnorm(length(yr.vec.run), mean=K.dat$Ksmooth.md[yr.st.sub:yr.en.sub], sd=K.dat$Ksmooth.sd[yr.st.sub:yr.en.sub]) # this iteration's realised K series

    ## reset popmat to original values
    popmat <- popmat.orig
    
    ## set up projection loop
    for (i in 1:t[e]) { # i projection loop
      
      ## reconstruct popmat with stochastic elements
      popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
      diag(popmat[2:(stages), ]) <- ((ifelse(sum(f.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
      
      # survival (+ catastrophic mortality at 50%)
      if ((rbinom(1, 1, cat.pr)) == 1) {
        diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
      popmat[stages,stages] <- 0
      
      ## project over interval
      f.mat[,i+1] <- popmat %*% f.mat[,i] # females
      
    } # end i projection-length loop
    
    n.sums.mat[e,] <- 2 * as.vector(colSums(f.mat)) # total population (F + M)
    
    setTxtProgressBar(pb, e)
  } # end e iterations loop
  
  # N confidence limits
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  n.min <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # minimum over all projected years
  q.ext.vec <- ifelse(n.min < (min.thresh), 1, 0)
  pr.ext.vec[f] <- sum(q.ext.vec)/iter

  # plot
  plot(yr.vec.run, (n.mn), type="l", xlab="ka", ylab="N", ylim=c((min(n.lo)),(max(n.up))))
  lines(yr.vec.run, (n.lo), lty=2, col="red")
  lines(yr.vec.run, (n.up), lty=2, col="red")
  
  if (f > 1) {
    plot(found.vec, pr.ext.vec, type="l", lty=2, col="red", xlab="founding N", ylab="Pr(ext)")
  }

} # end f

plot(found.vec, pr.ext.vec, type="l", lty=2, col="red", xlab="founding N", ylab="Pr(Qext)")

foundN.ext <- data.frame(found.vec, pr.ext.vec)
write.table(foundN.ext, "foundNext.csv", sep=",", row.names = F)
