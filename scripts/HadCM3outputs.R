###############################################################################################
## Hadley CM 3 hindcasted model outputs (NPP, temperature, precipitation, + their anomalies) ##
## Corey Bradshaw
## September 2023
###############################################################################################

## libraries
library(sp)
library(rgdal)
library(raster)
library(oceanmap)
library(insol)
library(OceanView)
library(abind)
library(pracma)
library(binford)
library(rgl)
library(scatterplot3d) 
library(spatstat)
library(spatialEco)
library(SpatialPack)

## functions 
# rescale a range
rscale <- function (x, nx1, nx2, minx, maxx) {
  nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
  return(nx)
}

# matrix rotation
rot.mat <- function(x) t(apply(x, 2, rev))

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



####################################################
## set grids
####################################################

## NPP (HadCM3)
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
writeRaster(nppH.entry, filename="nppHentry.grd", format="raster")

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


## NPP anomaly (HadCM3)
nppanomH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_NPP(anomaliesvalues).csv", header=T, sep=",") # 0.5 deg lat resolution

as.numeric(attr(table(nppanomH$Lat), "names")) # lats
as.numeric(attr(table(nppanomH$Lon), "names")) # lons

# Cyprus region
nppanomH.cypr <- nppanomH[cypr.keepH,]
nppanomH.cypr.entry <- nppanomH.cypr[,c(1,2,sub.entryH)]

coordinates(nppanomH.cypr.entry) = ~ Lon + Lat
proj4string(nppanomH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(nppanomH.cypr.entry) = TRUE
nppanomH.entry = raster(nppanomH.cypr.entry)
image(nppanomH.entry, col=rev(grey(1:100/100)))

# transform to array
lzH <- dim(nppanomH.cypr)[2] - 2
nppanomH.array <- array(data=NA, dim=c(dim(raster2matrix(nppanomH.entry)),lzH))
for (k in 3:(lzH+2)) {
  nppanomH.cypr.k <- nppanomH.cypr[,c(1,2,k)] 
  coordinates(nppanomH.cypr.k) = ~ Lon + Lat
  proj4string(nppanomH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(nppanomH.cypr.k) = TRUE
  nppanomH.k = raster(nppanomH.cypr.k)
  nppanomH.array[,,k-2] <- raster2matrix(nppanomH.k)
}
image((nppanomH.array[,,5]), col=rev(grey(1:100/100)))
dim(nppanomH.array)

# only Cyprus proper
nppanomH.cyp <- nppanomH[cyp.keepH,]

sub.entryH <- which(colnames(nppanomH.cyp) == paste("X",14000,sep=""))
nppanomH.cyp.entry <- nppanomH.cyp[,c(1,2,sub.entryH)]

coordinates(nppanomH.cyp.entry) = ~ Lon + Lat
proj4string(nppanomH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(nppanomH.cyp.entry) = TRUE
nppanomHcyp.entry = raster(nppanomH.cyp.entry)
image(nppanomHcyp.entry, col=rev(grey(1:100/100)))
image(nppanomH.entry, col=rev(grey(1:100/100)))

plot(nppanomH.entry)
writeRaster(nppanomH.entry, filename="nppanomHentry.grd", format="raster", overwrite=T)

# transform to array
lzH <- dim(nppanomH.cyp)[2] - 2
nppanomHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(nppanomHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  nppanomH.cyp.k <- nppanomH.cyp[,c(1,2,k)] 
  coordinates(nppanomH.cyp.k) = ~ Lon + Lat
  proj4string(nppanomH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(nppanomH.cyp.k) = TRUE
  nppanomHcyp.k = raster(nppanomH.cyp.k)
  nppanomHcyp.array[,,k-2] <- raster2matrix(nppanomHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((nppanomHcyp.array[,,ka.show]), col=rev(grey(1:100/100)))
image((nppanomH.array[,,ka.show]), col=rev(grey(1:100/100)))
par(mfrow=c(1,1))
dim(nppanomHcyp.array)

## NPP anomaly temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
nppanomHcyp.array.20pres <- nppanomHcyp.array[,,1:length(t1000Hvec)]
dim(nppanomHcyp.array.20pres)
image((nppanomHcyp.array.20pres[,,4]), col=rev(grey(1:100/100)))

cyp.nppanomH.mn <- cyp.nppanomH.lo <- cyp.nppanomH.up <- rep(NA,dim(nppanomHcyp.array.20pres)[3])
for (t in 1:dim(nppanomHcyp.array.20pres)[3]) {
  cyp.nppanomH.mn[t] <- mean(nppanomHcyp.array.20pres[,,t],na.rm=T)
  cyp.nppanomH.lo[t] <- quantile(nppanomHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.nppanomH.up[t] <- quantile(nppanomHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.nppanomH.mn, type="l", xlab="ka", ylab="NPP anomaly", ylim=c(min(cyp.nppanomH.lo), max(cyp.nppanomH.up)))
lines(t1000Hvec, cyp.nppanomH.lo, lty=2, col="brown")
lines(t1000Hvec, cyp.nppanomH.up, lty=2, col="brown")


## TEMPERATURE (HadCM3)
tempH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_AnnualMeanTemperature(absolutevalues).csv", header=T, sep=",") # 0.5 deg lat resolution

as.numeric(attr(table(tempH$Lat), "names")) # lats
as.numeric(attr(table(tempH$Lon), "names")) # lons

# Cyprus region
tempH.cypr <- tempH[cypr.keepH,]
tempH.cypr.entry <- tempH.cypr[,c(1,2,sub.entryH)]

coordinates(tempH.cypr.entry) = ~ Lon + Lat
proj4string(tempH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(tempH.cypr.entry) = TRUE
tempH.entry = raster(tempH.cypr.entry)
image(tempH.entry, col=colorRamps::blue2red(20))

# transform to array
lzH <- dim(tempH.cypr)[2] - 2
tempH.array <- array(data=NA, dim=c(dim(raster2matrix(tempH.entry)),lzH))
for (k in 3:(lzH+2)) {
  tempH.cypr.k <- tempH.cypr[,c(1,2,k)] 
  coordinates(tempH.cypr.k) = ~ Lon + Lat
  proj4string(tempH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(tempH.cypr.k) = TRUE
  tempH.k = raster(tempH.cypr.k)
  tempH.array[,,k-2] <- raster2matrix(tempH.k)
}
image((tempH.array[,,5]), col=colorRamps::blue2red(20))
dim(tempH.array)

# only Cyprus proper
tempH.cyp <- tempH[cyp.keepH,]

sub.entryH <- which(colnames(tempH.cyp) == paste("X",14000,sep=""))
tempH.cyp.entry <- tempH.cyp[,c(1,2,sub.entryH)]

coordinates(tempH.cyp.entry) = ~ Lon + Lat
proj4string(tempH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(tempH.cyp.entry) = TRUE
tempHcyp.entry = raster(tempH.cyp.entry)
image(tempHcyp.entry, col=colorRamps::blue2red(20))
image(tempH.entry, col=colorRamps::blue2red(20))

plot(tempH.entry)
writeRaster(tempH.entry, filename="tempHentry.grd", format="raster", overwrite=T)

# transform to array
lzH <- dim(tempH.cyp)[2] - 2
tempHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(tempHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  tempH.cyp.k <- tempH.cyp[,c(1,2,k)] 
  coordinates(tempH.cyp.k) = ~ Lon + Lat
  proj4string(tempH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(tempH.cyp.k) = TRUE
  tempHcyp.k = raster(tempH.cyp.k)
  tempHcyp.array[,,k-2] <- raster2matrix(tempHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((tempHcyp.array[,,ka.show]), col=colorRamps::blue2red(20))
image((tempH.array[,,ka.show]), col=colorRamps::blue2red(20))
par(mfrow=c(1,1))
dim(tempHcyp.array)

## TEMPERATURE temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
tempHcyp.array.20pres <- tempHcyp.array[,,1:length(t1000Hvec)]
dim(tempHcyp.array.20pres)
image((tempHcyp.array.20pres[,,4]), col=colorRamps::blue2red(20))

cyp.tempH.mn <- cyp.tempH.lo <- cyp.tempH.up <- rep(NA,dim(tempHcyp.array.20pres)[3])
for (t in 1:dim(tempHcyp.array.20pres)[3]) {
  cyp.tempH.mn[t] <- mean(tempHcyp.array.20pres[,,t],na.rm=T)
  cyp.tempH.lo[t] <- quantile(tempHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.tempH.up[t] <- quantile(tempHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.tempH.mn, type="l", xlab="ka", ylab="temperature", ylim=c(min(cyp.tempH.lo), max(cyp.tempH.up)))
lines(t1000Hvec, cyp.tempH.lo, lty=2, col="red")
lines(t1000Hvec, cyp.tempH.up, lty=2, col="red")


## TEMPERATURE anomaly (HadCM3)
tempanomH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_AnnualMeanTemperature(anomaliesvalues).csv", header=T, sep=",") # 0.5 deg lat resolution

as.numeric(attr(table(tempanomH$Lat), "names")) # lats
as.numeric(attr(table(tempanomH$Lon), "names")) # lons

# Cyprus region
tempanomH.cypr <- tempanomH[cypr.keepH,]
tempanomH.cypr.entry <- tempanomH.cypr[,c(1,2,sub.entryH)]

coordinates(tempanomH.cypr.entry) = ~ Lon + Lat
proj4string(tempanomH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(tempanomH.cypr.entry) = TRUE
tempanomH.entry = raster(tempanomH.cypr.entry)
image(tempanomH.entry, col=colorRamps::blue2red(20))

# transform to array
lzH <- dim(tempanomH.cypr)[2] - 2
tempanomH.array <- array(data=NA, dim=c(dim(raster2matrix(tempanomH.entry)),lzH))
for (k in 3:(lzH+2)) {
  tempanomH.cypr.k <- tempanomH.cypr[,c(1,2,k)] 
  coordinates(tempanomH.cypr.k) = ~ Lon + Lat
  proj4string(tempanomH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(tempanomH.cypr.k) = TRUE
  tempanomH.k = raster(tempanomH.cypr.k)
  tempanomH.array[,,k-2] <- raster2matrix(tempanomH.k)
}
image((tempanomH.array[,,5]), col=colorRamps::blue2red(20))
dim(tempanomH.array)

# only Cyprus proper
tempanomH.cyp <- tempanomH[cyp.keepH,]

sub.entryH <- which(colnames(tempanomH.cyp) == paste("X",14000,sep=""))
tempanomH.cyp.entry <- tempanomH.cyp[,c(1,2,sub.entryH)]

coordinates(tempanomH.cyp.entry) = ~ Lon + Lat
proj4string(tempanomH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(tempanomH.cyp.entry) = TRUE
tempanomHcyp.entry = raster(tempanomH.cyp.entry)
image(tempanomHcyp.entry, col=colorRamps::blue2red(20))
image(tempanomH.entry, col=colorRamps::blue2red(20))

plot(tempanomH.entry)
writeRaster(tempanomH.entry, filename="tempanomHentry.grd", format="raster", overwrite=T)

# transform to array
lzH <- dim(tempanomH.cyp)[2] - 2
tempanomHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(tempanomHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  tempanomH.cyp.k <- tempanomH.cyp[,c(1,2,k)] 
  coordinates(tempanomH.cyp.k) = ~ Lon + Lat
  proj4string(tempanomH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(tempanomH.cyp.k) = TRUE
  tempanomHcyp.k = raster(tempanomH.cyp.k)
  tempanomHcyp.array[,,k-2] <- raster2matrix(tempanomHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((tempanomHcyp.array[,,ka.show]), col=colorRamps::blue2red(20))
image((tempanomH.array[,,ka.show]), col=colorRamps::blue2red(20))
par(mfrow=c(1,1))
dim(tempanomHcyp.array)

## TEMPERATURE temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
tempanomHcyp.array.20pres <- tempanomHcyp.array[,,1:length(t1000Hvec)]
dim(tempanomHcyp.array.20pres)
image((tempanomHcyp.array.20pres[,,4]), col=colorRamps::blue2red(20))

cyp.tempanomH.mn <- cyp.tempanomH.lo <- cyp.tempanomH.up <- rep(NA,dim(tempanomHcyp.array.20pres)[3])
for (t in 1:dim(tempanomHcyp.array.20pres)[3]) {
  cyp.tempanomH.mn[t] <- mean(tempanomHcyp.array.20pres[,,t],na.rm=T)
  cyp.tempanomH.lo[t] <- quantile(tempanomHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.tempanomH.up[t] <- quantile(tempanomHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.tempanomH.mn, type="l", xlab="ka", ylab="temperature anomaly", ylim=c(min(cyp.tempanomH.lo), max(cyp.tempanomH.up)))
lines(t1000Hvec, cyp.tempanomH.lo, lty=2, col="red")
lines(t1000Hvec, cyp.tempanomH.up, lty=2, col="red")


## PRECIPITATION  (HadCM3)
prcpH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_AnnualPrecipitation(absolutevalues).csv", header=T, sep=",") # 0.5 deg lat resolution

as.numeric(attr(table(prcpH$Lat), "names")) # lats
as.numeric(attr(table(prcpH$Lon), "names")) # lons

# Cyprus region
prcpH.cypr <- prcpH[cypr.keepH,]
prcpH.cypr.entry <- prcpH.cypr[,c(1,2,sub.entryH)]

coordinates(prcpH.cypr.entry) = ~ Lon + Lat
proj4string(prcpH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(prcpH.cypr.entry) = TRUE
prcpH.entry = raster(prcpH.cypr.entry)
image(prcpH.entry, col=rev(grDevices::blues9))

# transform to array
lzH <- dim(prcpH.cypr)[2] - 2
prcpH.array <- array(data=NA, dim=c(dim(raster2matrix(prcpH.entry)),lzH))
for (k in 3:(lzH+2)) {
  prcpH.cypr.k <- prcpH.cypr[,c(1,2,k)] 
  coordinates(prcpH.cypr.k) = ~ Lon + Lat
  proj4string(prcpH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(prcpH.cypr.k) = TRUE
  prcpH.k = raster(prcpH.cypr.k)
  prcpH.array[,,k-2] <- raster2matrix(prcpH.k)
}
image((prcpH.array[,,5]), col=rev(grDevices::blues9))
dim(prcpH.array)

# only Cyprus proper
prcpH.cyp <- prcpH[cyp.keepH,]

sub.entryH <- which(colnames(prcpH.cyp) == paste("X",14000,sep=""))
prcpH.cyp.entry <- prcpH.cyp[,c(1,2,sub.entryH)]

coordinates(prcpH.cyp.entry) = ~ Lon + Lat
proj4string(prcpH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(prcpH.cyp.entry) = TRUE
prcpHcyp.entry = raster(prcpH.cyp.entry)
image(prcpHcyp.entry, col=rev(grDevices::blues9))
image(prcpH.entry, col=rev(grDevices::blues9))

plot(prcpH.entry)
writeRaster(prcpH.entry, filename="prcpHentry.grd", format="raster", overwrite=T)

# transform to array
lzH <- dim(prcpH.cyp)[2] - 2
prcpHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(prcpHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  prcpH.cyp.k <- prcpH.cyp[,c(1,2,k)] 
  coordinates(prcpH.cyp.k) = ~ Lon + Lat
  proj4string(prcpH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(prcpH.cyp.k) = TRUE
  prcpHcyp.k = raster(prcpH.cyp.k)
  prcpHcyp.array[,,k-2] <- raster2matrix(prcpHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((prcpHcyp.array[,,ka.show]), col=rev(grDevices::blues9))
image((prcpH.array[,,ka.show]), col=rev(grDevices::blues9))
par(mfrow=c(1,1))
dim(prcpHcyp.array)

## PRECIPITATION temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
prcpHcyp.array.20pres <- prcpHcyp.array[,,1:length(t1000Hvec)]
dim(prcpHcyp.array.20pres)
image((prcpHcyp.array.20pres[,,4]), col=rev(grDevices::blues9))

cyp.prcpH.mn <- cyp.prcpH.lo <- cyp.prcpH.up <- rep(NA,dim(prcpHcyp.array.20pres)[3])
for (t in 1:dim(prcpHcyp.array.20pres)[3]) {
  cyp.prcpH.mn[t] <- mean(prcpHcyp.array.20pres[,,t],na.rm=T)
  cyp.prcpH.lo[t] <- quantile(prcpHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.prcpH.up[t] <- quantile(prcpHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.prcpH.mn, type="l", xlab="ka", ylab="precipitation", ylim=c(min(cyp.prcpH.lo), max(cyp.prcpH.up)))
lines(t1000Hvec, cyp.prcpH.lo, lty=2, col="blue")
lines(t1000Hvec, cyp.prcpH.up, lty=2, col="blue")


## PRECIPITATION anomaly (HadCM3)
prcpanomH <- read.table("~/data/HadCM3/CyprusRegion(20ka)_AnnualPrecipitation(anomaliesvalues).csv", header=T, sep=",") # 0.5 deg lat resolution

as.numeric(attr(table(prcpanomH$Lat), "names")) # lats
as.numeric(attr(table(prcpanomH$Lon), "names")) # lons

# Cyprus region
prcpanomH.cypr <- prcpanomH[cypr.keepH,]
prcpanomH.cypr.entry <- prcpanomH.cypr[,c(1,2,sub.entryH)]

coordinates(prcpanomH.cypr.entry) = ~ Lon + Lat
proj4string(prcpanomH.cypr.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(prcpanomH.cypr.entry) = TRUE
prcpanomH.entry = raster(prcpanomH.cypr.entry)
image(prcpanomH.entry, col=rev(grDevices::blues9))

# transform to array
lzH <- dim(prcpanomH.cypr)[2] - 2
prcpanomH.array <- array(data=NA, dim=c(dim(raster2matrix(prcpanomH.entry)),lzH))
for (k in 3:(lzH+2)) {
  prcpanomH.cypr.k <- prcpanomH.cypr[,c(1,2,k)] 
  coordinates(prcpanomH.cypr.k) = ~ Lon + Lat
  proj4string(prcpanomH.cypr.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(prcpanomH.cypr.k) = TRUE
  prcpanomH.k = raster(prcpanomH.cypr.k)
  prcpanomH.array[,,k-2] <- raster2matrix(prcpanomH.k)
}
image((prcpanomH.array[,,5]), col=rev(grDevices::blues9))
dim(prcpanomH.array)

# only Cyprus proper
prcpanomH.cyp <- prcpanomH[cyp.keepH,]

sub.entryH <- which(colnames(prcpanomH.cyp) == paste("X",14000,sep=""))
prcpanomH.cyp.entry <- prcpanomH.cyp[,c(1,2,sub.entryH)]

coordinates(prcpanomH.cyp.entry) = ~ Lon + Lat
proj4string(prcpanomH.cyp.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(prcpanomH.cyp.entry) = TRUE
prcpanomHcyp.entry = raster(prcpanomH.cyp.entry)
image(prcpanomHcyp.entry, col=rev(grDevices::blues9))
image(prcpanomH.entry, col=rev(grDevices::blues9))

plot(prcpanomH.entry)
writeRaster(prcpanomH.entry, filename="prcpanomHentry.grd", format="raster", overwrite=T)

# transform to array
lzH <- dim(prcpanomH.cyp)[2] - 2
prcpanomHcyp.array <- array(data=NA, dim=c(dim(raster2matrix(prcpanomHcyp.entry)),lzH))
for (k in 3:(lzH+2)) {
  prcpanomH.cyp.k <- prcpanomH.cyp[,c(1,2,k)] 
  coordinates(prcpanomH.cyp.k) = ~ Lon + Lat
  proj4string(prcpanomH.cyp.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
  gridded(prcpanomH.cyp.k) = TRUE
  prcpanomHcyp.k = raster(prcpanomH.cyp.k)
  prcpanomHcyp.array[,,k-2] <- raster2matrix(prcpanomHcyp.k)
}
ka.show <- 21
par(mfrow=c(1,2))
image((prcpanomHcyp.array[,,ka.show]), col=rev(grDevices::blues9))
image((prcpanomH.array[,,ka.show]), col=rev(grDevices::blues9))
par(mfrow=c(1,1))
dim(prcpanomHcyp.array)

## PRECIPITATION temporal outputs 20 ka—present (HadCM3) (Cyprus only)
t1000Hvec <- 0:20
prcpanomHcyp.array.20pres <- prcpanomHcyp.array[,,1:length(t1000Hvec)]
dim(prcpanomHcyp.array.20pres)
image((prcpanomHcyp.array.20pres[,,4]), col=rev(grDevices::blues9))

cyp.prcpanomH.mn <- cyp.prcpanomH.lo <- cyp.prcpanomH.up <- rep(NA,dim(prcpanomHcyp.array.20pres)[3])
for (t in 1:dim(prcpanomHcyp.array.20pres)[3]) {
  cyp.prcpanomH.mn[t] <- mean(prcpanomHcyp.array.20pres[,,t],na.rm=T)
  cyp.prcpanomH.lo[t] <- quantile(prcpanomHcyp.array.20pres[,,t], probs=0.025, na.rm=T)
  cyp.prcpanomH.up[t] <- quantile(prcpanomHcyp.array.20pres[,,t], probs=0.975, na.rm=T)
}
plot(t1000Hvec, cyp.prcpanomH.mn, type="l", xlab="ka", ylab="precipitation anomaly", ylim=c(min(cyp.prcpanomH.lo), max(cyp.prcpanomH.up)))
lines(t1000Hvec, cyp.prcpanomH.lo, lty=2, col="blue")
lines(t1000Hvec, cyp.prcpanomH.up, lty=2, col="blue")


had.out <- data.frame(t1000Hvec, cyp.nppanomH.mn, cyp.nppanomH.up, cyp.nppanomH.lo,
                      cyp.tempanomH.mn, cyp.tempanomH.up, cyp.tempanomH.lo,
                      cyp.prcpanomH.mn, cyp.prcpanomH.up, cyp.prcpanomH.lo,
                      cyp.nppH.mn, cyp.nppH.up, cyp.nppH.lo,
                      cyp.tempH.mn, cyp.tempH.up, cyp.tempH.lo,
                      cyp.prcpH.mn, cyp.prcpH.up, cyp.prcpH.lo)
write.table(had.out, "HadTimeSeries.csv", row.names = F, sep=",")



