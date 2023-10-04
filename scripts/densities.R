## Pleistocene human densities
library(boot)

# import
setwd("~/Documents/Papers/Palaeo/Cyprus/data/densities")
dat <- read.table("densities.csv", header=T, sep=",")
dat

# just Mediterranean & Levant/Near East
#dat <- dat[which(dat$region != ""), ]

## bootstrap mean
iter <- 10000
dat.err <- dat[which(is.na(dat$up)==F),]
dat.noerr <- dat[which(is.na(dat$up)==T),]

# create runif for estimates with upper/lower
med.boot <- rep(NA, iter)
for (s in 1:iter) {
  
  med.it <- rep(NA, dim(dat.err)[1])
  for (i in 1:dim(dat.err)[1]) {
    med.it[i] <- runif(n=1, min=dat.err$lo[i], max=dat.err$up[i])
  } # end i
  
  meds <- c(med.it, dat.noerr$med)
  med.boot[s] <- mean(sample(meds, size=length(meds), replace=T))
  
} # end s

med.boot.med <- median(med.boot)
med.boot.up <- quantile(med.boot, probs=0.975)
med.boot.lo <- quantile(med.boot, probs=0.025)
med.boot.up
med.boot.lo

