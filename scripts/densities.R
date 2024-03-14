## Estimated densities of Homo sapiens hunter-gatherer in the Pleistocene 
## of Europe and the Near East

## resamples the available densities to estimate boostrapped means and confidence bounds

## required R libraries
library(boot)

# import data
dat <- read.table("~/data/densities.csv", header=T, sep=",")

# include only densities estimated for the Mediterranean & Levant/Near East
#dat <- dat[which(dat$region != ""), ]

## bootstrap mean
iter <- 10000 # number of iterations for bootstrap procedure
dat.err <- dat[which(is.na(dat$up)==F),]
dat.noerr <- dat[which(is.na(dat$up)==T),]

# create random uniform sampler (runif) for estimates with upper/lower values
med.boot <- rep(NA, iter)
for (s in 1:iter) {
  
  med.it <- rep(NA, dim(dat.err)[1])
  for (i in 1:dim(dat.err)[1]) {
    med.it[i] <- runif(n=1, min=dat.err$lo[i], max=dat.err$up[i])
  } # end i
  
  meds <- c(med.it, dat.noerr$med)
  med.boot[s] <- mean(sample(meds, size=length(meds), replace=T))
  
} # end s

med.boot.med <- median(med.boot) # median bootstrapped density
med.boot.up <- quantile(med.boot, probs=0.975) # upper bootstrapped density bound
med.boot.lo <- quantile(med.boot, probs=0.025) # lower bootstrapped density bound
med.boot.up
med.boot.lo

