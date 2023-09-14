library(rcarbon)
library(stringr)

# date calibrations
singDate.mn <- 11720
singDate.er <- 240
ID <- "Beta-40380"

calib.date <- calibrate(x=singDate.mn, errors=singDate.er, ids=ID, calCurves="intcal20", type="full", spdnormalised=T)
plot(calib.date, HPD=T, credMass=0.95)
summary(calib.date)
summ.calib.date <- summary(calib.date)
summ.calib.date$MedianBP

calib.dateU <- as.numeric(str_split(summ.calib.date$TwoSigma_BP_1, " to ", simplify=T))[1]
calib.dateU
calib.dateL <- as.numeric(str_split(summ.calib.date$TwoSigma_BP_1, " to ", simplify=T))[2]
calib.dateL

## CRIWM
library(Rexinct)

## all Cyprus dates
setwd("~/Documents/Papers/Palaeo/Cyprus/data")
cyp.ages <- read.table("cyprusages.csv", sep=",", header=T)
cypTS <- data.frame(cyp.ages$age, cyp.ages$err)

# save time series to working directory
write.table(cypTS, file = "cypTS.txt", row.names = FALSE, col.names = FALSE)

criwm(chrono_data = "cypTS.txt", signor_lipps = "arr", biased=F, cal_curve = "intcal20", cal_save=T, criwm_save = T)

