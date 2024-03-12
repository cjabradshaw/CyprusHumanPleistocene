#######################################################################
## calibrate oldest human archaeological evidence for Cyprus,
## and use Signor-Lipps correction to estimate window of first entry
## Corey Bradshaw
## March 2024
#######################################################################

# required R libraries
library(rcarbon)
library(stringr)
devtools::install_github("FredSaltre/CRIWM/Rextinct", force=T) # install 'Rextinct' if not already installed
library(Rexinct)

# date calibrations
# single date example
singDate.mn <- 11720
singDate.er <- 240
ID <- "Beta-40380"

# calibrate (in this case, using the IntCal20 calibration curve: Krapp et al. 2021 doi:10.1038/s41597-021-01009-3)
calib.date <- calibrate(x=singDate.mn, errors=singDate.er, ids=ID, calCurves="intcal20", type="full", spdnormalised=T)
plot(calib.date, HPD=T, credMass=0.95)
summary(calib.date)
summ.calib.date <- summary(calib.date)
summ.calib.date$MedianBP

# display upper and lower confidence bounds
calib.dateU <- as.numeric(str_split(summ.calib.date$TwoSigma_BP_1, " to ", simplify=T))[1]
calib.dateU
calib.dateL <- as.numeric(str_split(summ.calib.date$TwoSigma_BP_1, " to ", simplify=T))[2]
calib.dateL

## estimate the Signor-Lipps corrected arrival date based on the calibrated re-sampled inverse-weighted McInerney (CRIWM) method
## CRIWM described in detail in Herrando-Pérez & Saltré 2024 doi:10.1016/j.quageo.2023.101489

## all Cyprus radiocarbon dates (uncalibrated) from the Near East Radiocarbon Dates (NERD) database (Palmisano et al. 2022 doi:10.5334/joad.90) 
cyp.ages <- read.table("~/data/archdates/cyprusages.csv", sep=",", header=T)
cypTS <- data.frame(cyp.ages$age, cyp.ages$err)

# save time series to working directory
write.table(cypTS, file = "cypTS.txt", row.names = FALSE, col.names = FALSE)

## apply CRIWM algorithm
criwm(chrono_data = "cypTS.txt", signor_lipps = "arr", biased=F, cal_curve = "intcal20", cal_save=T, criwm_save = T)

