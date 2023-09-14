##########################################################################################
## Signor-Lipps corrected estimates of extinction of pygmy hippos & elephants in Cyprus ##
## Corey Bradshaw
## September 2023
##########################################################################################

library(dplyr) 
#devtools::install_github("FredSaltre/CRIWM/Rextinct")
library(Rexinct)

setwd("/Users/brad0317/Documents/Papers/Palaeo/Cyprus/data/palaeo")

# source functions
source("qualityRating.R")
source("EndRating.R")


## pygmy hippo Phanourios minor
phanourios <- read.table("phanourios.txt", header = T, sep="\t")
phanourios[is.na(phanourios)] <- 'na' #replaces all the missing data with "na"to avoid TRUE/FALSE errors
phanourios$C14_CNRatioValue <-as.numeric(as.character(phanourios$C14_CNRatioValue)) #makes CNRatioValue numeric
phanourios$C14_CNRatioValue[is.na(phanourios$C14_CNRatioValue)] <- '0' #replaces missing data for CNRatioValue with 0
phanourios$C14_NPercentage <-as.numeric(as.character(phanourios$C14_NPercentage)) #makes NPercentage nurmeric
phanourios$C14_NPercentage[is.na(phanourios$C14_NPercentage)] <- '0' #replaces missing data for NPercentage with 0

# quality rating
phanourios.rated <- FosSahulQR(phanourios)
phanourios.endrated <- EndRating(phanourios.rated)

phanourios.merged <- merge(phanourios, phanourios.endrated, by="AgeID")

phanourios.out <- data.frame(phanourios.merged$AgeID, phanourios.merged$Age, phanourios.merged$Precision, phanourios.merged$preQuality)
colnames(phanourios.out) <- c("ID", "age", "err", "rating")
phanourios.out
phanourios.good <- subset(phanourios.out, rating!="C")
phanourios.good
phanourios.sort <- phanourios.good[order(phanourios.good[,2], decreasing=F),]
phanourios.sort

phanouriosTS <- data.frame(phanourios.sort$age, phanourios.sort$err)

# save time series to working directory
write.table(phanouriosTS, file = "phanouriosTS.txt", row.names = FALSE, col.names = FALSE)

# CRIWM
criwm(chrono_data = "phanouriosTS.txt", signor_lipps = "ext", biased=F, radiocarbon="all", cal_curve = "intcal20", cal_save=T, criwm_save = T)


# pygmy elephant Elaphus cypriotes
elaphus <- read.table("elaphus.txt", header = T, sep="\t")
elaphus[is.na(elaphus)] <- 'na' #replaces all the missing data with "na"to avoid TRUE/FALSE errors
elaphus$C14_CNRatioValue <-as.numeric(as.character(elaphus$C14_CNRatioValue)) #makes CNRatioValue numeric
elaphus$C14_CNRatioValue[is.na(elaphus$C14_CNRatioValue)] <- '0' #replaces missing data for CNRatioValue with 0
elaphus$C14_NPercentage <-as.numeric(as.character(elaphus$C14_NPercentage)) #makes NPercentage nurmeric
elaphus$C14_NPercentage[is.na(elaphus$C14_NPercentage)] <- '0' #replaces missing data for NPercentage with 0

# quality rating
elaphus.rated <- FosSahulQR(elaphus)
elaphus.endrated <- EndRating(elaphus.rated)

elaphus.merged <- merge(elaphus, elaphus.endrated, by="AgeID")
elaphus.merged$qualRating <- ifelse(elaphus.merged$preRating == "C", "C", elaphus.merged$preQuality)

elaphus.out <- data.frame(elaphus.merged$AgeID, elaphus.merged$Age*1000, elaphus.merged$Precision*1000, elaphus.merged$qualRating)
colnames(elaphus.out) <- c("ID", "age", "err", "rating")
elaphus.out
elaphus.good <- subset(elaphus.out, rating=="B" | rating=="A")
elaphus.good
elaphus.sort <- elaphus.good[order(elaphus.good[,2], decreasing=F),]
elaphus.sort

elaphusTS <- data.frame(elaphus.sort$age, elaphus.sort$err)

# save time series to working directory
write.table(elaphusTS, file = "elaphusTS.txt", row.names = FALSE, col.names = FALSE)

# CRIWM
criwm(chrono_data = "elaphusTS.txt", signor_lipps = "ext", biased=F, radiocarbon=0, cal_save=T, criwm_save = T)
