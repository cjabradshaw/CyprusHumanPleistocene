##########################################################################################
## Signor-Lipps corrected estimates of extinction of dwarf hippos & elephants in Cyprus ##
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


## dwarf hippo Phanourios minor
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


# dwarf elephant Elaphus cypriotes
elaphas <- read.table("elaphas.txt", header = T, sep="\t")
elaphas[is.na(elaphas)] <- 'na' #replaces all the missing data with "na"to avoid TRUE/FALSE errors
elaphas$C14_CNRatioValue <-as.numeric(as.character(elaphas$C14_CNRatioValue)) #makes CNRatioValue numeric
elaphas$C14_CNRatioValue[is.na(elaphas$C14_CNRatioValue)] <- '0' #replaces missing data for CNRatioValue with 0
elaphas$C14_NPercentage <-as.numeric(as.character(elaphas$C14_NPercentage)) #makes NPercentage nurmeric
elaphas$C14_NPercentage[is.na(elaphas$C14_NPercentage)] <- '0' #replaces missing data for NPercentage with 0

# quality rating
elaphas.rated <- FosSahulQR(elaphas)
elaphas.endrated <- EndRating(elaphas.rated)

elaphas.merged <- merge(elaphas, elaphas.endrated, by="AgeID")
elaphas.merged$qualRating <- ifelse(elaphas.merged$preRating == "C", "C", elaphas.merged$preQuality)

elaphas.out <- data.frame(elaphas.merged$AgeID, elaphas.merged$Age*1000, elaphas.merged$Precision*1000, elaphas.merged$qualRating)
colnames(elaphas.out) <- c("ID", "age", "err", "rating")
elaphas.out
elaphas.good <- subset(elaphas.out, rating=="B" | rating=="A")
elaphas.good
elaphas.sort <- elaphas.good[order(elaphas.good[,2], decreasing=F),]
elaphas.sort

elaphasTS <- data.frame(elaphas.sort$age, elaphas.sort$err)

# save time series to working directory
write.table(elaphasTS, file = "elaphasTS.txt", row.names = FALSE, col.names = FALSE)

# CRIWM
criwm(chrono_data = "elaphasTS.txt", signor_lipps = "ext", biased=F, radiocarbon=0, cal_save=T, criwm_save = T)
