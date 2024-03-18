#========================================================================================
#
#
# STEP-BY-STEP guideline to generate the Signor-Lipps-corrected estimate of human arrival
# F.SALTRE - 18/03/2024
#
#
#========================================================================================

# description:
#*************
#
# We accessed the 10 oldest archaeological sites in Cyprus and their oldest dated artefacts (and their standard deviations) from the Near East Radiocarbon Dates (NERD) database 
# The only site we did not include was Vretsia Roudias because there are no available radiocarbon dates from this site
# however, there are optically stimulated luminescence dates available for Vretsia Roudias, the oldest of which is 12.8 ± 1.6 ka 
# To these 10 radiocarbon dates (see Results), we applied the Rextinct package in R first to calibrate the dates from radiocarbon years to years before present based on the IntCal20 Northern Hemisphere calibration curve  
# We then estimate a Signor-Lipps-corrected window of human arrival using the calibrated re-sampled inverse-weighted McInerney (CRIWM) resampling algorithm (22), the most advanced and up-to-date method available that incorporates dating uncertainty in its calculations. 
# The Signor-Lipps correction based on resampling approaches generally requires > 8 dates to provide robust estimates of initial or terminal phenomena from dated archaeological or palaeontological evidence


##+++++++++++++++++++++++++++++++++++++++++++++
##+ STEP 0: INSTALLING RELEVANT LIBRARY
##+ ++++++++++++++++++++++++++++++++++++++++++++
#installing the package Rextinct + loading the package
devtools::install_github("FredSaltre/CRIWM/Rextinct", force=T)
library(Rextinct)

##+++++++++++++++++++++++++++++++++++++++++++++
##+ STEP 1: IMPORT THE ARCHAEOLOGICAL DATA
##+ ++++++++++++++++++++++++++++++++++++++++++++
# downloading oldest radiocarbon dates from 10 sites around Cyprus FROM GITHUB 
# URL to the raw CSV file on GitHub
url <- "https://raw.githubusercontent.com/FredSaltre/CyprusHumanPleistocene/main/data/archdates/cyprusages.csv"

# Specify the full path including the filename where you want to save the file
#dest_file <- "/Users/fredsaltre/Desktop/GlobalEcologyLab/Corey/Cyprus/ArrivalEstimate/NonSpatial_Estimates(CRIWM)/cyprusages.csv"
dest_file <- " xxxxxx/xxxxxx/xxx/cyprusages.csv" #user to add a directory to save the file the end of the parth MUST be /cyprusages.csv

# Download the file
download.file(url, dest_file)

##+++++++++++++++++++++++++++++++++++++++++++++
##+ STEP 2: FORMATTING THE DATASET TO USE CRIWM
##+ ++++++++++++++++++++++++++++++++++++++++++++
# Loading the whole .csv file 
Full.dat <- read.csv("cyprusages.csv",header=T,sep=",",dec=".",na.strings="na")

# creating a dataframe with only the data that will be then save as a Textfile to be used by CRIWM
# must contain two columns without row names or column headings, 
# column 1 includes dates, column 2 includes errors. 
# Unlimited number of dates (rows). All dates must be Before Present, where Present = 1950.

mat<-data.frame(Full.dat$age,Full.dat$err) # dataframe with 2 columns that are: 1) uncalibrated dates and 2) dating errors

#saving the dataframe as a textfile in the specified directory named "dest_file"
write.table(mat, file = "Cyprus_TimeSeries.txt", row.names = FALSE, col.names = FALSE)

##+++++++++++++++++++++++++++++++++++++++++++++
##+ STEP 3: RUN CRIWM
##+ ++++++++++++++++++++++++++++++++++++++++++++
# we run CRIWM using the following argument
# chrono_data = "Cyprus_TimeSeries.txt" : Text file read from working directory
# signor_lipps = "arr" : CRIWM accounts for the Signor-Lipps effect at oldest (arrival) end of a time series of dated observations
# biased = F : return the unbiased estimates (independant of alpha parameter)
# radiocarbon = "all" :  assumes that all records in the time series are non-calibrated-radiocarbon dated observations.
# cal_curve = "intcal20" : uses intcal20 calibration curve for calibrating non-calibrated radiocarbon dates


out<-criwm(chrono_data = "Cyprus_TimeSeries.txt", signor_lipps = "arr", biased = F, radiocarbon = "all", cal_curve = "intcal20")

#display the result of arrival time in BP
out$criwm[2,]

