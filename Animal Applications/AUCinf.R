require(Simcyp)
library(RSQLite)
library(PKNCA)
library(gridExtra)
library(Rmisc)
library(plyr)
library(dplyr)
library(moments)


dir <- getwd()
setwd(dir)
Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V22\\Screens\\SystemFiles",22,species = SpeciesID$Dog, verbose = FALSE) 

Simcyp::SetWorkspace("Dogx10IV.wksz")
Simulate(database = "local_copy.db")
conn <- RSQLite::dbConnect(RSQLite::SQLite(), "local_copy.db")
nSubjects <- GetParameter("idNumberGroupMembers", CategoryID$SimulationData, 0)
nPop <- GetParameter(SimulationParameterID$Pop, CategoryID$SimulationData, 0)
values <- array(data.frame())
for (i in 1:nPop) {
  AUC <- Simcyp::GetAUCFrom_DB(ProfileID$CsysFull, CompoundID$Substrate, i, allDoses = TRUE)
  
  auc <- split(AUC , f = AUC$Inhibition)
  AUC <- data.frame(auc[["0"]])
  AUC <- tail(AUC, 1)
  AUC <- AUC[["AUCinf"]]

  values <- c(values, AUC)
  values <- as.numeric(values)
}

Mean <- mean(values)
GeoMean <- geomean(values)
Median <- median(values)
GeoCV <- geocv(values)/100
SD <- sd(values)
conf10 <- CI(log(values), ci=0.90)
perc10geo <- exp(conf10[['lower']])
perc90geo <- exp(conf10[['upper']])
perc5 <- unname(quantile(values, probs = 0.05))
perc95 <- unname(quantile(values, probs = 0.95))
Skewness <- skewness(values)
CV <- SD/Mean
rangeLow <- min(values)
rangeHigh <- max(values)
Fold <- rangeHigh/rangeLow 

RSQLite::dbDisconnect(conn)
stats <- data.frame(Mean, Median, GeoMean, GeoCV, perc10geo, perc90geo, perc5, perc95, Skewness, CV, rangeLow, rangeHigh, Fold, SD)












