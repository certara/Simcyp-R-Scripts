require(Simcyp)
library(RSQLite)
library(reshape2)
library(PKNCA)
library(Rmisc)
library(plyr)
library(dplyr)
library(readr)
library(caTools)

dir <- getwd()
setwd(dir)
Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V23\\Screens\\SystemFiles",23,species = SpeciesID$Dog, verbose = FALSE) 
Simcyp::SetWorkspace("Dogx10IV.wksz")
Simulate(database = "local_copy.db")
conn <- RSQLite::dbConnect(RSQLite::SQLite(), "local_copy.db")

Time <- Simcyp::GetProfile_DB(ProfileID$nTimeSubFull, compound=-1 , inhibition = FALSE, conn, individual = 1) # Time data
nPop <- GetParameter(SimulationParameterID$Pop, CategoryID$SimulationData, 0)
values <- array(data.frame(), dim = nSubjects)
for (i in 1:nPop) {
  Profile <- GetProfile_DB(ProfileID$CsysFull , individual=i, compound=CompoundID$Substrate, conn, inhibition = FALSE) # substrate conc in plasma in the absence of drug interaction for individuals

  values[[i]] <- c(values[[i]], Profile)
}

df <- data.frame(sapply(values,c))
Profile <- data.frame(Time, df)

#Enter values for AUC to be calculated in hours
From <- 184 
To <- 192

AUCprof <- Profile %>% filter(between(Time, From, To))
Time <-  AUCprof$Time
AUCprof <- AUCprof[ -c(1) ]
AUCprof <- as.list(as.data.frame(AUCprof))

for(i in 1:length(AUCprof)){
  AUCprof[[i]] <- trapz(Time, AUCprof[[i]])
}

AUC <- mean(as.numeric((AUCprof)))