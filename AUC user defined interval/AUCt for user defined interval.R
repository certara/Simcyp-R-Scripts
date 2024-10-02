require(Simcyp)
library("RSQLite")
library(PKNCA)
library(gridExtra)
library(XML)
library(magrittr)
library(plyr)
library(DBI)
library(dplyr)
library(pracma)



if (exists("debugSource"))
   Simcyp::Initialise("C://Program files//Simcyp simulator v23//screens//systemfiles", 23)
setwd("~/V23 Cobicistat")


Simcyp::SetWorkspace("Mikus_2020_150mgCOBI_midazolam-DDI.wksz")
Simulate(database = "MikusMDZ.db")
conn <- RSQLite::dbConnect(RSQLite::SQLite(), "MikusMDZ.db")

nPop <- GetParameter(SimulationParameterID$Pop, CategoryID$SimulationData, 0)
nTrials <- GetParameter("idNumberGroup", CategoryID$SimulationData, 
                        0)

auc_values <- array(data.frame(), dim = nTrials)
cmax_values <- array(data.frame(), dim = nTrials)

for (i in 1:nPop) {
   nTrial <- GetIndividualValue_DB(i, "idGroupNo", conn = conn)
   Profile <- GetProfile_DB(ProfileID$CsysFull , individual=i, compound=CompoundID$Substrate, conn, inhibition = FALSE) # substrate conc in plasma in the absence of drug interaction for individuals
   Profile_inh <- GetProfile_DB(ProfileID$CsysFull , individual=i, compound=CompoundID$Substrate, conn, inhibition = TRUE) # substrate conc in plasma in the absence of drug interaction for individuals
   Time <- Simcyp::GetProfile_DB(ProfileID$nTimeSubFull, compound=-1 , inhibition = FALSE, conn, individual = 1) # Time data
   
   df1 <- data.frame(sapply(Profile,c))
   Profile1 <- data.frame(Time, df1)
   
   From <- 26.5 
   To <- 28.5
   
   prof1 <- Profile1 %>% filter(between(Time, From, To))
   Time <-  prof1$Time
   prof1 <- prof1[ -c(1) ]
   AUCprof1 <- as.list(as.data.frame(prof1))
   
   for(i in 1:length(prof1)){
      AUCprof1[[i]] <- trapz(Time, AUCprof1[[i]])
   }
   
   Time <- Simcyp::GetProfile_DB(ProfileID$nTimeSubFull, compound=-1 , inhibition = FALSE, conn, individual = 1) # Time data
   df2 <- data.frame(sapply(Profile_inh,c))
   Profile2 <- data.frame(Time, df2)
   prof2 <- Profile2 %>% filter(between(Time, From, To))
   Time <-  prof2$Time
   prof2 <- prof2[ -c(1) ]
   AUCprof2 <- as.list(as.data.frame(prof2))
   
   for(i in 1:length(prof2)){
      AUCprof2[[i]] <- trapz(Time, AUCprof2[[i]])
   }    
   
   auc_ratio <- mapply(FUN = `/`, AUCprof2, AUCprof1, SIMPLIFY = FALSE)
   auc_ratio <- t(data.frame(auc_ratio))
   auc_values[[nTrial]] <- c(auc_values[[nTrial]], auc_ratio)
   
   
   colMax <- function(data) sapply(data, max, na.rm = TRUE)
   Cmax1 <- data.frame(colMax(prof1))
   Cmax2 <- data.frame(colMax(prof2))
   cmax_ratio <- Cmax2/Cmax1
   cmax_ratio <- unlist(cmax_ratio)
   
   cmax_values[[nTrial]] <- c(cmax_values[[nTrial]], cmax_ratio)
   
}

cmax <- sapply(cmax_values, geomean)
AUC <- sapply(auc_values, geomean)
cmax_mean <- geomean(cmax)
cmax_high <- max(cmax)
cmax_low <- min(cmax)
AUC_mean <- geomean(AUC)
AUC_high <- max(AUC)
AUC_low <- min(AUC)
name <-c("Mikus et al. 2020 (B)*")
observed.cmax <- c(NA)
observed.AUCratio <- c(8.79)
diffCmax <- vector()
diffAUC <- vector()
predicted.cmax <- vector()
predicted.AUCratio <- vector()
diffCmax <- c(diffCmax, sprintf("%0.2f", cmax_mean/observed.cmax))
diffAUC <- c(diffAUC, sprintf("%0.2f", AUC_mean/observed.AUCratio))
predicted.cmax <- c(predicted.cmax, sprintf("%0.2f( %0.2f to %0.2f)", 
                                            cmax_mean, cmax_low, cmax_high))
predicted.AUCratio <- c(predicted.AUCratio, sprintf("%0.2f( %0.2f to %0.2f)", 
                                                    AUC_mean, AUC_low, AUC_high))
RSQLite::dbDisconnect(conn)

df1 <- data.frame(name, observed.cmax, observed.AUCratio, 
                  predicted.cmax, predicted.AUCratio, diffCmax, diffAUC)
colnames(df1) <- c("Study", "Cmax Ratio", "AUC Ratio", "Cmax Ratio", "AUC Ratio", "Cmax Ratio", "AUC ratio")


