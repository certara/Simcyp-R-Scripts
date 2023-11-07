require(Simcyp)
require(ggplot2)
library(RSQLite)
library(PKNCA)
library(gridExtra)
library(Rmisc)
library(plyr)
library(dplyr)
library(moments)
library(cowplot)


dir <- getwd()
setwd("C:/Users/lcurry/OneDrive - Certara/Documents/R animal applications")
Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V22\\Screens\\SystemFiles",22,species = SpeciesID$Dog, verbose = FALSE) 
Simcyp::SetWorkspace("Dogx10IV.wksz")
Simulate()
conn <- RSQLite::dbConnect(RSQLite::SQLite(), "local_copy.db")

Time <- Simcyp::GetProfile_DB(ProfileID$nTimeSubFull, compound=-1 , inhibition = FALSE, conn, individual = 1) # Time data
nSubjects <- GetParameter("idNumberGroupMembers", CategoryID$SimulationData, 0)
nPop <- GetParameter(SimulationParameterID$Pop, CategoryID$SimulationData, 0)
values <- array(data.frame(), dim = nSubjects)
for (i in 1:nPop) {
  Profile <- GetProfile_DB(ProfileID$CsysFull , individual=i, compound=CompoundID$Substrate, conn, inhibition = FALSE) # substrate conc in plasma in the absence of drug interaction for individual 1
  
  values[[i]] <- c(values[[i]], Profile)
}

df <- data.frame(sapply(values,c))
p95.obs <- apply(df, 1, quantile, probs = c(0.05, 0.95))
p95.mean <- apply(df, 1, mean)
"perc05" <- p95.obs[1, ]
"perc95" <- p95.obs[2, ]
Mean <- rowMeans(df)

Profile <- data.frame(Time, Mean, perc05, perc95)
theme_set(theme_cowplot())

plot <- ggplot(Profile, aes(x = Time)) + 
  geom_line(aes(y = Mean), colour = "black") +
  geom_line(aes(y = `perc05`), colour = "red", linetype = "dashed") + 
  geom_line(aes(y = `perc95`), colour = "red", linetype = "dashed") +
  xlab('Time (h)') + 
  ylab('Concentration (mg/L)')+
  theme(text = element_text(size = 12)) 

