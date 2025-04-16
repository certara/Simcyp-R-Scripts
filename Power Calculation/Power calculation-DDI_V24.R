############################################################################################
# AO:16/04/2024                                                                            #
#                                Simcyp Power Calculation                                  #
############################################################################################


rm(list=ls()) 
library("Simcyp")
library("RSQLite")


# Set script to source file location
path_user <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path_user)


#2. Initialise Simcyp

#Intialise the system files path
Simcyp::Initialise(species = SpeciesID$Human, verbose = FALSE)



SetWorkspace("Lorezapam+Probenecid_v24.wksz") #Enter name of workspace

Simulate(database= "Lorezapam+Probenecid.db") 

conn<- RSQLite::dbConnect(SQLite(), "Lorezapam+Probenecid.db")

nPop<-GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate) # total population size 
nPop1<-GetParameter(SimulationParameterID$Poppercent1,CategoryID$SimulationData, CompoundID$Substrate) # 1st pop size
nPop2<-GetParameter(SimulationParameterID$Poppercent2,CategoryID$SimulationData, CompoundID$Substrate) # 2nd pop size
nPop3<-GetParameter(SimulationParameterID$Poppercent3,CategoryID$SimulationData, CompoundID$Substrate) # 3rd pop size

# AUC Extraction from simulated DDI study (Lorezapam+Probenecid) ------------------------------------
# AO: This process takes a while. Therefore, load the pre-saved data file.

if(file.exists("MultipopulationAUCData.RData")){
  load("MultipopulationAUCData.RData")
  
} else {
  
  # Substrate AUC-----------
  
  AUC_SubDeets<- GetAUCFrom_DB(ProfileID$CsysFull ,CompoundID$Substrate,individual = c(1:nPop),conn, allDoses = TRUE)
  # Sort the AUC_SubDeets dataframe by Individual
  AUC_SubDeets <- AUC_SubDeets[order(AUC_SubDeets$Individual),]
  
  # Inhibtion =0 i.e substrate pk WITHOUT the presence of inhibitor
  AUCSubPop<- AUC_SubDeets %>% filter(Dose==-1, Inhibition==0) %>% pull(AUC)    #OVERALL
  
  # Inhibtion =1 i.e substrate pk WITH the presence of inhibitor
  AUCSubDDIPop<- AUC_SubDeets %>% filter(Dose==-1, Inhibition==1) %>% pull(AUC)  #OVERALL
  
 }

#idAUC=AUC first dose 
#4.Get AUC values for each individual for previously run workspace
AUC.sub<-AUCSubPop
AUC.sub.healthy<-AUC.sub[1:1000]
AUC.sub.renal.36<-AUC.sub[1001:2000]
AUC.sub.renal.l3<-AUC.sub[2001:3000]

AUC.sub.inb<-AUCSubDDIPop
AUC.sub.inb.healthy<-AUC.sub.inb[1:1000]
AUC.sub.inb.renal.36<-AUC.sub.inb[1001:2000]
AUC.sub.inb.renal.l3<-AUC.sub.inb[2001:3000]

# 
# AUC.1<-AUC.sub.renal.l3         #Sub AUC
# AUC.2<-AUC.sub.inb.renal.l3     #Interaction AUC
AUC.1<-AUC.sub.healthy        #Sub HV
AUC.2<-AUC.sub.renal.l3    #Sub RI.l3
AUC.3<-AUC.sub.renal.36    #Sub RI.l3

#Plot graphs of population values
plot(density(AUC.1),col=1, lwd=2, main="", xlab="AUC", ylim=c(0, 6), xlim=c(0,0.8))
lines(density(AUC.2),col=2, lwd=2)
lines(density(AUC.3),col=3, lwd=2)
legend("topright",c("AUC Sub HV","AUC Sub RI.l3","AUC Sub RI.36"), col=c(1,2,3), lwd=2, lty=1)#, "AUC Sub RI 36" ,3
# legend("topright",c("AUC substrate","AUC substrate+inhibitor"), col=c(1,2), lwd=2, lty=1)
abline(v=0.285444, lty=2)

#2.Calculate mean and variance of the AUC in both populations. 
# population 1 has a mean m.pop1 and variance v.pop1
# population 2 has a mean m.pop2 and variance v.pop2

m.pop1<-mean(AUC.1)
m.pop2<-mean(AUC.2)
m.pop3<-mean(AUC.3)

v.pop1<-var(AUC.1)
v.pop2<-var(AUC.2)
v.pop3<-var(AUC.3)

sd.pop1<-sqrt(v.pop1)
sd.pop2<-sqrt(v.pop2)
sd.pop3<-sqrt(v.pop3)


#3.Define Sample sizes
#Here it assumes that sampling equal numbers from both populations

#Define vector of sample sizes
#sample<-c(2,3,4,5,10,20,50)
sample.1<-c(10,20,30,40,50)  # c(2,3,4,5,6) #sample size in population 1
#sample.2<-c(2,4,6,8,10)# sample size in population 2
sample.2<-sample.1


# IF m.pop2>m.pop1, USE METHOD 1
# IF m.pop2<m.pop1, USE METHOD 2



#METHOD 1: If m.pop2>m.pop1
#4a. Define critical value(for an alpha of 0.05) of population 1 by sample size
alpha<-0.05   
c.value<-matrix(0,1,length(sample.1)) # 

for (i in 1:length(sample.1)){
   c.value[1,i]<-qnorm(1-alpha, m.pop1,sqrt(v.pop1/sample.1[i]))  # Critical region/ area.
  # Here we are getting the z (crit) value where we observe 1-aplha probability.  
  
}


#5.Calculating power for each sample size

power.PK<-matrix(0,1,length(sample.2))

for (i in 1:length(sample.2)){
  power.PK[1,i]<-1-pnorm(c.value[1,i],m.pop2,sqrt(v.pop2/sample.2[i]))
}

power<-power.PK*100


# -----------------------------------------------------------------------------------------------------------------

#METHOD 2: If m.pop2<m.pop1
#4b. Define critical value(for an alpha of 0.05) of population 1 by sample size
alpha<-0.05   
c.value<-matrix(0,1,length(sample.1)) # 

for (i in 1:length(sample.1)){
  c.value[1,i]<-qnorm(alpha, m.pop1,sqrt(v.pop1/sample.1[i]))  # Critical region/ area.
  
}

#5.Calculating power for each sample size

power.PK<-matrix(0,1,length(sample.2))

for (i in 1:length(sample.2)){
  power.PK[1,i]<-pnorm(c.value[1,i],m.pop2,sqrt(v.pop2/sample.2[i]))
}

power<-power.PK*100

#Graph of power by sample size
plot(sample.1,power, type="l", col="red",lwd=2, xlab="Sample size", ylab="Power(%)")
abline(h=8000,lty=2, lwd=3) # 80% power
