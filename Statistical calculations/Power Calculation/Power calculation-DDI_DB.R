############################################################################################
# AO:06/06/2023                                                                            #
#                                Simcyp Power Calculation                                  #
# https://youtu.be/VX_M3tIyiYk                                                             #
# https://youtu.be/Rsc5znwR5FA                                                             #
# Poster : https://www.certara.com/app/uploads/2017/06/Cain_2017_PAGE_Lorezapam.pdf        #
############################################################################################


rm(list=ls()) 
library("Simcyp")
library("RSQLite")


# Set script to source file location
path_user <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path_user)


#2. Initialise Simcyp

#Intialise the system files path
Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V22\\Screens\\SystemFiles",22,species = SpeciesID$Human, verbose = FALSE)



SetWorkspace("Lorezapam+Probenecid_v22.wksz") #Enter name of workspace

Simulate(database= "Lorezapam+Probenecid_v22.db") #This command is needed every time a new workspace is imported

conn<- RSQLite::dbConnect(SQLite(), "Lorezapam+Probenecid_v22.db")

nPop<-GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate) # total population size 
nPop1<-GetParameter(SimulationParameterID$Poppercent1,CategoryID$SimulationData, CompoundID$Substrate) # 1st pop size
nPop2<-GetParameter(SimulationParameterID$Poppercent2,CategoryID$SimulationData, CompoundID$Substrate) # 2nd pop size
nPop3<-GetParameter(SimulationParameterID$Poppercent3,CategoryID$SimulationData, CompoundID$Substrate) # 3rd pop size

# Data Extraction from Lorezapam+Probenecid_v22.db file ------------------------------------
# AO: This process takes a while. However, load the pre-saved data file and skip running line 44-62. 
# Load data needed. 
# save.image("MultiCpopulationAUCData.RData")
load("MultipopulationAUCData.RData")


##  START Skip-------------------------------------------
AUCSubPop<- vector()
AUCSubDDIPop<- vector()

for (i in 1:nPop){
# Substrate AUC-----------
  AUC_SubDeets<- Simcyp::GetAUCFrom_DB(ProfileID$CsysFull ,CompoundID$Substrate,individual = i,conn, allDoses = TRUE)

# Inhibtion =0 i.e substrate pk WITHOUT the presence of inhibitor
indSubj_AUC<-  AUC_SubDeets$AUC[AUC_SubDeets$Inhibition==0][1]  #OVERALL

# Inhibtion =1 i.e substrate pk WITH the presence of inhibitor
indSubj_AUC_DDI<- AUC_SubDeets$AUC[AUC_SubDeets$Inhibition==1][1] #OVERALL

AUCSubPop<-c(AUCSubPop, indSubj_AUC)
AUCSubDDIPop<-c(AUCSubDDIPop, indSubj_AUC_DDI)

##  END Skip-------------------------------------------


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
abline(v=0.285444, lty=2)  # Plot C.value

#2.Calculate mean and variance of the AUC in both populations. 
#population 1 has a mean m.pop1 and variance v.pop1
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
sample.1<-c(10,20,30,40,50)  # c(2,3,4,5,6) #sample size in population 1
sample.2<- 2*sample.1# sample size in population 2


# IF m.pop2>m.pop1, USE METHOD 1
# IF m.pop2<m.pop1, USE METHOD 2 <<< this is the Case for this script.



#METHOD 1: If m.pop2>m.pop1
#4a. Define critical value(for an alpha of 0.05) of population 1 by sample size
alpha<-0.05   
c.value<-matrix(0,1,length(sample.1)) # 

for (i in 1:length(sample.1)){
  # Understanding the difference in qnorm, pnorm, dnorm 
  # https://thomasleeper.com/Rcourse/Tutorials/distributions.html#:~:text=The%20pnorm%20function%20provides%20the,at%20a%20specified%20cumulative%20density.
  c.value[1,i]<-qnorm(1-alpha, m.pop1,sqrt(v.pop1/sample.1[i]))  # Critical region/ area.
  # Here we are getting the z (crit) value where we observe 1-aplha probability.  
  ## AO:NOTES
  #----------
  #assuming z0 is a possible value in your data. e.g in a dice, what is the probability of getting value 3
  #in context the sample size
  #probability= area under the graph at point z0.
  #pnorm= probability of Z<=z0. It only computes the area BELOW
  #qnorm is the inverse of Pnorm
  # so you if you are setting the probability/ area under the curve (set to 1-alpha ).then you are finding what z0 value is 
  
  
}
View(c.value)
# AO:
#Here we see as the number of samples increases, the closer Z0 aka C.value is to the mean of Pop1 (true mean)
#Meaning the SD/ deviation of the AUC data is reducing. Meaning more certainty if the sample size increases.
#Meaning the sample mean is less affected by extreme values. 


#5.Calculating power for each sample size

power.PK<-matrix(0,1,length(sample.2))

for (i in 1:length(sample.2)){
  power.PK[1,i]<-1-pnorm(c.value[1,i],m.pop2,sqrt(v.pop2/sample.2[i]))
  ## AO:NOTES
  #----------
  # here we are calculating pnorm. meaning we have the Z0 aka c.value now, 
  # So we are calculating the probability i.e area under the left hand side of the density curve then take away from 1.
  # pnorm= probability of Z<=z0. It only computes area BELOW
  # 1-pmorn gives you the are above which is power calc 
}

power<-power.PK*100
power

# -----------------------------------------------------------------------------------------------------------------

#METHOD 2: If m.pop2<m.pop1
#4b. Define critical value(for an alpha of 0.05) of population 1 by sample size
alpha<-0.05   
c.value<-matrix(0,1,length(sample.1)) # 

for (i in 1:length(sample.1)){
  # Understanding the difference in qnorm, pnorm, dnorm 
  # https://thomasleeper.com/Rcourse/Tutorials/distributions.html#:~:text=The%20pnorm%20function%20provides%20the,at%20a%20specified%20cumulative%20density.
  c.value[1,i]<-qnorm(alpha, m.pop1,sqrt(v.pop1/sample.1[i]))  # Critical region/ area.
  
}

#5.Calculating power for each sample size

power.PK<-matrix(0,1,length(sample.2))

for (i in 1:length(sample.2)){
  power.PK[1,i]<-pnorm(c.value[1,i],m.pop2,sqrt(v.pop2/sample.2[i]))
  ## AO:NOTES
  #----------
  # here we are calculating pnorm. meaning we have the Z0 aka c.value now, 
  # So we are calculating the probability i.e area under the left hand side of the density curve then take away from 1.
  # pnorm= probability of Z<=z0. It only computes area BELOW
  # 1-pmorn gives you the are above which is power calc 
}

power<-power.PK*100
power

#Graph of power by sample size
plot(sample.1,power, type="l", col="red",lwd=2, xlab="Sample size", ylab="Power(%)")
abline(h=8000,lty=2, lwd=3) # 80% power
