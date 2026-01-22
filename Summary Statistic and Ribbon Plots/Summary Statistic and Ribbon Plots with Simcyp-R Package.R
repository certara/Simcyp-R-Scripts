# Extract Profile Statistics and Plot ----
# Clear Global Environment ----
rm(list=ls())

# 1. Load Packages -----
# Required packages for Simcyp-R
library("Simcyp")
library("RSQLite")
# Others
library("ggplot2")


# 2. Initialise Simcyp engine ----

Simcyp::Initialise(species = SpeciesID$Human)


# 3. Set script to source file location ----

path_user <-Simcyp::ScriptLocation()
setwd(path_user)


# 4. Set workspace ----

SetWorkspace("Default DDI.wksz")


# 5. Run a simulation ----

Simulate(database = "NameYourFile.db")

# Connect to simulated database file ----

conn<- RSQLite::dbConnect(SQLite(), "NameYourFile.db")


# 6. Extract  concentration-Time data for the full population ----

# Population size
nSub<- GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate)

# Create an empty vector to store profiles for each subject
Group<-vector()
Time<-vector()
Csys<- vector()
CsysInh<- vector()


for (i in 1:nSub ){   # no of subject = 100

   # Time Profile
   Current_Time<-Simcyp::GetProfile_DB(ProfileID$nTimeSub, compound=-1 , inhibition = FALSE, conn, individual = i) # Time data

   # Csys without Inhibition
   Current_Csys<-GetProfile_DB(ProfileID$Csys , individual=i, compound=CompoundID$Substrate, conn, inhibition = FALSE) # substrate conc in plasma in the absence of drug interaction
   # Csys with Inhibition
   Current_CsysInh<-GetProfile_DB(ProfileID$Csys , individual=i, compound=CompoundID$Substrate, conn, inhibition = TRUE) # substrate conc in plasma in the presence of drug interaction

   # combine
   Group <- c(Group, rep(i,length(Current_Csys)))
   Time <- c(Time,Current_Time)
   Csys<- c(Csys, Current_Csys)
   CsysInh<-  c(CsysInh, Current_CsysInh)

}

# Combine all the data and save as a dataframe
SimCSys<- data.frame(factor(Group), Time, Csys, CsysInh)
colnames(SimCSys) <- c('IND','Time_h','CSys_mgL', 'CSysInh_mgL')


# 7. Plot individual profiles -----
colors <- c("CSys" = "red", "CSysInh" = "blue")

ggplot(data = SimCSys) +
   geom_line(aes(x=Time_h, y=CSys_mgL,  group = IND, color = "CSys"))+
   geom_line(aes(x=Time_h, y=CSysInh_mgL,  group = IND, color = "CSysInh"))+
   scale_x_continuous("Time(h)", limits = c(0,24)) +
   scale_y_continuous("Csys (mg/L)")+ #Note: Simcyp-R reports profile outputs as mg/L even if not used in the workspace
   ggtitle("Individual Systemic Concentration Time Profile")+
   labs(color = "Profile") +
   scale_colour_manual(values = colors)


# 8. Extract the Mean, Lower and Upper percentile of the simulated Systemic concentration Profile ----
# Population level
ProfilePopSummary<- GetProfileStats_DB(ProfileID$Csys,CompoundID$Substrate, Inhibition=FALSE, Upper= 95, Lower=5, Trial=FALSE, conn)
View(ProfilePopSummary)

ProfilePopSummaryInh<- GetProfileStats_DB(ProfileID$Csys,CompoundID$Substrate, Inhibition=TRUE, Upper= 95, Lower=5, Trial=FALSE, conn)
View(ProfilePopSummaryInh)

# Trial level
ProfileTrialSummary<- GetProfileStats_DB(ProfileID$Csys,CompoundID$Substrate, Inhibition=FALSE, Upper= 95, Lower=5, Trial=TRUE, conn)
View(ProfileTrialSummary)

ProfileTrialSummaryInh<- GetProfileStats_DB(ProfileID$Csys,CompoundID$Substrate, Inhibition=TRUE, Upper= 95, Lower=5, Trial=TRUE, conn)
View(ProfileTrialSummaryInh)

# 9. Create a Ribbon Plot -----
# Sub Profile Mean, 5th and 95th percentile
PlotProfilePercentile_DB(ProfileId = ProfileID$Csys, CompoundId = CompoundID$Substrate, Inhibition = FALSE, Upper = 95, Lower = 5, conn, Col="red")

# Sub+Inh Profile Mean, 5th and 95th percentile
PlotProfilePercentile_DB(ProfileId = ProfileID$Csys, CompoundId = CompoundID$Substrate, Inhibition = TRUE, Upper = 95, Lower = 5, conn, Col=c("red","blue"))

# Turn off Ribbon plot
PlotProfilePercentile_DB(ProfileId = ProfileID$Csys, CompoundId = CompoundID$Substrate, Inhibition = TRUE, Upper = 95, Lower = 5, conn, Col=c("red","blue"), ShowRibbon=FALSE)


# 10. Other Summary stats functions ----

# Demographic statistic:GetIndividualValueStats_DB ----

# Compound Statistic:GetCompoundResultStats_DB ----

# PD statistic:GetPDResultStats_DB ----

# 11. Disconnect the database connection ----
RSQLite::dbDisconnect(conn)


# 12. Disconnect the Simcyp Simulator Engine from R session ----
Simcyp::Uninitialise()

# END ----

