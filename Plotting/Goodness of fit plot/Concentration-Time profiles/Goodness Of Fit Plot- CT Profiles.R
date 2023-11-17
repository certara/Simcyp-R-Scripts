#################################################################################################
#                                                                                               #
#                                Goodness of Fit plot                                           #
# Using data from Fluvoxamine V22 compound summary file (Fig 1 and 2)                           #
# 3 studies were used in this example  deBree1983,  Devries1993 and Fig2_Fleishaker1994_MD_V22  #
#                                                                                               #
#################################################################################################

# Clear Global Environment 
rm(list=ls()) 

# To skip running the simulation and producing the database files yourself,
# Load "GOFObsandPredData.RData" below. To do this, run line 15-32: 
# Which loads the packages and set the script location. 
# Then skip to and run from Goodness of fit plot section (line 136)
# ------------------------------------------------------------------


# load packages
library("Simcyp")
library("RSQLite")
library("ggplot2")
library("scales")  # To transform the plot axis
library("tidyverse")


# Set script to source file location
path_user <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path_user)

# Load pre-saved data
load("GOFObsandPredData.RData")

#. Initialise Simcyp Engine
Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V22\\Screens\\SystemFiles",22,species = SpeciesID$Human, verbose = FALSE)

#load("GOFObsandPredData.RData")


# ------------------------- Create Observed Dataframe ----------------------------
# Data 1: deBree1983
ObsTimeDeBree1983<- c(0.75,1.50,2.00,3.00, 5.00,8.00,12.00,24.00 )
ObsCsys_ngmLDeBree1983<- c(5.605,27.264,28.736,41.745,
                           41.118,35.400,33.245,19.618)



# Data 2: Devries1993
ObsTimeDevries1993<- c(0.00,0.71,1.98,3.02,4.11,6.24,8.27,10.19,
                       12.11,16.13,24.31,32.38,48.13,72.18)
ObsCsys_ngmLDevries1993<- c(0.085,7.318,21.444,30.550,33.187,34.376,29.352,26.371,
                            23.730,21.087,14.184,9.834,5.390,2.213)



# Data 3: Fleishaker1994
ObsTimeFleishaker1994<- c(1.00,2.00,4.00,8.00,12.00,16.00,20.00,24.00,
                          48.00,72.00,96.00,120.00,144.00,168.00,192.00,216.00,
                          217.00,218.00,219.00,220.00,221.00,222.00,224.00,225.00,
                          226.00,227.00,228.00,229.00,230.00,231.00,234.00,236.00,
                          240.00,245.00,252.00)
ObsCsys_ngmLFleishaker1994<- c(2.463,9.852,19.704,22.168,18.227,15.271,13.301,11.330,
                               16.256,19.212,34.483,41.872,50.246,52.709,58.621,54.680,
                               59.606,69.951,84.729,88.177,86.695,84.729,88.177,84.729,
                               81.281,77.340,96.355,72.906,70.443,70.936,65.025,61.576,
                               58.621,56.650,50.739)





# ----------------------- Obtain Predicted Data and Save into a data frame -----------------------
# Data 1: deBree1983
# Set workspace and run simulation
SetWorkspace("Fig1A_Debree1983_100mg_single_oral_V22.wksz") #Enter name of workspace
Simulate(database= "Fig1A_Debree1983_100mg_single_oral_V22.db") #This command is needed every time a new workspace is imported
conn<- RSQLite::dbConnect(SQLite(), "Fig1A_Debree1983_100mg_single_oral_V22.db")
# Extract CT profile
CsysDB<- GetProfileAtSpecifiedTimes_DB(ObsTimeDeBree1983, ProfileID$Csys,individual = 1,inhibition=FALSE,CompoundID$Substrate,conn) #Get the Csys data for the specified time points from the data base 
PredTimeDeBree1983<-CsysDB$x
PredCsys_mgLDeBree1983<-CsysDB$y
PredCsys_ngmLDeBree1983<-PredCsys_mgLDeBree1983*1000 # unit transformation
# Disconnect database file.
RSQLite::dbDisconnect(conn)


# Data 2: Devries1993
# Set workspace and run simulation
SetWorkspace("Fig1B_Devries1993_100mg_single_oral_V22.wksz") 
Simulate(database= "Fig1B_Devries1993_100mg_single_oral_V22.db") 
conn<- RSQLite::dbConnect(SQLite(), "Fig1B_Devries1993_100mg_single_oral_V22.db")
# Extract CT profile
CsysDB<- GetProfileAtSpecifiedTimes_DB(ObsTimeDevries1993, ProfileID$Csys,individual = 1,inhibition=FALSE,CompoundID$Substrate,conn) #Get the Csys data for the specified time points from the data base 
PredTimeDevries1993<-CsysDB$x
PredCsys_mgLDevries1993<-CsysDB$y
PredCsys_ngmLDevries1993<-PredCsys_mgLDevries1993*1000 #unit transformation
# Disconnect database file.
RSQLite::dbDisconnect(conn)



# Data 3: Fleishaker1994
# Set workspace and run simulation
SetWorkspace("Fig2_Fleishaker1994_MD_V22.wksz") 
Simulate(database= "Fig2_Fleishaker1994_MD_V22.db") 
conn<- RSQLite::dbConnect(SQLite(), "Fig2_Fleishaker1994_MD_V22.db")
# Extract CT profile
CsysDB<- GetProfileAtSpecifiedTimes_DB(ObsTimeFleishaker1994, ProfileID$Csys,individual = 1,inhibition=FALSE,CompoundID$Substrate,conn) #Get the Csys data for the specified time points from the data base 
PredTimeFleishaker1994<-CsysDB$x
PredCsys_mgLFleishaker1994<-CsysDB$y
PredCsys_ngmLFleishaker1994<-PredCsys_mgLFleishaker1994*1000 #unit transformation
# Disconnect database file.
RSQLite::dbDisconnect(conn)


# ------- Combine both observed and predicted data into a new dataframe for GOF plot ------------
# Dataframe should include the following columns: Study, Observed, Predicted
DataDeBree1983<- data.frame(Study=factor(rep("deBree1983", length(ObsTimeDeBree1983))), 
                            Observed= ObsCsys_ngmLDeBree1983 , 
                            Predicted= PredCsys_ngmLDeBree1983 ) 

DataDevries1993<- data.frame(Study=factor(rep("Devries1993", length(ObsTimeDevries1993))), 
                             Observed= ObsCsys_ngmLDevries1993 , 
                             Predicted= PredCsys_ngmLDevries1993 ) 

DataFleishaker1994<- data.frame(Study=factor(rep("Fleishaker1994", length(ObsTimeFleishaker1994))), 
                                Observed= ObsCsys_ngmLFleishaker1994 , 
                                Predicted= PredCsys_ngmLFleishaker1994 ) 

GOFData<- rbind(DataDeBree1983, DataDevries1993, DataFleishaker1994)

# save.image("GOFObsandPredData.RData")



# ------------------------------ Goodness of fit plot -----------------------

breaks=c(0.1,1,10,100)

xLabelName<- "Observed Csys [ng/mL]"
yLabelName<- "Predicted Csys [ng/mL]"
GraphTitle<- "Goodness of fit plot- Fluvoxamine"

# Graph Limits
round_up<- function(x) { 10^ceiling(log10(x))}
round_down<- function(x) 10^floor(log10(x))
Limits<- c(     #May need to update!!!
  min(c(GOFData$Observed , GOFData$Predicted ), na.rm = T)+1,  
  round_up(max(c(GOFData$Observed , GOFData$Predicted), na.rm = T))+100)


 # ---------------------------GOF and Guest plot -------------------------

Plot1 <- ggplot(data=GOFData, aes(x=Observed, y=Predicted)) + #,xaxs="i", yaxs="i"
  geom_point(aes(color=factor(Study)), size=2) + 
  geom_abline(slope=1,intercept=0) +  # Line of identity 
  geom_abline(slope=1.25,intercept=0,color="gray") +  #  1.25 fold deviation
  geom_abline(slope=1/1.25,intercept=0,color="gray") +   # 1.25 fold deviation
  geom_abline(slope=2,intercept=0,color="black",linetype="dashed") +  # 2 fold deviation 
  geom_abline(slope=1/2,intercept=0,color="black",linetype="dashed") +  # 2 fold deviation 
  coord_trans(x="log10", y="log10") +
  scale_x_continuous(expand=c(0,0), name=xLabelName, breaks=breaks, labels=breaks, limits=Limits) +
  scale_y_continuous(expand=c(0,0), name=yLabelName, breaks=breaks, labels=breaks, limits=Limits) +
  # limits=c(breaks[[1]], breaks[[length(breaks)]])
  theme(aspect.ratio = 1, 
        # axis.line = element_line(size = 0.3, colour = "black", linetype=1),
        # axis.ticks = element_line(size = 1, color="black"), 
        plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5), 
        panel.border = element_rect(color = "black", fill = NA,size = 1))+ 
  ggtitle(GraphTitle)+
  theme(  # Legend position
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_line(colour="gray", size=0.2), #element_blank(), #remove major gridlines
     # panel.grid.minor = element_line(colour="gray"), #element_blank(), #remove minor gridlines
    legend.position = c(0, 1), #legend.position = "bottom"
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6), 
    legend.background = element_rect(fill='transparent'),
    legend.key=element_blank()) +
  labs(color=NULL)   # remove legend title

print(Plot1)



# plot fold errors as a ribbon ----------------------------------------------
Error1.25 <- data.frame(Observed = seq(min(Limits), max(Limits), by= 0.01)) %>%  # 1.25 fold deviation
  mutate(Predictedmax = Observed*1.25) %>% mutate(Predictedmin = Observed/1.25)


Plot2<- ggplot() + #,xaxs="i", yaxs="i"
  geom_point(data= GOFData, aes(x=Observed, y=Predicted, color=factor(Study)), size=2) + 
  geom_abline(slope=1,intercept=0) +  # Line of identity 
  geom_ribbon(data=Error1.25, aes(x=Observed,ymin=Predictedmin,ymax=Predictedmax),alpha=0.5,fill="grey")+
  geom_abline(slope=2,intercept=0,color="black",linetype="dashed") +  # 2 fold deviation 
  geom_abline(slope=1/2,intercept=0,color="black",linetype="dashed") +  # 2 fold deviation 
  coord_trans(x="log10", y="log10") +
  scale_x_continuous(expand=c(0,0), name=xLabelName, breaks=breaks, labels=breaks, limits=Limits) +
  scale_y_continuous(expand=c(0,0), name=yLabelName, breaks=breaks, labels=breaks, limits=Limits) +
  # limits=c(breaks[[1]], breaks[[length(breaks)]])
  theme(aspect.ratio = 1, 
        # axis.line = element_line(size = 0.3, colour = "black", linetype=1),
        # axis.ticks = element_line(size = 1, color="black"), 
        plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5), 
        panel.border = element_rect(color = "black", fill = NA,size = 1))+ 
  ggtitle(GraphTitle)+
  theme(  # Legend position
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    # panel.grid.major = element_line(colour="gray", size=0.2), #element_blank(), #remove major gridlines
    # panel.grid.minor = element_line(colour="gray"), #element_blank(), #remove minor gridlines
    legend.position = c(0, 1), #legend.position = "bottom"
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6), 
    legend.background = element_rect(fill='transparent'),
    legend.key=element_blank()) +
  labs(color=NULL)   # remove legend title
Plot2



