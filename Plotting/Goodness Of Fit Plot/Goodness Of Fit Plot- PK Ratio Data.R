################################################################################################
#                                                                                              #
#                             Goodness of Fit plot with Guest Criteria                         #
#                                           Fluvoxamine                                        #
#                                                                                              #   
################################################################################################
#  Clear Global envrionment
rm(list=ls())

# Load packages
library("Simcyp")
library("RSQLite")
library("tidyverse") 

#Initialise the system files path
Simcyp::Initialise(species = SpeciesID$Human, verbose = FALSE) 


# Set script to source file location 
path_user <-Simcyp::ScriptLocation()
setwd(path_user)


# Enter workspace/study name manually. 
SimcypWksz<- c("CulmM2005_Caffeine_Fluvoxamine_DDI.wksz",
               "Christensen_caffeine_DDI.wksz",
               "Madsen_2001_FLUVOXDDI_75mg.wksz",
               "Madsen_2001_FLUVOXDDI_150mg.wksz")  


# Extract summary stats of Predicted AUC and Cmax ratio ---------------------------------------------
# Empty vector for each workspace
ForestData<- NULL  # 


for (Wks in 1:length(SimcypWksz)){ 
  
  Workspace<- SimcypWksz[Wks]
  SetWorkspace(file.path("V23 workspace/",Workspace))
  
  
  ###### Run simulation and save to database
  #fs::path_ext_remove() "removes the last extension and returns the rest of the path".
  DBfilename <- paste(fs::path_ext_remove(paste(Workspace)) ,".db",sep="")  
  DBfilepath <- file.path(path_user, DBfilename) # File path to save the database results to
  Simulate(database = DBfilepath)
  conn <- RSQLite::dbConnect(SQLite(),DBfilepath)
  
  # Get the AUC and Cmax ratio for the last dose and for the full population  !!
  PredRatio<- GetForestData_DB(Alpha = 0.05, Upper = 95, Lower = 5, conn, Last_Dose = TRUE, AUC_Type = "AUCt")
  
  ForestData<- data.frame(rbind(ForestData, PredRatio))
  
  RSQLite::dbDisconnect(conn)      #Detach database connection
  message( paste( SimcypWksz[Wks], "is now complete :) "  ) )
  
}

# Predicted PK  Ratio
PKData <- split( ForestData , f = ForestData$PKparameter)
Predicted.CmaxRatio <- data.frame(PKData[["Cmax_ratio"]])
Predicted.AUCRatio <- data.frame(PKData[["AUCt_ratio"]])


# Observed PK ratio
# Note this is ordered and matched to the SimcypWksz above!
observed.Cmaxratio <- c(1.4,2.90,NA,NA)
observed.AUCratio <- c(13.71,6.14,1.25,1.93)

# Compound Name; order matters!
# make sure they are grouped by compound names. 
Compound <- c("Caffeine",      
              "Caffeine",
              "Tolbutamide",
              "Tolbutamide")

# Prep data for GOF plot. combine observed and Predicted data
GOFDataCmax<- data.frame(Compound=Compound, Observed= observed.Cmaxratio , 
           Predicted= Predicted.CmaxRatio$Mean ) 

GOFDataAUC<- data.frame(Compound=Compound, Observed= observed.AUCratio , 
                         Predicted= Predicted.AUCRatio$Mean ) 

GOFData<-GOFDataAUC # change PK here!
# GOFData<-GOFDataCmax  


# ------------------------------ Goodness of fit plot -----------------------

# Guest criteria parameters
GuestData<- function(breaks){
  
  xLow <- breaks[1]
  xHigh <- breaks[length(breaks)]
  
  stepSize <- xLow/10
  
  RobsLow <- seq(xLow,1,stepSize)
  LimXlow <- (1+2*((1/RobsLow)-1))/(1/RobsLow)
  
  LimLowUpper <- RobsLow*LimXlow
  LimLowLower <- RobsLow/LimXlow
  
  RobsHigh <- seq(from=1,to=xHigh,by=stepSize)
  LimXhigh <- (1+2*(RobsHigh-1))/RobsHigh
  
  LimHighUpper <- RobsHigh*LimXhigh
  LimHighLower <- RobsHigh/LimXhigh
  
  Robs <- c(RobsLow,RobsHigh)
  LimUpper <- c(LimLowUpper,LimHighUpper)
  LimLower <- c(LimLowLower, LimHighLower)
  
  guestUpper <- data.frame(Robs, LimUpper)
  guestLower <- data.frame(Robs, LimLower)
  
  return(list( Upper= guestUpper, Lower=guestLower))
}

breaks=c(0.1,1,10,100, 1000)

xLabelName<- "Observed AUC Ratio"       # UPDATE!
yLabelName<- "Predicted AUC Ratio"       # UPDATE!
GraphTitle<- "Goodness of fit plot- Fluvoxamine"

# Graph Limits
round_up<- function(x) { 10^ceiling(log10(x))}
round_down<- function(x) 10^floor(log10(x))
Limits<- c(
  round_down(min(c(GOFData$Observed , GOFData$Predicted ), na.rm = T)),  #
  round_up(max(c(GOFData$Observed , GOFData$Predicted), na.rm = T)))


# ---------------------------GOF and Guest plot -------------------------

Plot1 <- ggplot(data=GOFData, aes(x=Observed, y=Predicted)) + #,xaxs="i", yaxs="i"
  geom_point(aes(color=factor(Compound)), size=3) + 
  geom_line(data=GuestData(breaks)$Upper, aes(x=Robs,y=LimUpper), color="red") +
  geom_line(data=GuestData(breaks)$Lower, aes(x=Robs,y=LimLower), color="red") +
  geom_abline(slope=1,intercept=0) +  # Line of identity 
  # geom_ribbon(data=Fold_Error_Guest,aes(x=x,ymin=ymin1,ymax=ymax1),alpha=0.5,fill="grey") #  1.25 fold deviation
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
    panel.grid.major = element_line(colour="gray", size=0.2), #element_blank(), #remove major gridlines
    # panel.grid.minor = element_line(colour="gray"), #element_blank(), #remove minor gridlines
    legend.position = c(0, 1), #legend.position = "bottom"
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6), 
    legend.background = element_rect(fill='transparent'),
    legend.key=element_blank()) +
  labs(color=NULL)   # remove legend title

print(Plot1)


ggsave("GOF-Fluvoxamine-AUCRatio.pdf")


# END 
