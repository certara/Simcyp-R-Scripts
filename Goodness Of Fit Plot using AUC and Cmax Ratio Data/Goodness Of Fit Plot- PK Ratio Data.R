################################################################################################
#                                                                                              #
#                      Goodness of Fit plot with Guest Criteria                                #
#                              PK paramater ratio data                                         #
#                         Caffeine/Tolbutamide + Fluvoxamine                                   #
#                                                                                              #
################################################################################################

# Clear global environment
rm(list=ls())


# Load packages
library("SimcypR")
library("RSQLite")
library("ggplot2")


# Set script to source file location
path_user <- ScriptLocation()
setwd(path_user)


# Here we check if "GOFObsandPredPKData.RData" file exists in your working folder, then load it.
# Otherwise, we create it.
# -------------------------------------------------------------------------------------------
if(file.exists("GOFObsandPredPKData.RData")){

  # Load pre-saved data
  load("GOFObsandPredPKData.RData")
  message("Loaded GOFObsandPredPKData.RData file")
  
} else {

  # We will:
  # .Load the Human Simcyp simulator
  # .In a loop, run a simulation and extract the summary stats of the Predicted
  #  AUC and Cmax ratios of 4 studies. Store the results in an object called ForestData.
  # .Create the Observed data and save into a dataframe for comparison.


# Initialise the system files path
InitialiseHuman(verbose = FALSE)


# Enter workspace/study name manually
SimcypWksz <- c("CulmM2005_Caffeine_Fluvoxamine_DDI.wksz",
                "Christensen_caffeine_DDI.wksz",
                "Madsen_2001_FLUVOXDDI_75mg.wksz",
                "Madsen_2001_FLUVOXDDI_150mg.wksz")


# Extract summary stats of Predicted AUC and Cmax ratio ---------------------------------------------
# Empty vector for each workspace
ForestData <- NULL  


for (Wks in 1:length(SimcypWksz)){

  # Set workspace
  Workspace <- SimcypWksz[Wks]
  SetWorkspace(file.path("Workspaces", Workspace))

  # Run simulation and save to database
  #fs::path_ext_remove() "removes the last extension and returns the rest of the path"
  DBfilename <- paste(fs::path_ext_remove(paste(Workspace)), ".db", sep="")
  DBfilepath <- file.path(path_user, DBfilename) # File path to save the database results to
  Simulate(database = DBfilepath)
  conn <- RSQLite::dbConnect(SQLite(), DBfilepath)

  # Get the AUC and Cmax ratio for the last dose and for the full population
  PredRatio <- GetForestData_DB(Alpha = 0.05, Upper = 95, Lower = 5,
                                conn = conn, Last_Dose = TRUE,
                                AUC_Type = AUCType$AUCt,
                                Tag = ProfileID$CsysFull)

  ForestData <- data.frame(rbind(ForestData, PredRatio))

  RSQLite::dbDisconnect(conn)      # Detach database connection
  message(paste(SimcypWksz[Wks], "is now complete :) "))

}



# Predicted PK  Ratio
PKData <- split(ForestData, f = ForestData$PKparameter)
Predicted.CmaxRatio <- data.frame(PKData[["Cmax_ratio"]])
Predicted.AUCRatio <- data.frame(PKData[["AUCt_ratio"]])


# Observed PK ratio
# Note this is ordered and matched to the SimcypWksz above!
Observed.Cmaxratio <- c(1.4,2.90,NA,NA)
Observed.AUCratio <- c(13.71,6.14,1.25,1.93)

# Prep data for GOF plot
# Combine Observed and Predicted data and add labels of Victim Compound name.
GOFDataCmax <- data.frame(Compound = Predicted.CmaxRatio$VictimCompound,
                          Observed = Observed.Cmaxratio,
                          Predicted = Predicted.CmaxRatio$Mean )

GOFDataAUC <- data.frame(Compound = Predicted.AUCRatio$VictimCompound,
                         Observed = Observed.AUCratio,
                         Predicted = Predicted.AUCRatio$Mean )

save(GOFDataCmax, GOFDataAUC, file = "GOFObsandPredPKData.RData")
}

# ------------------------------ Goodness of fit plot -----------------------

breaks = c(0.1,1,10,100,1000)

GOFData <- GOFDataAUC   # To plot Cmax data , UPDATE this entry to: GOFData<- GOFDataCmax
xLabelName <- "Observed AUC Ratio"       # UPDATE!
yLabelName <- "Predicted AUC Ratio"       # UPDATE!
GraphTitle <- "Goodness of fit plot - Fluvoxamine"


# PlotObsvsPred draws the log10-log10 axes, the line of identity, the 1.25- and
# 2-fold deviation lines and the Guest criteria curves.

Plot1 <- PlotObsvsPred(
  Observed = data.frame(Study = GOFData$Compound, Observed = GOFData$Observed),
  Predicted = data.frame(Study = GOFData$Compound, Predicted = GOFData$Predicted),
  DataType = DataType$PKParam,
  title = GraphTitle,
  x_label = xLabelName,
  y_label = yLabelName,
  breaks_vec = breaks,
  limits = NULL,
  show_identity = TRUE,
  show_2fold = TRUE,
  show_1.25fold = TRUE,
  show_guest = TRUE,
  color_by_study = TRUE
)

print(Plot1)


ggsave("GOF-Fluvoxamine-AUCRatio.pdf", plot = Plot1, width = 7, height = 7, units = "in")


# ---------------- Optional: Launch Goodness of fit shiny app ------------------
# Interactive Goodness of fit shiny application for PK parameter data.
# Uncomment to run:
# RunShinyApp(Tag = appID$GoFPKParameters)


# END
