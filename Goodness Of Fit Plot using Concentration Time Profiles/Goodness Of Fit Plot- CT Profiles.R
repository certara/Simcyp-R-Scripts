################################################################################################
#                                                                                              #
#                                Goodness of Fit plot                                          #
#                          Concentration-time profile data                                     #
#         Using data from Fluvoxamine V22 compound summary file (Fig 1 and 2)                  #
#   3 studies were used in this example  DeBree1983,  DeVries1993 and Fleishaker1994           #
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


# Here we check if "GOFObsandPredData.RData" file exists in your working folder, then load it.
# Otherwise, we create it.
# -------------------------------------------------------------------------------------------
if(file.exists("GOFObsandPredData.RData")){

  # Load pre-saved data
  load("GOFObsandPredData.RData")
  message("Loaded GOFObsandPredData.RData file")

} else {

  # We will:
  # .Load the Human Simcyp simulator
  # .Run the simulation for each study
  # .Extract the Observed CT from each database and pair it with the Predicted CT data
  # .Combine both Observed and Predicted data into a new dataframe for GOF plotting

# Initialise Simcyp Engine
InitialiseHuman(verbose = FALSE)


# Enter workspace/study name manually
SimcypWksz <- c("Fig1A_Debree1983_100mg_single_oral.wksz",
                "Fig1B_Devries1993_100mg_single_oral.wksz",
                "Fig2_Fleishaker1994_MD.wksz")

# Study label for each workspace above; the order matters!
StudyNames <- c("DeBree1983",
                "DeVries1993",
                "Fleishaker1994")



# ----------------- Obtain Observed and Predicted Data and save into a data frame ----------------
# Dataframe should include the following columns: Study, Observed, Predicted
GOFData <- NULL


# Observed data settings
ObsScale <- 1  # Observed data stored in ng/mL
ObsDVID  <- 1  # DVID for Sub Plasma


for (Std in 1:length(SimcypWksz)){

  # Set workspace
  Workspace <- SimcypWksz[Std]
  SetWorkspace(file.path("Workspaces", Workspace))
  
  # Run simulation and save to database
  #fs::path_ext_remove() "removes the last extension and returns the rest of the path"
  DBfilename <- paste(fs::path_ext_remove(Workspace), ".db", sep="")
  DBfilepath <- file.path(path_user, DBfilename) # File path to save the database results to
  Simulate(database = DBfilepath)
  conn <- RSQLite::dbConnect(SQLite(), DBfilepath)

  # Check Predicted profile unit
  Profile_Unit <- GetProfileUnit(ProfileID$Csys, conn)
  message(paste(StudyNames[Std], "predicted profile unit:", Profile_Unit))

  # Check Observed data info held in the workspace
  PEData <- ReadPEData(NULL, scale = 1, sourceType = FileType$database, conn = conn)
  message(paste(StudyNames[Std], "observed DV info:"))
  print(PEData$DVInfo)

   
  # Extract Predicted Csys at the Observed times and pair with the Observed data
  GoFProfile <- GetObsvsPredProfileData_DB(Observations = NULL,
                                           ObsSourceType = FileType$database,
                                           Tag = ProfileID$Csys,
                                           Inhibition = FALSE,
                                           Compound = CompoundID$Substrate,
                                           conn = conn,
                                           Scale = ObsScale)

  # Keep only the DV of interest
  GoFProfile <- GoFProfile[GoFProfile$DVID == ObsDVID, ]
  if(nrow(GoFProfile) == 0){
    stop(paste("No observed data with DVID", ObsDVID, "found in", DBfilename,
               "- check the DV info printed above and update ObsDVID."))
  }

  # Create dataframe with Observed and mean Predicted concentration data
  GOFData <- rbind(GOFData,
                   data.frame(Study = factor(StudyNames[Std]),
                              Observed = GoFProfile$ObsConc,
                              Predicted = GoFProfile$SimMean*1000) # Transform to ng/mL
                   ) 

  RSQLite::dbDisconnect(conn)   # Disconnect database file
  message(paste(StudyNames[Std], "is now complete :) "))

}

save(GOFData, file = "GOFObsandPredData.RData")

}


# ------------------------------ Goodness of fit plot -----------------------

breaks=c(0.1,1,10,100)

xLabelName <- "Observed Csys [ng/mL]"
yLabelName <- "Predicted Csys [ng/mL]"
GraphTitle <- "Goodness of fit plot - Fluvoxamine"

# Graph Limits
Limits <- c(min(c(GOFData$Observed, GOFData$Predicted), na.rm=TRUE)+1,
           10^ceiling(log10(max(c(GOFData$Observed, GOFData$Predicted), na.rm=TRUE)))+100)


# PlotObsvsPred draws the log10-log10 axes, the line of identity and the
# 1.25-fold and 2-fold deviation lines.

Plot1 <- PlotObsvsPred(
  Observed = data.frame(Study = GOFData$Study, Observed = GOFData$Observed),
  Predicted = data.frame(Study = GOFData$Study, Predicted = GOFData$Predicted),
  DataType = DataType$Profile,
  title = GraphTitle,
  x_label = xLabelName,
  y_label = yLabelName,
  breaks_vec = breaks,
  limits = Limits,
  show_identity = TRUE,
  show_2fold = TRUE,
  show_1.25fold = TRUE,
  color_by_study = TRUE
)

print(Plot1)


ggsave("GOF-Fluvoxamine-CTProfiles.pdf", plot = Plot1, width = 7, height = 7, units = "in")


# ---------------- Optional: Launch Goodness of fit shiny app ------------------
# Interactive Goodness of fit shiny application for concentration-time profile data.
# Uncomment to run:
# RunShinyApp(Tag = appID$GoFProfiles)


# END
