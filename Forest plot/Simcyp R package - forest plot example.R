# Simcyp R package - forest plot example

# This script will show you: 
# 1. How to run the Simcyp Simulator from R
# 2. How to extract the necessary data for forest plots
# 3. How to use those data to make a forest plot


# Housekeeping ----------------------------------------------------------------
library(Simcyp)
library(RSQLite)
library(tidyverse)
library(PKNCA)
library(tictoc)

# Note: You will need to have the Simcyp Simulator open to simulate. I'm using
# version 23 here.

# Setting up working directory
path_user <- Simcyp::ScriptLocation()


# Simulating ------------------------------------------------------------------

setwd(path_user)

# Initialize the Simulator. Even if you're not simulating, you'll need to do
# this to use the Simcyp package. We're based in the UK, so, USA folks, please
# note the spelling here! :)
Simcyp::Initialise(filePath = "C:/Program Files/Simcyp Simulator V23/Screens/SystemFiles",
                   requestedVersion = 23, 
                   species = SpeciesID$Human) 

# If you've already run this script in the past -- maybe you're only interested
# in tweaking the appearance of your forest plot, for example -- then you don't
# need to re-simulate everything. This script is set up to just load the results
# if you've already simulated previously.

if(file.exists("SV-Atazanavir Forest Data.RData")){
   load("SV-Atazanavir Forest Data.RData")
} else {
   
   # We will run 5 DDI simulations and store their results; specifically, we
   # need the AUC and Cmax ratios for the forest plot.
   
   # Get workspace names
   SimcypWksz <- list.files(pattern = "wksz")
   
   # Creating lists to store data
   ForestData <- list()
   SimTime <- list() # This is optional; it just tracks how long each iteration takes.
   
   for (Wks in SimcypWksz){ 
      
      tic(msg = paste("Simulation", Wks))
      
      # Setting what workspace to simulate
      SetWorkspace(Wks)
      # This will give you useful summary information on your simulation. Check
      # that you're simulating what you think you're simulating! :)
      
      # Run the simulation and save to database. For saving our database files,
      # we'll use the same file name as the workspace but with the ".db" extension.
      DBfilename <- sub("wksz", "db", Wks)
      Simulate(database = DBfilename)
      
      # Make a connection with the database file using RSQLite
      conn <- RSQLite::dbConnect(SQLite(), DBfilename)
      
      # Extract the population statistics of the predicted AUC and Cmax ratios 
      ForestData[[Wks]] <- GetForestData_DB(Alpha = 0.1,
                                            Upper = 95,
                                            Lower = 5, 
                                            conn, 
                                            Last_Dose = TRUE, 
                                            AUC_Type = "AUCt")
      
      SimTime[[Wks]] <- toc(log = TRUE)
      
      # Disconnect from the database
      RSQLite::dbDisconnect(conn)
      
   }
   
   # Putting all the forest data into a single data.frame
   ForestData <- bind_rows(ForestData)
   
   # Check the times for each iteration
   map(SimTime, "toc")
   
   save(ForestData, SimTime, SimcypWksz, 
        file = "SV-Atazanavir Forest Data.RData")
   
}

# Setting up observed data ---------------------------------------------------

# If you would like to compare your simulated data with observed in the forest
# plot, you'll need to set up a data.frame with your observed data and note
# which simulation file to compare the observed data to. Let's check out how
# those data are set up because we'll set up the observed data in the same way:
glimpse(ForestData)

# Columns we'll need to match in the observed data: 

# File
# We'll need to specify which simulation file should be matched to the observed
# data, so we'll include a column titled "File" to do that. Please note that
# there is no file extension, though.

# PKparameter
# In the forest data we extracted from those simulations, you'll see that we
# have named the PK parameters either "AUCt_ratio" or "Cmax_ratio". We'll use
# those exact same names for the PK parameters in the observed data so that the
# forest_plot function knows which parameters match. 

# GeoMean
# What we plan to graph are the geometric means and geometric confidence
# intervals (this is the most typical set of statistics for forest plots), so
# we'll need the matching data in the observed data.frame. We don't have
# observed confidence interval data, though, so we'll skip those. If we had
# them, they would be shown as error bars on the plot. 

# What if you don't have all the observed data? That's fine. Skip whatever you
# don't have and that simulation just won't have any observed data points. 

ObsRatios <- data.frame(File = c("Neely_raltegravirDDI",
                                 "Zhu_2010_raltegravirDDI",
                                 "Krishna_raltegravirDDI",
                                 "Iwamoto_2008_raltegravirDDI",
                                 "Mummaneni_Clarith_ATZ_DDI"), 
                        PKparameter = c(rep("AUCt_ratio", 5), 
                                        rep("Cmax_ratio", 5)), 
                        GeoMean = c(1.72, 1.54, 1.67, 1.72, 1.94, # AUC ratios
                                    1.37, 1.39, 1.16, 1.53, 1.50)) # Cmax ratios


# Making the forest plot -------------------------------------------------------

# Please do check out the help file. :)
?PlotForestDDI

# Let's make draft forest plot with mostly default parameters.
PlotForestDDI(SimForestData = ForestData, 
              ObsForestData = ObsRatios)

# As you can see, the y axis is labeled according to the simulation, which does
# not make for the prettiest of graph labels. We can make nicer labels either by
# by specifying which simulation needs which label with a named character
# vector, where the names are the simulation files without the file extension.
# Note that the y axis labels *must* be unique. Also note that you can include a
# carriage return with "\n".
YAxisLabels<- c("Neely_raltegravirDDI" = "Raltegravir\nNeely et al. (2010)",
                "Zhu_2010_raltegravirDDI" = "Raltegravir\nZhu et al. (2010)",
                "Krishna_raltegravirDDI" = "Raltegravir\nKrishna et al. (2016)",
                "Iwamoto_2008_raltegravirDDI" = "Raltegravir\nIwamoto et al. (2008)",
                "Mummaneni_Clarith_ATZ_DDI" = "Clarithromycin\nMummaneni et al. (2002)")

# And the revised forest plot with the nicer labels:
PlotForestDDI(SimForestData = ForestData, 
              ObsForestData = ObsRatios,
              y_axis_column = YAxisLabels)

# There are a number of ways that you can customize your forest plot, and we
# hope you'll play around with the function a bit with assistance from the help
# file. Here are a few more options for an example.

PlotForestDDI(SimForestData = ForestData, 
              ObsForestData = ObsRatios,
              y_axis_column = YAxisLabels, 
              Mean_type = "Geometric",
              Variability_type = "CI",
              
              # You can change the set of colors to one of the sets already
              # included (see the help file) or specify your own colors. We'll
              # use one of the built-in sets.
              color_set = "yellow to red",
              
              # By default, the PlotForestDDI function orders the simulations
              # from strongest inhibition at the top to strongest induction at
              # the bottom, but you can specify that you want the order to be as
              # it is in the source data.frame instead. 
              y_order = "as is", 
              
              x_axis_title = "Geometric Mean Ratio (90% CI)",
              x_axis_limits = c(1, 2),
              graph_title = "Predicted Geometric Mean AUC and Cmax Ratios for DDI studies",
              graph_title_size = 14 ,
              legend_position = "right", 
              
              # We recommend playing around with the output graph dimensions to
              # make sure things are clear and not smooshed together.
              fig_height = 7, 
              fig_width = 8,
              save_graph = "Forest plot SV-Atazanavir.png")


# Finishing up -----------------------------------------------------------------

Uninitialise()   # Uninitialise the Simcyp engine 



