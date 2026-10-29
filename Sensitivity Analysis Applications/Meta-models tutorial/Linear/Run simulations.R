
#This is an optional script in the workflow and runs Simcyp models corresponding to each design point and extracts and store PK summary measures

#####################################################################################################################################

## 1. Initialization
## Here we read in the R packages for analysis and read in the Simcyp workspace to be studied in GSA. 
## We run an initial simulation to check that we can run a simulation, write a database, and make changes to model parameters.

## Clear memory and remove all objects: 
rm(list=ls())

## Load Required Packages for the analysis
library("SimcypR")
library("RSQLite")


## Initialise the version V25 engine for a human workspace
SimcypR::Initialise(species = SpeciesID$Human, requestedVersion = 25, verbose = FALSE) 

## Set script to source file location
path_user = SimcypR::ScriptLocation()
setwd(path_user)

## Import the Simcyp workspace that fully defines the PBPK model for it's baseline parameterization. 
## In the below the workspace is in the same directory as the script.
name.file = c("Theophylline_V25")
SimcypR::SetWorkspace(paste0(name.file, ".wksz"))

## Create a database file name.
DBfilename = paste0(name.file, ".db")
## File path to save the database results to. The database is saved to the same directory as the workspace and script.
DBfilepath = file.path(path_user, DBfilename)                 

## Run a simulation and save the result to database. 
SimcypR::Simulate(database = DBfilepath)

## RSQLite database connection: Needed to extract values from a database file. 
conn = RSQLite::dbConnect(SQLite(),DBfilename)

## Check that model parameter may be overwritten and read
SimcypR::SetParameter(Tag='idFu1', Category=CategoryID$Compound, CompoundID$Substrate, Value=0.6)
SimcypR::GetParameter('idFu1', Category=CategoryID$Compound, CompoundID$Substrate)

##Disconnect from the database.
RSQLite::dbDisconnect(conn)  

#####################################################################################################################################

## 2. Design
## Here we read tabulate information about the model parameters under study and create the designs needed for Morris Test and the Sobol-Saltelli variance-based GSA.
## We additionally generate a space-filling design for exploring the bounding behavior of the model.

## Five parameters are studied in this work as detailed below.
## Vmax
## VSS
## ka
##B/P
##fu

## A mix of normal and uniform distributions is assumed here to demonstrate flexibility of approach - this readily extends.
## For normal distributions the lower and upper bounds correspond to 3 standard deviations

names = c("Vmax", 
          "VSS", 
          "ka", 
          "B/P",
          "fu")
id    = c("/Workspace/Compounds/0/CYPRoutes/1/Vmax", 
          "idEnteredVss",
          "idka",
          "idBP",
          "idFu1")


lowB  = c(0.05, 0.05, 0.2, 0.55, 0.45)
upB   = c(50, 10, 40, 2, 0.55)


# This dataframe contain parameter names, their Simcyp tags, lower and upper lower bounds.
para.input = data.frame(names, id, lowB, upB)
n.param = 5

## Here we read in the design  run the model and extract and store output for subsequent analysis.
## Strip out column 1, which is a point label

DesignModel = read.csv('Theo Design.csv')
DesignModel = DesignModel[,2:6]


#####################################################################################################################################


## 3. Run the PBPK models and collate output

## Design point 1
## Overwrite default parameters. In this loop we overwrite each of the five parameters using the first line of the design.

for (j in 1:n.param){
  SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(DesignModel[1,j]))
}

## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
DBfilepath <- tempfile()

## Runs simulation and connect to database
Simulate(database = DBfilepath)
conn = RSQLite::dbConnect(SQLite(), DBfilepath)


#Output times for the simulation.
Time    = SimcypR::GetProfile_DB(ProfileID$nTimeSub, Compound=-1 , individual = 1, conn, Inhibition = FALSE)

##Extract plasma concentration
Plasma <- SimcypR::GetProfile_DB(ProfileID$Csys, Compound=CompoundID$Substrate, individual=1, conn, Inhibition = FALSE) 

Store_Plasma <- data.frame(Plasma)


RSQLite::dbDisconnect(conn)

unlink(DBfilepath)

# We loop through the remaining design points.
# In each case we overwrite parameters, run the simulation, extract outputs and append the profiles to the output dataframes.

points = dim(DesignModel)[1]
for(i in 2:points){
  
  for (j in 1:n.param){
    SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(DesignModel[i,j]))
  }
  
  ## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
  DBfilepath <- tempfile()
  
  ## Runs simulation and connect to database
  Simulate(database = DBfilepath)
  conn = RSQLite::dbConnect(SQLite(), DBfilepath)
  
  ## Vascular-Kidney
  Plasma <- SimcypR::GetProfile_DB(ProfileID$Csys, Compound=CompoundID$Substrate, individual=1, conn, Inhibition = FALSE) 
  

  #Keep track of the number of simulations run
  Store_Plasma <- cbind(Store_Plasma, Plasma)
    print(i)
  
}

## Entire set of simulations for Sobol variance-based analysis is complete and stored in four dataframes.
## Dataframes could be output and worked with directly but the files are too large for a practical example.
## The code below processes model output and extracts AUC, Cmax and Tmax for each summary measure.

#Initialize empty vectors for results
AUC = c()
Cmax = c()
Tmax = c()

for(i in 1:points)
{
  run = Store_Plasma[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}


## Results is a dataframe with three output variables corresponding to the design
## We write as .csv, which can be read in for the meta model-based GSA
Results = data.frame(AUC, Cmax, Tmax)

write.csv(Results, 'Results.csv')


#####################################################################################################################################