## In this exercise we demonstrate the application of two methods for global sensitivity analysis (GSA): Morris test and variance-based.
## We additionally demonstrate the use of space-filling sampling to examine the bounding behavior of the model using an efficient design.
## In GSA we focus on several outputs from a large molecule workspace that cannot be analysed using GSA in the simulator itself.
## We use the SimcypR package to interface with the Simcyp Simulator. 
## All analysis is conducted using R packages: the Simcyp Simulator is used solely as an ODE-solver within this exercise.

#This is the first of 2 scripts and deals entirely with running simulations required for SA.

#####################################################################################################################################

## 1. Initialization
## Here we read in the R packages for analysis and read in the Simcyp workspace to be studied in GSA. 
## We run an initial simulation to check that we can run a simulation, write a database, and make changes to model parameters.

## Clear memory and remove all objects: 
rm(list=ls())

## Load Required Packages for the analysis
library("SimcypR")
library("RSQLite")
library("reshape")
library("spacefillr")
library("sensitivity")

## Initialise the version V25 engine for a human workspace
SimcypR::Initialise(species = SpeciesID$Human, requestedVersion = 25, verbose = FALSE) 

## Set script to source file location
path_user = Simcyp::ScriptLocation()
setwd(path_user)

## Import the Simcyp workspace that fully defines the PBPK model for it's baseline parameterization. 
## In the below the workspace is in the same directory as the script.
name.file = c("Large molecule GSA_V25")
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
## fu: /Workspace/Compounds/0/idFu1
## Kup adipose: /Workspace/Compounds/0/BioKupAdipose
## Kre adipose: /Workspace/Compounds/0/BioKreAdipose
## nuclease clearance: /Workspace/Compounds/0/OligosNucleaseClearance
## RISC binding: /Workspace/Compounds/0/RiscComplexBindingKd

## A mix of normal and uniform distributions is assumed here to demonstrate flexibility of approach - this readily extends.
## For normal distributions the lower and upper bounds correspond to 3 standard deviations

names = c("fu", 
          "Kup_adip", 
          "Kre_adip", 
          "CLnu",
          "KD")
id    = c("idFu1", 
          "BioKupAdipose",
          "BioKreAdipose",
          "OligosNucleaseClearance",
          "RiscComplexBindingKd")

central = c(0.25,  0.447,  0.00023386,  2.91,  0.000435)
lowB  = c(0.1, 0.0447, 0.000023386, 1.41, 0.0000435)
upB   = c(0.4,   4.47,   0.0023386,  3.41,   0.00435)
mu = c(0.2, NA, NA, 2.91, NA)
sigma = c(0.05, NA, NA, 0.5, NA)   
distribution = c("normal", "uniform", "uniform", "normal", "uniform")

# This dataframe contain parameter names, their Simcyp tags, default values, distribution parameters (where applicable), upper and lower bounds and a named probability distribution.
para.input = data.frame(names, id, central, mu, sigma, lowB, upB, distribution)
n.param = 5

## Here we use the morris function from the R sensitivity package.
## We use the decoupling feature of the package which allows separation of input and output
## Through decoupling we create a morris test object, which takes initial input and creates a  design.
## The design may executed in an external software package, and results from analysis are subsequently passed to the object.
## We execute the design in a later step of this script by running each line of the design in Simcyp and computing output.
## Morris Test indices are calculated once a vector of results has been passed back to the morris test object. 

## In the below, the command model = NULL activates decoupling.
## The argument factors specifies the number of parameters, and r is the number of trajectories.
## In the design arguments a discrete grid of 9 values is assumed for each parameter (levels argument in Simcyp). 
## The grid jump follows the levels/2 proposal of the Morris (1991) paper. Use ?morris for further details and options.
MorrisObject <- sensitivity::morris(model = NULL, factors = n.param , r = 20, design = list(type = "oat", levels = 8, grid.jump = 4))

##The design X on the unit hypercube is extracted from the object here
MorrisDesignUnit = MorrisObject$X

## Here we rescale the design on the unit hypercube to the Simcyp model parameter space.
## We use just upper and lower bounds here, ignoring information on probability distributions.
## MorrisDesignModel below is a 120 by 5 matrix of design points. Each line is a set of parameters at which the PBPK model must be run.
MorrisDesignModel = MorrisDesignUnit
for(i in 1:n.param){
  
  MorrisDesignModel[,i] = qunif(MorrisDesignUnit[,i], para.input$lowB[i] , para.input$upB[i])
}

colnames(MorrisDesignModel) = para.input$names

## Here we use the sobolSalt function from the R sensitivity package.
## We use the decoupling feature of the package which allows separation of input and output
## Through decoupling we create a sobolSalt test object, which takes initial input and creates a  design.
## The design may executed in an external software package, and results from analysis are subsequently passed to the object.
## We execute the design in a later step of this script by running each line of the design in Simcyp and computing output.
##variance-based sensitivity indices are calculated once a vector of results has been passed back to the sobolSalt test object. 

## In the below, the command model = NULL activates decoupling.
## As input the function takes two matrices of values simulated from uniform distributions and a design is created.
## The sample scheme 'A' used here is consistent with Simcyp. Bootstrapping is used to derive confidence intervals at no additional cost 

m = 1000
X1 <- data.frame(matrix(runif(n.param * m), nrow = m))
X2 <- data.frame(matrix(runif(n.param * m), nrow = m))

Sobol_object = sensitivity::sobolSalt(model = NULL, X1, X2, scheme="A", nboot = 10)

# This contains the 6000 design points for Sobol SA
SobolDesignUnit = Sobol_object$X
SobolDesignModel =  SobolDesignUnit 


## Here we rescale the design on the unit hypercube to the Simcyp model parameter space.
## In this case we use the probability distributions for rescaling.
## SobolDesignModel below is a 6000 by 5 matrix of design points. Each line is a set of parameters at which the PBPK model must be run.
for(i in 1:n.param){
  if(para.input$distribution[i] == "normal")
  {
    SobolDesignModel[,i] = qnorm(SobolDesignUnit[,i], para.input$mu[i] , para.input$sigma[i])
  }
  
  if(para.input$distribution[i] == "uniform")
  {
    SobolDesignModel[,i] = qunif(SobolDesignUnit[,i], para.input$lowB[i] , para.input$upB[i])
  }
}

colnames(SobolDesignModel) = para.input$names



## Here we use the generate_sobol_owen_set from the spacefillr package to generate a space filling design to efficiency sample parameter space
## This function takes two arguments, the number of points and the number of parameters.
## The resulting design is on the unit hypercube
n.sobol=2^7
design_unit = spacefillr::generate_sobol_owen_set(n.sobol, n.param)  # AO: do i apply the Kidney:Liver Ratio here or below in the next step (for loop)?!?


## Here we rescale the design on the unit hypercube to the Simcyp model parameter space.
## We use just upper and lower bounds here, ignoring information on probability distributions.
SpaceFillingDesign = design_unit
for(i in 1:n.param){
  SpaceFillingDesign[,i] = qunif(design_unit[,i], para.input$lowB[i] , para.input$upB[i])
}

colnames(SpaceFillingDesign) = para.input$names

#####################################################################################################################################


## 3. Run the PBPK models and collate output
## Here we read each line from the various designs, run the model and extract and store output for subsequent analysis.
## Running models and storing output follows the same looping approach for each of the three analyses.
## In this example we extract the vascular, interstitial, intracellular and endothelial concentrations in the kidney.
## The profiles for conjugated and unconjugated concentrations are available in output - we store the total concentrations (the sum) in this script
## Profile references may in general be found for GetProfileList_DB

## Beginning with simulations for Morris Test. We store the entire profiles

## Design point 1
## Overwrite default parameters. In this loop we overwrite each of the five parameters using the first line of the design.

for (j in 1:n.param){
  SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(MorrisDesignModel[1,j]))
}

## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
DBfilepath <- tempfile()

## Runs simulation and connect to database
Simulate(database = DBfilepath)
conn = RSQLite::dbConnect(SQLite(), DBfilepath)

#Output times for the simulation.
Time    = SimcypR::GetProfile_DB(ProfileID$nTimeSub, compound=-1 , individual = 1, conn, Inhibition = FALSE)

## Vascular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Vascular<- Unconjugated+Conjugated

## Endothelial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Endothelial <- Unconjugated+Conjugated

## Interstitial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Interstitial <- Unconjugated+Conjugated

# Intracellular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Intracellular <- Unconjugated+Conjugated

## Initialise dataframes for each of these four outputs and store the full concentration profiles for each output
Vascular_Kid <- data.frame(Vascular)
Endothelial_Kid<- data.frame(Endothelial)
Interstitial_Kid<- data.frame(Interstitial)
Intracellular_Kid<- data.frame(Intracellular)

RSQLite::dbDisconnect(conn)

unlink(DBfilepath)

# We loop through the remaining design points.
# In each case we overwrite parameters, run the simulation, extract outputs and append the profiles to the output dataframes.

points = dim(MorrisDesignModel)[1]
for(i in 2:points){
  
  for (j in 1:n.param){
    SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(MorrisDesignModel[i,j]))
  }
  
  ## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
  DBfilepath <- tempfile()
  
  ## Runs simulation and connect to database
  Simulate(database = DBfilepath)
  conn = RSQLite::dbConnect(SQLite(), DBfilepath)
  
  ## Vascular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Vascular<- Unconjugated+Conjugated
  
  ## Endothelial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
  Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Endothelial <- Unconjugated+Conjugated
  
  ## Interstitial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
  Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Interstitial <- Unconjugated+Conjugated
  
  # Intracellular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Intracellular <- Unconjugated+Conjugated
  
  ## We use cbind here to append the next profile to the results database.
  Vascular_Kid <- cbind(Vascular_Kid, Vascular)
  Endothelial_Kid<- cbind(Endothelial_Kid, Endothelial)
  Interstitial_Kid<- cbind(Interstitial_Kid, Interstitial)
  Intracellular_Kid<- cbind(Intracellular_Kid, Intracellular)
  
  RSQLite::dbDisconnect(conn)
  unlink(DBfilepath)
  
  #Keep track of the number of simulations run
  print(i)
  
}

## Entire set of simulations for Morris Test is complete and stored in four dataframes.
## The design and outputs are all written. The output files are ready for an analysis that may be repeated.
## These write statements assume that the folder Morris already exists in this directory.

write.csv(MorrisDesignUnit, 'Morris/MorrisDesignUnit.csv')
write.csv(MorrisDesignModel, 'Morris/MorrisDesignModel.csv')
write.csv(Vascular_Kid, 'Morris/Vascular_Kid_Morris.csv')
write.csv(Endothelial_Kid, 'Morris/Endothelial_Kid_Morris.csv')
write.csv(Interstitial_Kid, 'Morris/Interstitial_Kid_Morris.csv')
write.csv(Intracellular_Kid, 'Morris/Intracellular_Kid_Morris.csv')
write.csv(Time, 'Morris/Time.csv')


## Moving to simulations for variance based analysis (Sobol method). We store the entire profiles and then post process. 
## The process is identical to running the earlier simulations, differing only in the design points we loop thorough

## Design point 1
## Overwrite default parameters. In this loop we overwrite each of the five parameters using the first line of the design.

for (j in 1:n.param){
  SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(SobolDesignModel[1,j]))
}

## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
DBfilepath <- tempfile()

## Runs simulation and connect to database
Simulate(database = DBfilepath)
conn = RSQLite::dbConnect(SQLite(), DBfilepath)

#Output times for the simulation.
Time    = SimcypR::GetProfile_DB(ProfileID$nTimeSub, compound=-1 , individual = 1, conn, Inhibition = FALSE)

## Vascular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Vascular<- Unconjugated+Conjugated

## Endothelial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Endothelial <- Unconjugated+Conjugated

## Interstitial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Interstitial <- Unconjugated+Conjugated

# Intracellular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Intracellular <- Unconjugated+Conjugated

## Initialise dataframes for each of these four outputs and store the full concentration profiles for each output
Vascular_Kid <- data.frame(Vascular)
Endothelial_Kid<- data.frame(Endothelial)
Interstitial_Kid<- data.frame(Interstitial)
Intracellular_Kid<- data.frame(Intracellular)


RSQLite::dbDisconnect(conn)

unlink(DBfilepath)

# We loop through the remaining design points.
# In each case we overwrite parameters, run the simulation, extract outputs and append the profiles to the output dataframes.

points = dim(SobolDesignModel)[1]
for(i in 2:points){
  
  for (j in 1:n.param){
    SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(SobolDesignModel[i,j]))
  }
  
  ## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
  DBfilepath <- tempfile()
  
  ## Runs simulation and connect to database
  Simulate(database = DBfilepath)
  conn = RSQLite::dbConnect(SQLite(), DBfilepath)
  
  ## Vascular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Vascular<- Unconjugated+Conjugated
  
  ## Endothelial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
  Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Endothelial <- Unconjugated+Conjugated
  
  ## Interstitial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
  Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Interstitial <- Unconjugated+Conjugated
  
  # Intracellular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Intracellular <- Unconjugated+Conjugated
  
  ## We use cbind here to append the next profile to the results database.
  Vascular_Kid <- cbind(Vascular_Kid, Vascular)
  Endothelial_Kid<- cbind(Endothelial_Kid, Endothelial)
  Interstitial_Kid<- cbind(Interstitial_Kid, Interstitial)
  Intracellular_Kid<- cbind(Intracellular_Kid, Intracellular)
  
  RSQLite::dbDisconnect(conn)
  unlink(DBfilepath)
  
  #Keep track of the number of simulations run
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
  run = Vascular_Kid[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}

Results_Vascular = data.frame(AUC, Cmax, Tmax)

#Initialize empty vectors for results
AUC = c()
Cmax = c()
Tmax = c()

for(i in 1:points)
{
  run = Endothelial_Kid[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}

Results_Endothelial = data.frame(AUC, Cmax, Tmax)


#Initialize empty vectors for results
AUC = c()
Cmax = c()
Tmax = c()

## We consider just Sobol indices applied to profiles_Interstitial here, but the other outputs are trivial adaptions.
## We loop through and pull out AUC using the flux package and calculate a Cmax and Tmax for each point using simple calculations.
## We append each new point to our existing vector of results.
for(i in 1:points)
{
  run = Interstitial_Kid[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}

Results_Interstitial = data.frame(AUC, Cmax,Tmax)


#Initialize empty vectors for results
AUC = c()
Cmax = c()
Tmax = c()

for(i in 1:points)
{
  run = Intracellular_Kid[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}

Results_Intracellular = data.frame(AUC, Cmax,Tmax)



## The design and outputs are all written. The output files are ready for an analysis that may be repeated.
## These write statements assume that the folder Sobol already exists in this directory.


write.csv(SobolDesignUnit, 'Sobol/SobolDesignUnit.csv')
write.csv(SobolDesignModel, 'Sobol/SobolDesignModel.csv')

write.csv(Results_Interstitial, 'Sobol/Results_Interstitial.csv')
write.csv(Results_Intracellular, 'Sobol/Results_Intracellular.csv')
write.csv(Results_Vascular, 'Sobol/Results_Vascular.csv')
write.csv(Results_Endothelial, 'Sobol/Results_Endothelial.csv')


write.csv(Time, 'Sobol/Time.csv')

#We additionally store the matrices that were used as input for building the Sobol design.
SobolD_X1 = Sobol_object$X1
SobolD_X2 = Sobol_object$X2

write.csv(SobolD_X1, 'Sobol/SobolD_X1.csv')
write.csv(SobolD_X2, 'Sobol/SobolD_X2.csv')





## Finally we move simulations for a space filling design. We store again the entire profiles. 
## The process is identical to running previous simulations, differing only in the design points we loop thorough

## Design point 1
## Overwrite default parameters. In this loop we overwrite each of the five parameters using the first line of the design.

for (j in 1:n.param){
  SimcypR::SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(SpaceFillingDesign[1,j]))
}

## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
DBfilepath <- tempfile()

## Runs simulation and connect to database
Simulate(database = DBfilepath)
conn = RSQLite::dbConnect(SQLite(), DBfilepath)

#Output times for the simulation.
Time    = SimcypR::GetProfile_DB(ProfileID$nTimeSub, compound=-1 , individual = 1, conn, Inhibition = FALSE)

## Vascular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Vascular<- Unconjugated+Conjugated

## Endothelial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Endothelial <- Unconjugated+Conjugated

## Interstitial-Kidney
Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Interstitial <- Unconjugated+Conjugated

# Intracellular-Kidney
Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 

## The sum of conjugated and conjugated
Intracellular <- Unconjugated+Conjugated

## Initialise dataframes for each of these four outputs and store the full concentration profiles for each output
Vascular_Kid <- data.frame(Vascular)
Endothelial_Kid<- data.frame(Endothelial)
Interstitial_Kid<- data.frame(Interstitial)
Intracellular_Kid<- data.frame(Intracellular)

RSQLite::dbDisconnect(conn)

unlink(DBfilepath)

# We loop through the remaining design points.
# In each case we overwrite parameters, run the simulation, extract outputs and append the profiles to the output dataframes.

points = dim(SpaceFillingDesign)[1]
for(i in 2:points){
  
  for (j in 1:n.param){
    SetParameter(Tag=para.input$id[j], Category=CategoryID$Compound, CompoundID$Substrate, Value=as.numeric(SpaceFillingDesign[i,j]))
  }
  
  ## We write simulations to a temporary database. For looping operations this avoids problems with sticky memory.  
  DBfilepath <- tempfile()
  
  ## Runs simulation and connect to database
  Simulate(database = DBfilepath)
  conn = RSQLite::dbConnect(SQLite(), DBfilepath)
  
  ## Vascular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7696833052713", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701128020009", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Vascular<- Unconjugated+Conjugated
  
  ## Endothelial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696849829929", individual=1, conn, Inhibition = FALSE) 
  Conjugated <-SimcypR::GetProfile_DB("7701144797225", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Endothelial <- Unconjugated+Conjugated
  
  ## Interstitial-Kidney
  Unconjugated <-SimcypR::GetProfile_DB("7696866607145", individual=1, conn, Inhibition = FALSE)
  Conjugated <-SimcypR::GetProfile_DB("7701161574441", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Interstitial <- Unconjugated+Conjugated
  
  # Intracellular-Kidney
  Unconjugated<-SimcypR::GetProfile_DB("7697185374249", individual=1, conn, Inhibition = FALSE) 
  Conjugated<-SimcypR::GetProfile_DB("7701480341545", individual=1, conn, Inhibition = FALSE) 
  
  ## The sum of conjugated and conjugated
  Intracellular <- Unconjugated+Conjugated
  
  ## We use cbind here to append the next profile to the results database.
  Vascular_Kid <- cbind(Vascular_Kid, Vascular)
  Endothelial_Kid<- cbind(Endothelial_Kid, Endothelial)
  Interstitial_Kid<- cbind(Interstitial_Kid, Interstitial)
  Intracellular_Kid<- cbind(Intracellular_Kid, Intracellular)
  
  RSQLite::dbDisconnect(conn)
  unlink(DBfilepath)
  
  #Keep track of the number of simulations run
  print(i)
  
}


write.csv(design_unit, 'Space Filling/design_unit.csv')
write.csv(SpaceFillingDesign, 'Space Filling/SpaceFillingDesign.csv')
write.csv(Vascular_Kid, 'Space Filling/Vascular_Kid_SF.csv')
write.csv(Endothelial_Kid, 'Space Filling/Endothelial_Kid_SF.csv')
write.csv(Interstitial_Kid, 'Space Filling/Interstitial_Kid_SF.csv')
write.csv(Intracellular_Kid, 'Space Filling/Intracellular_Kid_SF.csv')
write.csv(Time, 'Space Filling/Time.csv')

SpaceFillingDesign = design_unit

#####################################################################################################################################