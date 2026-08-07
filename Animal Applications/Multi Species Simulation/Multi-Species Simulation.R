# Load Packages
library(SimcypR)
library(RSQLite)

# 01.Housekeeping ------
# Setting up working directory
path_user <- SimcypR::ScriptLocation()
setwd(path_user)

# Species <- c("Dog", "Monkey", "Mouse", "Other", "Rat")
# Species <- c("Other")
Species <- c("Dog", "Monkey", "Mouse", "Rat")

# Create empty dataframe to store output
PK_Profile_Combined <- data.frame()
ST_PK_Profile_Combined <- data.frame()

# Initialize the Simulator to latest release (V26) and to the Animal Species.
SimcypR::InitialiseAnimal(verbose = FALSE)

for(i in 1:length(Species)) {

folder <- paste0("/GenericAnimal/", Species[i])

# Extract All Simcyp workspace in the folder
SimcypWksz <- list.files(path = paste0(path_user, folder), pattern = "wksz")

for(j in 1:length(SimcypWksz)) {
  SimcypR::SetWorkspace(paste0(path_user,folder,"/",SimcypWksz[j]))

  # Create DB file name
  DBfilename <- paste(fs::path_ext_remove(SimcypWksz[j]), ".db", sep="")
  DBfilepath <- paste0(path_user, folder, "/", DBfilename)

  #  Run a simulation
  Simulate(database = DBfilepath)

  #  Connect to the simulated DB file
  conn <- RSQLite::dbConnect(SQLite(), DBfilepath)

  # Extract Full conc time profiles
  Csys_Full <- GetProfile_DB(ProfileID$Csys, Compound=CompoundID$Substrate, individual=1, conn)
  Time_Full <- GetProfile_DB(ProfileID$nTimeSub, Compound=-1 , Inhibition = FALSE, conn, individual = 1) # Bug - no error with compound=-1

  PK_Profile <- data.frame(WORKSPACE=SimcypWksz[j], SPECIES=Species[i], TIME = Time_Full, CSYS = Csys_Full)
  PK_Profile_Combined <- rbind(PK_Profile_Combined, PK_Profile)

  #  Extract Profiles at specified times ----
  #  Specify Sampling time (ST)
  Sampling_Time <- c(0,0.15,0.5,1,2,3,4,5,8,12,26,20,24)

  ST_Csys <- GetProfileAtSpecifiedTimes_DB(Sampling_Time, ProfileID$Csys, individual = 1, Inhibition=FALSE, CompoundID$Substrate, conn)
  Time <- ST_Csys$x
  Csys <- ST_Csys$y

  ST_PK_Profile <- data.frame(WORKSPACE=SimcypWksz[j], SPECIES=Species[i], TIME = Time, CSYS = Csys)
  ST_PK_Profile_Combined <- rbind(ST_PK_Profile_Combined, ST_PK_Profile)

  # Disconnect database file.
  RSQLite::dbDisconnect(conn)

}

}

SimcypR::Uninitialise()

View(ST_PK_Profile_Combined)
