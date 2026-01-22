
# Load Packages
library(Simcyp)
library(RSQLite)


# 01.Housekeeping ------
# Setting up working directory
path_user <- Simcyp::ScriptLocation()
setwd(path_user)


# Initialize the Simulator to latest release and to the Animal Species.
Simcyp::Initialise( species = SpeciesID$Monkey)

# Extract All Simcyp workspace in the working folder
SimcypWksz <- list.files(pattern = "wksz")

#  Set a workspace and Run a simulation ----
# Set workspace
SetWorkspace(SimcypWksz)

# Create DB file name
DBfilename <- sub("wksz", "db", SimcypWksz)   # update file extension

#  Run a simulation
Simulate(database = DBfilename)

#  Connect to the simulated DB file
conn <- RSQLite::dbConnect(SQLite(),DBfilename)


#  Extract Profiles at specified times (ST)----
#  Specify Sampling time
Sampling_Time<-c(0,0.15,0.5,1,2,3,4,5,8,12,26,20,24)

ST_Csys<- GetProfileAtSpecifiedTimes_DB(Sampling_Time, ProfileID$Csys,individual = 1,inhibition=FALSE,CompoundID$Substrate,conn)
Time<-ST_Csys$x
Csys<-ST_Csys$y

# Check if the times match
Time_Comb<- data.frame(Sampling_Time, Time)
Time_Comb


plot(Time, Csys, lwd=2, col=2, main= "Systemic concentration-Monkey")


# Overlay full sampled time predictions ----
# Full conc time profiles
Csys_Full<-GetProfile_DB(ProfileID$Csys, compound=CompoundID$Substrate,  individual=1, conn)
Time_Full<- GetProfile_DB(ProfileID$nTimeSub, compound=-1 , inhibition = FALSE, conn, individual = 1)


# Overlay full data
lines(Time_Full, Csys_Full)


# Disconnect from the database
RSQLite::dbDisconnect(conn)

# Finishing up ------------------
Uninitialise()   # Uninitialise the Simcyp engine
