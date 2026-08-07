rm( list = ls() )

library("SimcypR")


#######  PREPARATIONS - combinations of solver-parameters to examine

Tot_iter <- 9 # {RK5, Livermore, CVODES} x {rtol = 0,1, 0.01, 0.001}

solvName <- c("RK5", "Livermore", "CVODES")
solv     <- c(1,1,1,2,2,2,3,3,3) # 1:RK5, 2:Livermore, 3:CVODES
rtol     <- c(1e-8, 1e-9, 1e-10, 1e-8, 1e-9, 1e-10, 1e-8, 1e-9, 1e-10)
atol     <- rtol * 10^6

res.df   <- data.frame(matrix(ncol = 6, nrow = Tot_iter,
                              dimnames = list(NULL, c("method","rtol","atol",
                                                      "mean_time","sd_time","subjects_negC"))))

res.df$method <- solvName[solv[1:Tot_iter]]
res.df$rtol   <- rtol
res.df$atol   <- atol

#######  SIMULATIONS - generation of databases

SimcypR::InitialiseHuman(verbose = FALSE)

path_user <- SimcypR::ScriptLocation()
setwd(path_user)

SetWorkspace("minPBPK_ADAM_CLiv.wksz") # 10 x 10 workspace

db_file  <- array( dim = Tot_iter )

for( iter in 1 : Tot_iter ){
   db_file[iter] <- paste0("sim_",sprintf("%004d",iter),".db")

   SimcypR::SetSolver( solv[iter],
                      rtol = rtol[iter], atol = atol[iter], maxsteps = 1e+07 )

   Simulate( database = db_file[iter] )
}


#######  ANALYSIS of databases - speed and accuracy (non-negativity)
library("RSQLite")

for( iter in 1 : Tot_iter ){
   conn <- RSQLite::dbConnect( SQLite(), db_file[iter] )

   # SIMULATION DURATIONS
   PopResults    <- dbGetQuery( conn, "SELECT * FROM PopResults9" )
   sim_duration  <- PopResults$SimulationDuration

   res.df$mean_time[iter] <- mean( sim_duration )
   res.df$sd_time[  iter] <-   sd( sim_duration )

   # NO. OF SUBJECTS WITH A NEGATIVE CONC
   negC <- vector()
   for(k in 1:100){
      Csys <- GetProfile_DB(ProfileID$Csys, Compound = CompoundID$Substrate, individual = k, conn, Inhibition = FALSE)
      negC[k] <- as.numeric( 0 < sum( ( Csys < 0) ) ) # there is a neg conc or not
   }

   res.df$subjects_negC[iter] <- sum( negC )

   RSQLite::dbDisconnect( conn )
}

res.df


SimcypR::Uninitialise()

