library("SimcypR")

## 1-Trial, 100-Subject simulations

SimcypR::InitialiseHuman(verbose = FALSE)
SetWorkspace("Alfentanil_1x100.wksz") # 1 Trial x 100 Subjects
Simulate( database = "tempDB.db" ) # finishes in about a minute

## Analysis of correlation b/w SA parameters - idFu1, idBP,  idka,  idUserQgut, idEnteredVss
##    corresponding SIMULATED parameters are - fuAdj, bpAdj, kaAdj, idQgutAdj,  vdAdj
##    PROBLEM - finding latter ones is not easy 

library("RSQLite")

conn <- RSQLite::dbConnect( SQLite(), "tempDB.db" )

sim_fuAdj   <- SimcypR::GetAllCompoundResults_DB( ResultID$fuAdj,   CompoundID$Substrate, conn)
sim_bpAdj   <- SimcypR::GetAllCompoundResults_DB( ResultID$bpAdj,   CompoundID$Substrate, conn)
sim_kaAdj   <- SimcypR::GetAllCompoundResults_DB( ResultID$kaAdj,   CompoundID$Substrate, conn)
sim_QgutAdj <- SimcypR::GetAllCompoundResults_DB( ResultID$QgutAdj, CompoundID$Substrate, conn)
sim_vdAdj   <- SimcypR::GetAllCompoundResults_DB( ResultID$vdAdj,   CompoundID$Substrate, conn)

RSQLite::dbDisconnect( conn )
SimcypR::Uninitialise()

DF <- data.frame(sim_bpAdj, sim_fuAdj, sim_kaAdj, sim_QgutAdj, sim_vdAdj)
pairs( DF ) # a collection of scatter plots between SA parameters
corMat <- cor( DF ) # Correlation Matrix - this is what we want, but see below

## Safeguard - In case corMat is not Positive-Definite (PD), search the nearest PD mat

library(Matrix)

EV <- eigen( corMat )$values # PD means all eigenvalues are positive
if ( ! all( (EV >= 0) ) ) {
  corMatRev <- nearPD( corMat, corr = TRUE)$mat
}
