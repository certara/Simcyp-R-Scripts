# ==================================================================
# AUC and Cmax ratio (with DDI / without DDI) for a user-defined
# time interval, using SimcypR V26
# ==================================================================

# ---- 1. Load packages -------------------------------------------
library(SimcypR)   # talks to the Simcyp Simulator
library(RSQLite)   # reads the simulation database
library(pracma)    # trapz() = area under curve

# ---- 2. Start Simcyp, run the workspace, open the database ------
setwd(SimcypR::ScriptLocation())      # work in the folder of this script

InitialiseHuman(verbose = FALSE)      # if Simcyp is not found, add the SystemFiles path
SetWorkspace("Mikus_2020_150mgCOBI_midazolam-DDI.wksz")
Simulate(database = "MikusMDZ.db")
conn <- dbConnect(SQLite(), "MikusMDZ.db")

# ---- 3. Choose the time interval (hours) ------------------------
From <- 26.5
To   <- 28.5

# ---- 4. Get the data for ALL individuals (one call each) --------
# Each result has 2 columns: Individual, and a column holding that
# individual's time (or concentration) values as a list.

# Time points
time_pop <- GetProfiles_DB(ProfileID$nTimeSubFull, individual = -1,
                           Compound = -1, Inhibition = FALSE, conn = conn)

# Substrate concentration WITHOUT the inhibitor
csys_base <- GetProfiles_DB(ProfileID$CsysFull, individual = -1,
                            Compound = CompoundID$Substrate,
                            Inhibition = FALSE, conn = conn)

# Substrate concentration WITH the inhibitor
csys_inh <- GetProfiles_DB(ProfileID$CsysFull, individual = -1,
                           Compound = CompoundID$Substrate,
                           Inhibition = TRUE, conn = conn)

# Which trial does each individual belong to?
pk_all <- GetPKParameters_DB(ProfileID$CsysFull, CompoundID$Substrate,
                             individual = -1, conn = conn)
trial_map <- unique(pk_all[, c("Individual", "Trial")])

dbDisconnect(conn)

# Safety check: all three tables must list individuals in the same order
stopifnot(all(time_pop$Individual == csys_base$Individual),
          all(time_pop$Individual == csys_inh$Individual))

# Trial number for each individual (same order as the tables above)
Trial <- trial_map$Trial[match(csys_base$Individual, trial_map$Individual)]

# ---- 5. AUC and Cmax ratio for each individual ------------------
nInd       <- nrow(csys_base)     # number of individuals
AUC_ratio  <- numeric(nInd)       # empty vectors to fill in
Cmax_ratio <- numeric(nInd)

for (i in 1:nInd) {
  
  # this individual's time and concentrations
  time      <- time_pop[[2]][[i]]
  conc_base <- csys_base[[2]][[i]]
  conc_inh  <- csys_inh[[2]][[i]]
  
  # keep only the points inside the interval
  keep <- time >= From & time <= To
  time      <- time[keep]
  conc_base <- conc_base[keep]
  conc_inh  <- conc_inh[keep]
  
  # AUC by the trapezoidal rule, then the ratio (with / without inhibitor)
  AUC_base <- trapz(time, conc_base)
  AUC_inh  <- trapz(time, conc_inh)
  AUC_ratio[i] <- AUC_inh / AUC_base
  
  # Cmax = highest concentration in the interval, then the ratio
  Cmax_ratio[i] <- max(conc_inh) / max(conc_base)
}

# ---- 6. Summary statistics with CalcStats ----------------------
# Put the individual ratios in one table (CalcStats needs a "Trial" column)
ratios <- data.frame(Individual = csys_base$Individual,
                     Trial      = Trial,
                     AUC_ratio  = AUC_ratio,
                     Cmax_ratio = Cmax_ratio)

# (a) Geometric mean of each trial.
#     Trial = TRUE gives one row per trial plus a "Population" row at the bottom.
#     The include... options are switched off because they do not apply in trial mode.
trial_stats <- CalcStats(ratios,
                         parameters      = c("AUC_ratio", "Cmax_ratio"),
                         statistic       = "geometric_mean",
                         Trial           = TRUE,
                         includeCV       = FALSE, includeSD   = FALSE,
                         includeSkewness = FALSE, includeRange = FALSE,
                         includeFold     = FALSE)

trial_stats <- trial_stats[trial_stats$Trial != "Population", ]   # keep trial rows only

# (b) Summarise the trial values: geometric mean, lowest and highest trial.
#     (needs at least 2 trials)
across_trials <- CalcStats(trial_stats,
                           parameters      = c("AUC_ratio", "Cmax_ratio"),
                           statistic       = "geometric_mean",
                           includeSkewness = FALSE)

# Pick out the rows we need by their label in the "Statistics" column
AUC_mean  <- across_trials$AUC_ratio[across_trials$Statistics == "Geomean"]
AUC_low   <- across_trials$AUC_ratio[across_trials$Statistics == "Min"]
AUC_high  <- across_trials$AUC_ratio[across_trials$Statistics == "Max"]
cmax_mean <- across_trials$Cmax_ratio[across_trials$Statistics == "Geomean"]
cmax_low  <- across_trials$Cmax_ratio[across_trials$Statistics == "Min"]
cmax_high <- across_trials$Cmax_ratio[across_trials$Statistics == "Max"]

# ---- 7. Compare with the observed study value -------------------
name              <- "Mikus et al. 2020 (B)*"
observed.cmax     <- NA       # not reported
observed.AUCratio <- 8.79

diffCmax <- sprintf("%0.2f", cmax_mean / observed.cmax)
diffAUC  <- sprintf("%0.2f", AUC_mean  / observed.AUCratio)

predicted.cmax     <- sprintf("%0.2f( %0.2f to %0.2f)", cmax_mean, cmax_low, cmax_high)
predicted.AUCratio <- sprintf("%0.2f( %0.2f to %0.2f)", AUC_mean,  AUC_low,  AUC_high)

# ---- 8. Summary table -------------------------------------------
df1 <- data.frame(name, observed.cmax, observed.AUCratio,
                  predicted.cmax, predicted.AUCratio, diffCmax, diffAUC)
colnames(df1) <- c("Study", "Cmax Ratio", "AUC Ratio", "Cmax Ratio",
                   "AUC Ratio", "Cmax Ratio", "AUC ratio")

df1