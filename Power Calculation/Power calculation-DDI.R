# SimcypR package - power calculation (DDI) example
#
# This script will show you:
# 1. How to run a Simcyp simulation and store results to a database
# 2. How to extract substrate AUC with and without interaction
# 3. How to run a power analysis across sample sizes


# Housekeeping ----------------------------------------------------------------
rm(list = ls())
library(SimcypR)
library(RSQLite)
library(tidyverse)

# The licensed Simcyp Simulator must be installed for the requested major
# version (V26 here). You do not need the Simulator GUI open.

# Setting up working directory
path_user <- SimcypR::ScriptLocation()
if (!nzchar(path_user)) path_user <- getwd()
setwd(path_user)


# Initialise ------------------------------------------------------------------

# Initialise the Simulator. Even if you are not simulating, you need to do this
# to use SimcypR.
try(SimcypR::Uninitialise(), silent = TRUE)
SimcypR::InitialiseHuman(requestedVersion = 26, verbose = FALSE)


# Simulate --------------------------------------------------------------------

workspace_file <- "Lorezapam+Probenecid.wksz"
db_file <- "Lorezapam+Probenecid.db"

ok <- SimcypR::SetWorkspace(workspace_file)
if (!isTRUE(ok)) {
  stop("SetWorkspace failed for ", workspace_file,
       " — check the path and that the workspace major version matches V26.")
}

# Skip re-simulation if the database is already in this folder.
if (file.exists(db_file)) {
  message("Using existing database: ", db_file)
} else {
  SimcypR::Simulate(database = db_file)
}

# Make a connection with the database file using RSQLite.
# Makes sure the database is fully written before extraction.
conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_file)
try(RSQLite::dbExecute(conn, "PRAGMA wal_checkpoint(TRUNCATE);"),
    silent = TRUE)


# Workspace population settings -----------------------------------------------

n_pop <- SimcypR::GetParameter(
  Tag = SimulationParameterID$Pop,
  Category = CategoryID$SimulationData,
  Compound = CompoundID$Substrate
)
n_pop1 <- SimcypR::GetParameter(
  Tag = SimulationParameterID$Poppercent1,
  Category = CategoryID$SimulationData,
  Compound = CompoundID$Substrate
)
n_pop2 <- SimcypR::GetParameter(
  Tag = SimulationParameterID$Poppercent2,
  Category = CategoryID$SimulationData,
  Compound = CompoundID$Substrate
)
n_pop3 <- SimcypR::GetParameter(
  Tag = SimulationParameterID$Poppercent3,
  Category = CategoryID$SimulationData,
  Compound = CompoundID$Substrate
)


# AUC extraction --------------------------------------------------------------

# Option to load a pre-saved AUC file instead of extracting from the .db:
# load("MultipopulationAUCData.RData")

pk_sub <- SimcypR::GetPKParameters_DB(
  Tag = ProfileID$CsysFull,
  Compound = CompoundID$Substrate,
  individual = -1,
  conn = conn,
  allDoses = TRUE
)
pk_sub <- pk_sub[order(pk_sub$Individual), ]

# Inhibition = 0: substrate PK without inhibitor
# Inhibition = 1: substrate PK with inhibitor
# Dose = -1: overall AUC across the simulation
auc_sub_pop <- pk_sub %>%
  dplyr::filter(Dose == -1, Inhibition == 0) %>%
  dplyr::pull(AUC)

auc_sub_ddi_pop <- pk_sub %>%
  dplyr::filter(Dose == -1, Inhibition == 1) %>%
  dplyr::pull(AUC)

RSQLite::dbDisconnect(conn)


# Split mixed populations (this workspace: 1000 subjects each) ----------------

auc_sub <- auc_sub_pop
auc_sub_healthy <- auc_sub[1:1000]
auc_sub_renal_36 <- auc_sub[1001:2000]
auc_sub_renal_l3 <- auc_sub[2001:3000]

auc_sub_inh <- auc_sub_ddi_pop
auc_sub_inh_healthy <- auc_sub_inh[1:1000]
auc_sub_inh_renal_36 <- auc_sub_inh[1001:2000]
auc_sub_inh_renal_l3 <- auc_sub_inh[2001:3000]

# DDI comparison (commented):
# auc_1 <- auc_sub_renal_l3         # substrate AUC
# auc_2 <- auc_sub_inh_renal_l3     # interaction AUC

# Population comparison used below:
auc_1 <- auc_sub_healthy     # substrate HV
auc_2 <- auc_sub_renal_l3    # substrate RI <30
auc_3 <- auc_sub_renal_36    # substrate RI 30-60

# Density of population AUC values
plot(density(auc_1), col = 1, lwd = 2, main = "", xlab = "AUC",
     ylim = c(0, 7), xlim = c(0, 0.8))
lines(density(auc_2), col = 2, lwd = 2)
lines(density(auc_3), col = 3, lwd = 2)
legend("topright",
       c("AUC Sub HV", "AUC Sub RI.l3", "AUC Sub RI.36"),
       col = c(1, 2, 3), lwd = 2, lty = 1)


# Means and variances ---------------------------------------------------------

m_pop1 <- mean(auc_1)
m_pop2 <- mean(auc_2)
m_pop3 <- mean(auc_3)

v_pop1 <- var(auc_1)
v_pop2 <- var(auc_2)
v_pop3 <- var(auc_3)

sd_pop1 <- sqrt(v_pop1)
sd_pop2 <- sqrt(v_pop2)
sd_pop3 <- sqrt(v_pop3)


# Sample sizes (equal n from each population) ---------------------------------

sample_1 <- c(10, 20, 30, 40, 50)
sample_2 <- sample_1

alpha <- 0.05
c_value <- matrix(0, 1, length(sample_1))
power_pk <- matrix(0, 1, length(sample_2))

# METHOD 1 if m_pop2 > m_pop1; METHOD 2 if m_pop2 < m_pop1
if (m_pop2 >= m_pop1) {
  # Critical value of population 1 (upper tail, alpha = 0.05)
  for (i in seq_along(sample_1)) {
    c_value[1, i] <- qnorm(1 - alpha, m_pop1,
                           sqrt(v_pop1 / sample_1[i]))
  }
  for (i in seq_along(sample_2)) {
    power_pk[1, i] <- 1 - pnorm(c_value[1, i], m_pop2,
                                sqrt(v_pop2 / sample_2[i]))
  }
} else {
  # Critical value of population 1 (lower tail, alpha = 0.05)
  for (i in seq_along(sample_1)) {
    c_value[1, i] <- qnorm(alpha, m_pop1,
                           sqrt(v_pop1 / sample_1[i]))
  }
  for (i in seq_along(sample_2)) {
    power_pk[1, i] <- pnorm(c_value[1, i], m_pop2,
                            sqrt(v_pop2 / sample_2[i]))
  }
}

power <- power_pk * 100

plot(sample_1, power, type = "l", col = "red", lwd = 2,
     xlab = "Sample size", ylab = "Power(%)")
abline(h = 80, lty = 2, lwd = 3)  # 80% power


# Finishing up ----------------------------------------------------------------

SimcypR::Uninitialise()
