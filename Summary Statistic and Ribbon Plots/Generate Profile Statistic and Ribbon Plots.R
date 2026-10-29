# ==================================================================
# Summary statistics and ribbon plots with SimcypR (V26)
#
# Teaching goals
#   1. Run a DDI workspace and open the .db
#   2. Extract plasma Csys for the FULL population in one call
#      (GetProfiles_DB)
#   3. Plot individual spaghetti profiles (ggplot2)
#   4. Pull population / trial summary stats (GetProfileStats_DB)
#   5. Draw ribbon plots with PlotRibbon_DB
# ==================================================================

# Clear Global Environment ----
rm(list = ls())

# 1. Load packages -------------------------------------------------
library(SimcypR)   # Simcyp Simulator interface
library(RSQLite)   # read the simulation .db
library(ggplot2)   # individual CT plot


# 2. Initialise Simcyp engine (V26 Human) --------------------------
# For a one-shot script we Initialise only. In Shiny apps, call
# Uninitialise() first so a leftover engine is cleared.
SimcypR::InitialiseHuman(requestedVersion = 26, verbose = FALSE)


# 3. Work in the folder of this script -----------------------------
setwd(SimcypR::ScriptLocation())      # work in the folder of this script


# 4. Load workspace and run simulation -----------------------------
# Workspace major version must match the engine (here: V26).
SimcypR::SetWorkspace("Default DDI.wksz")

SimcypR::Simulate(database = "NameYourFile.db")

# Open the results database. After Simulate(), SQLite may still be in
# WAL mode — checkpoint so population tables are fully readable.
conn <- RSQLite::dbConnect(RSQLite::SQLite(), "NameYourFile.db")
try(
  DBI::dbExecute(conn, "PRAGMA wal_checkpoint(TRUNCATE);"),
  silent = TRUE
)


# 5. Extract concentration–time data for the FULL population -------
#
# Prefer GetProfiles_DB(individual = -1) over looping GetProfile_DB.
#   - GetProfile_DB  → one subject per call
#   - GetProfiles_DB → multiple subjects, or -1 for the full population
#
# Time (nTimeSub) is the same plot grid for every subject, so one
# GetProfile_DB call (individual = 1) is enough.

# Shared time vector (hours) for the plot grid
time_h <- SimcypR::GetProfile_DB(
  Tag        = ProfileID$nTimeSub,
  individual = 1L,
  Compound   = -1L,
  Inhibition = FALSE,
  conn       = conn
)

# Substrate Csys WITHOUT inhibitor (all subjects)
csys_pop <- SimcypR::GetProfiles_DB(
  Tag        = ProfileID$Csys,
  individual = -1L,
  Compound   = CompoundID$Substrate,
  Inhibition = FALSE,
  conn       = conn
)

# Substrate Csys WITH inhibitor (DDI arm)
csys_inh_pop <- SimcypR::GetProfiles_DB(
  Tag        = ProfileID$Csys,
  individual = -1L,
  Compound   = CompoundID$Substrate,
  Inhibition = TRUE,
  conn       = conn
)

# Safety: same subjects, same order, in both arms
stopifnot(
  identical(csys_pop$Individual, csys_inh_pop$Individual)
)

# Expand list-columns → long data.frame for ggplot
# (one row per subject × time point)
n_time  <- length(time_h)
n_sub   <- nrow(csys_pop)
# Profile column is the second column (name is usually "Csys")
csys_col     <- names(csys_pop)[names(csys_pop) != "Individual"][1L]
csys_inh_col <- names(csys_inh_pop)[names(csys_inh_pop) != "Individual"][1L]

sim_csys <- data.frame(
  IND         = factor(rep(csys_pop$Individual, each = n_time)),
  Time_h      = rep(time_h, times = n_sub),
  CSys_mgL    = unlist(csys_pop[[csys_col]], use.names = FALSE),
  CSysInh_mgL = unlist(csys_inh_pop[[csys_inh_col]], use.names = FALSE)
)

# Optional: confirm population size from the workspace setting
# n_sub_ws <- SimcypR::GetParameter(
#   Tag         = SimulationParameterID$Pop,
#   Category    = CategoryID$SimulationData,
#   SubCategory = CompoundID$Substrate
# )


# 6. Plot individual profiles (spaghetti) --------------------------
# Note: SimcypR reports profile concentrations as mg/L even when the
# workspace display unit differs.
#
# Okabe–Ito colours (colourblind-friendly): reuse for ggplot + ribbons
col_alone <- "#0072B2"   # blue  — substrate alone
col_ddi   <- "#D55E00"   # vermillion — with inhibitor
col_alt   <- "#009E73"   # bluish green — alternate single-arm plot

colors <- c("CSys" = col_alone, "CSysInh" = col_ddi)

ggplot2::ggplot(data = sim_csys) +
  ggplot2::geom_line(
    ggplot2::aes(x = Time_h, y = CSys_mgL, group = IND, colour = "CSys"),
    alpha = 0.35
  ) +
  ggplot2::geom_line(
    ggplot2::aes(x = Time_h, y = CSysInh_mgL, group = IND, colour = "CSysInh"),
    alpha = 0.35
  ) +
  ggplot2::scale_x_continuous("Time (h)", limits = c(0, 24)) +
  ggplot2::scale_y_continuous("Csys (mg/L)") +
  ggplot2::ggtitle("Individual systemic concentration–time profile") +
  ggplot2::labs(colour = "Profile") +
  ggplot2::scale_colour_manual(values = colors) +
  ggplot2::theme_bw()


# 7. Population / trial summary statistics -------------------------
# Trial = FALSE → Profile stats are calculated for the whole population
# Trial = TRUE  → Profile stats are calculated for per trial

# Population level — substrate alone
profile_pop_summary <- SimcypR::GetProfileStats_DB(
  Tag        = ProfileID$Csys,
  Compound   = CompoundID$Substrate,
  Inhibition = FALSE,
  Alpha      = 0.05,
  Upper      = 95,
  Lower      = 5,
  Trial      = FALSE,
  conn       = conn
)
View(profile_pop_summary)

# Population level — substrate with inhibitor
profile_pop_summary_inh <- SimcypR::GetProfileStats_DB(
  Tag        = ProfileID$Csys,
  Compound   = CompoundID$Substrate,
  Inhibition = TRUE,
  Alpha      = 0.05,
  Upper      = 95,
  Lower      = 5,
  Trial      = FALSE,
  conn       = conn
)
View(profile_pop_summary_inh)

# Trial level — substrate alone
profile_trial_summary <- SimcypR::GetProfileStats_DB(
  Tag        = ProfileID$Csys,
  Compound   = CompoundID$Substrate,
  Inhibition = FALSE,
  Alpha      = 0.05,
  Upper      = 95,
  Lower      = 5,
  Trial      = TRUE,
  conn       = conn
)
View(profile_trial_summary)

# Trial level — substrate with inhibitor
profile_trial_summary_inh <- SimcypR::GetProfileStats_DB(
  Tag        = ProfileID$Csys,
  Compound   = CompoundID$Substrate,
  Inhibition = TRUE,
  Alpha      = 0.05,
  Upper      = 95,
  Lower      = 5,
  Trial      = TRUE,
  conn       = conn
)
View(profile_trial_summary_inh)


# 8. Ribbon plots (package function) -------------------------------

# Arithmetic mean + 5th–95th percentile ribbon (substrate alone)
SimcypR::PlotRibbon_DB(
  Tag         = ProfileID$Csys,
  Compound    = CompoundID$Substrate,
  Inhibition  = FALSE,
  Ribbon_type = "Percentiles",
  Lower       = 5,
  Upper       = 95,
  Mean_type   = "Arithmetic",
  Col         = col_alone,
  Legend      = "right",
  conn        = conn
)

# Geometric mean + CI ribbon (5% alpha), log y-scale
SimcypR::PlotRibbon_DB(
  Tag         = ProfileID$Csys,
  Compound    = CompoundID$Substrate,
  Inhibition  = FALSE,
  Ribbon_type = "CI",
  Alpha       = 0.05,
  Mean_type   = "Geometric",
  ylog        = TRUE,
  Col         = col_alt,
  Legend      = "right",
  conn        = conn
)

# DDI: substrate alone + with inhibitor (two colours)
SimcypR::PlotRibbon_DB(
  Tag         = ProfileID$Csys,
  Compound    = CompoundID$Substrate,
  Inhibition  = TRUE,
  Ribbon_type = "Percentiles",
  Lower       = 5,
  Upper       = 95,
  Mean_type   = "Arithmetic",
  Col         = c(col_alone, col_ddi),
  Legend      = "right",
  conn        = conn
)

# Same DDI plot with the ribbon turned off (central lines only)
SimcypR::PlotRibbon_DB(
  Tag         = ProfileID$Csys,
  Compound    = CompoundID$Substrate,
  Inhibition  = TRUE,
  Ribbon_type = "Percentiles",
  Lower       = 5,
  Upper       = 95,
  Mean_type   = "Arithmetic",
  Col         = c(col_alone, col_ddi),
  Legend      = "right",
  ShowRibbon  = FALSE,
  conn        = conn
)


# 9. Other summary-stat helpers (explore later) --------------------
# Demographic stats:  GetIndividualValueStats_DB
# Compound results:   GetCompoundResultStats_DB
# PD stats:           GetPDResultStats_DB
# PK table summary:   CalcStats() on GetPKParameters_DB output


# 10. Disconnect ---------------------------------------------------
RSQLite::dbDisconnect(conn)

# Detach the Simulator engine from this R session.
# Do this only after all ProfileID / CompoundID extracts are done —
# Uninitialise() removes those enums.
SimcypR::Uninitialise()

# END ----
