# ==================================================================
# Extract profiles at specified times — Animal
#
# Teaching goals
#   1. Initialise the V26 Animal engine
#   2. Load Default_Animal.wksz, simulate, open the .db
#   3. Plot the full systemic CT profile with PlotProfile_DB
#   4. Extract Csys at chosen sampling times using GetProfileAtSpecifiedTimes_DB
#   5. Overlay those points on the PlotProfile_DB ggplot object
# ==================================================================

# Clear Global Environment ----
rm(list = ls())

# 1. Load packages -------------------------------------------------
library(SimcypR)   # Simcyp Simulator interface
library(RSQLite)   # read the simulation .db
library(ggplot2)   # overlay sampling points on the profile plot


# 2. Work in the folder of this script -----------------------------
path_user <- SimcypR::ScriptLocation()
if (!nzchar(path_user)) path_user <- getwd()
setwd(path_user)


# 3. Initialise Simcyp engine (V26 Animal) -------------------------
SimcypR::InitialiseAnimal(requestedVersion = 26, verbose = FALSE)


# 4. Load workspace and run simulation -----------------------------
# Workspace major version must match the engine (here: V26).
# Default_Animal.wksz is a Dog workspace.
workspace_file <- "Default_Animal.wksz"
db_file <- "Default_Animal.db"

ok <- SimcypR::SetWorkspace(workspace_file)
if (!isTRUE(ok)) {
  stop(
    "SetWorkspace failed for ", workspace_file,
    " — check the path and that the workspace major version matches V26."
  )
}

SimcypR::Simulate(database = db_file)

# Open the results database. After Simulate(), SQLite may still be in
# WAL mode — checkpoint so tables are fully readable.
conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_file)
try(
  DBI::dbExecute(conn, "PRAGMA wal_checkpoint(TRUNCATE);"),
  silent = TRUE
)


# 5. Plot the full systemic concentration–time profile -------------
# PlotProfile_DB reads the .db and draws the profile. No need to
# call GetProfile_DB for the full series when you only want a plot.
#
# Here we keep the package plot, then add sampling points with ggplot2.
p_full <- SimcypR::PlotProfile_DB(
  Tag        = ProfileID$Csys,
  Compound   = CompoundID$Substrate,
  Inhibition = FALSE,
  conn       = conn,
  Title      = "Systemic concentration — Dog",
  Xlabel     = "Time (h)",
  Ylabel     = "Concentration (mg/L)"
)

print(p_full)


# 6. Extract profile at specified sampling times -------------------
# Requested sampling times (hours). If an exact time is not in the
# simulation output, V26 returns the value at the nearest time.
sampling_time <- c(0, 0.15, 0.5, 1, 2, 3, 4, 5, 8, 12, 20, 24)

# V26: returns a numeric vector of concentrations (not list(x, y)).
csys_at_st <- SimcypR::GetProfileAtSpecifiedTimes_DB(
  times      = sampling_time,
  Tag        = ProfileID$Csys,
  individual = 1L,
  Compound   = CompoundID$Substrate,
  Inhibition = FALSE,
  conn       = conn
)

sampling_df <- data.frame(
  Time_h   = sampling_time,
  CSys_mgL = csys_at_st
)

print(sampling_df)


# 7. Overlay sampling points on the full-profile plot --------------
# Add ggplot2 points onto the PlotProfile_DB object.
# inherit.aes = FALSE avoids clashing with package plot aesthetics.
p_overlay <- p_full +
  ggplot2::geom_point(
    data        = sampling_df,
    mapping     = ggplot2::aes(x = Time_h, y = CSys_mgL),
    inherit.aes = FALSE,
    colour      = "#D55E00",
    size        = 2.5
  ) +
  ggplot2::labs(
    caption = "Orange points: Csys at specified sampling times (Dog)"
  )

print(p_overlay)

# Optional alternative (package Obs argument, no extra ggplot layer):
# obs_df <- data.frame(
#   SubjectID = "Dog",
#   Time      = sampling_time,
#   DV        = csys_at_st
# )
# SimcypR::PlotProfile_DB(
#   Tag = ProfileID$Csys, Compound = CompoundID$Substrate,
#   Inhibition = FALSE, conn = conn, Obs = obs_df
# )


# 8. Disconnect and Uninitialise -----------------------------------
RSQLite::dbDisconnect(conn)

SimcypR::Uninitialise()
