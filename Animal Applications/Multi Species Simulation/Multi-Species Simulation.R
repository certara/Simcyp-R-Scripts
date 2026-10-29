# ==================================================================
# Multi-species simulation — Animal (SimcypR V26)
#
# Goals
#   1. Initialise the V26 Animal engine
#   2. Load one workspace (Default_Animal.wksz — Dog by default)
#   3. Simulate Dog, then LoadPopulation for Monkey / Mouse / Rat
#   4. Save each modified species workspace (+ modification log)
#   5. Overlay species profiles with PlotProfile_DB (LoadPopulation path)
#   6. Reuse Dog.db; reload the three saved workspaces, re-simulate,
#      and overlay again — same graph expected as the LoadPopulation path
# ==================================================================

# Clear Global Environment ----
rm(list = ls())

# 1. Load packages -------------------------------------------------
library(SimcypR)   # Simcyp Simulator interface
library(RSQLite)   # open each species .db for plotting
library(ggplot2)   # used by PlotProfile_DB


# 2. Work in the folder of this script -----------------------------
path_user <- SimcypR::ScriptLocation()
if (!nzchar(path_user)) path_user <- getwd()
setwd(path_user)


# 3. Initialise Simcyp engine (V26 Animal) -------------------------
# V26: one Animal engine. Switch species with LoadPopulation()
SimcypR::InitialiseAnimal(requestedVersion = 26, verbose = FALSE)


# 4. Load the default Animal workspace (Dog) -----------------------
# Workspace major version must match the engine (here: V26).
workspace_file <- "Default_Animal.wksz"

ok <- SimcypR::SetWorkspace(workspace_file)
if (!isTRUE(ok)) {
  stop(
    "SetWorkspace failed for ", workspace_file,
    " — check the path and that the workspace major version matches V26."
  )
}

# Colours — one per species
species_cols <- c(
  Dog    = "#0072B2",  # blue
  Monkey = "#D55E00",  # vermillion
  Mouse  = "#009E73",  # bluish green
  Rat    = "#E69F00"   # orange
)

# Populations to load after the workspace Dog default.
# Names come from AnimalPopulationFile (after InitialiseAnimal).
# Inspect with: names(AnimalPopulationFile)
species_pops <- list(
  Monkey = AnimalPopulationFile$sim_Monkey,
  Mouse  = AnimalPopulationFile$Sim_Mouse,
  Rat    = AnimalPopulationFile$Sim_Rat
)


# 5. Helpers -------------------------------------------------------
simulate_species_db <- function(db_file) {
  SimcypR::Simulate(
    database        = db_file,
    verbose         = FALSE,
    displayProgress = TRUE
  )

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_file)
  try(
    DBI::dbExecute(conn, "PRAGMA wal_checkpoint(TRUNCATE);"),
    silent = TRUE
  )
  conn
}

plot_csys <- function(conn, species_label, title, overlay = NULL) {
  SimcypR::PlotProfile_DB(
    Tag        = ProfileID$Csys,
    Compound   = CompoundID$Substrate,
    Inhibition = FALSE,
    conn       = conn,
    Col        = species_cols[[species_label]],
    DBLabel    = species_label,
    Overlay    = overlay,
    Title      = title,
    Xlabel     = "Time (h)",
    Ylabel     = "Concentration (mg/L)"
  )
}


# 6. Dog (workspace default) — first plot --------------------------
# In-memory LoadPopulation path (no saved .wksz yet for Dog).
conn_dog <- simulate_species_db("Dog.db")

p_loadpop <- plot_csys(
  conn          = conn_dog,
  species_label = "Dog",
  title         = "Multi Species Simulation"
)

RSQLite::dbDisconnect(conn_dog)


# 7. Monkey, Mouse, Rat — LoadPopulation, SaveWorkspace, Overlay plots ---
# Keep the same workspace loaded. LoadPopulation swaps the animal
# population. SaveWorkspace writes a new .wksz (and a modification
# log when logChanges = TRUE) into this script folder. Simulate then
# writes a separate .db per species. PlotProfile_DB(Overlay = ...)
# stacks each new profile.
saved_workspaces <- character()

for (species_label in names(species_pops)) {
  pop_path <- species_pops[[species_label]]

  loaded <- SimcypR::LoadPopulation(Population = 1, path = pop_path)
  if (!isTRUE(loaded)) {
    stop(
      "LoadPopulation failed for ", species_label,
      " — check AnimalPopulationFile after InitialiseAnimal()."
    )
  }

  # Persist the modified in-memory workspace for this species.
  # Path must differ from the currently loaded Default_Animal.wksz.
  # logChanges = TRUE also writes a modification-log CSV beside it.
  wksz_out <- file.path(path_user, paste0(species_label, ".wksz"))
  saved <- SimcypR::SaveWorkspace(path = wksz_out, logChanges = TRUE)
  if (!isTRUE(saved)) {
    stop("SaveWorkspace failed for ", wksz_out)
  }
  saved_workspaces[[species_label]] <- wksz_out
  message("Saved workspace (+ modification log): ", wksz_out)

  conn_sp <- simulate_species_db(paste0(species_label, ".db"))

  p_loadpop <- plot_csys(
    conn          = conn_sp,
    species_label = species_label,
    title         = "Multi Species Simulation",
    overlay       = p_loadpop
  )

  RSQLite::dbDisconnect(conn_sp)
}

print(p_loadpop)


# 8. Reload saved workspaces and overlay again ---------------------
# Reload the three saved species workspaces, simulate, and overlay.
# The four-species graph should match p_loadpop.
conn_dog_reuse <- RSQLite::dbConnect(RSQLite::SQLite(), "Dog.db")
try(
  DBI::dbExecute(conn_dog_reuse, "PRAGMA wal_checkpoint(TRUNCATE);"),
  silent = TRUE
)

p_saved_wksz <- plot_csys(
  conn          = conn_dog_reuse,
  species_label = "Dog",
  title         = "Multi Species Simulation- Method 2"
)

RSQLite::dbDisconnect(conn_dog_reuse)

for (species_label in names(saved_workspaces)) {
  wksz_path <- saved_workspaces[[species_label]]

  ok_ws <- SimcypR::SetWorkspace(wksz_path)
  if (!isTRUE(ok_ws)) {
    stop(
      "SetWorkspace failed for ", wksz_path,
      " — check the path and workspace vs engine version."
    )
  }

  conn_ws <- simulate_species_db(
    paste0("from_wksz_", species_label, ".db")
  )

  p_saved_wksz <- plot_csys(
    conn          = conn_ws,
    species_label = species_label,
    title         = "Multi Species Simulation- Method 2",
    overlay       = p_saved_wksz
  )

  RSQLite::dbDisconnect(conn_ws)
}

## Compare the two overlays. ----
# Print both plots — they should match (SaveWorkspace preserved the species swap).
print(p_loadpop)
print(p_saved_wksz)


# 9. Uninitialise --------------------------------------------------
SimcypR::Uninitialise()
