# Multi-Species Simulation

This example uses the **SimcypR** package with the **V26 Animal** engine. It starts from one workspace (`Default_Animal.wksz`, Dog by default), loads Monkey, Mouse, and Rat with `LoadPopulation`, **saves each modified workspace** (with a modification log), overlays the systemic CT profiles, then reloads those saved workspaces and shows that the profiles match.

## This script shows

1. How to initialise SimcypR for V26 Animal (`InitialiseAnimal`).
2. How to load one Animal workspace (Dog default).
3. How to switch species with `LoadPopulation` and `AnimalPopulationFile$*`.
4. How to **save** each modified species workspace with `SaveWorkspace(..., logChanges = TRUE)` (writes `.wksz` + modification-log CSV in the script folder).
5. How to overlay species profiles with `PlotProfile_DB(..., Overlay = …)` and a different colour per species.
6. How to reuse the existing `Dog.db`, reload the three saved workspaces, re-simulate those species, and overlay again — same result as the in-memory path (both overlay plots printed for comparison).

## Files

| File | Role |
|---|---|
| `Multi-Species Simulation.R` | Example script |
| `Default_Animal.wksz` | Starting Animal workspace (Dog) |

**Created when you run the script** (in the script folder):

| Output | Role |
|---|---|
| `Monkey.wksz`, `Mouse.wksz`, `Rat.wksz` | Saved species workspaces |
| Matching modification-log CSV(s) | Change log from `logChanges = TRUE` |
| `Dog.db`, `Monkey.db`, `Mouse.db`, `Rat.db` | Databases from the LoadPopulation path |
| `from_wksz_Monkey.db`, `from_wksz_Mouse.db`, `from_wksz_Rat.db` | Databases from the reloaded saved workspaces (Dog reuses `Dog.db`) |

## How to run this example

1. Ensure the Simcyp Simulator **V26** and matching **SimcypR** package are installed.
2. Keep `Default_Animal.wksz` in the same folder as the script.
3. Open `Multi-Species Simulation.R` and run the script.

## Notes

- On V26, use `InitialiseAnimal(requestedVersion = 26)`. Do not pass `SpeciesID$Rat` (etc.) into `Initialise()`.
- After initialise, inspect available populations with `names(AnimalPopulationFile)`.
- This example uses:
  - Dog — workspace default (no `LoadPopulation`)
  - Monkey — `AnimalPopulationFile$sim_Monkey`
  - Mouse — `AnimalPopulationFile$Sim_Mouse`
  - Rat — `AnimalPopulationFile$Sim_Rat`
- `SaveWorkspace` path must differ from the currently loaded workspace path.
- The two overlay plots should look the same: saving then reloading the workspaces preserves the species change. (`PlotProfile_DB` is ggplot-based, so `par(mfrow)` cannot arrange the two plots in one device.)
