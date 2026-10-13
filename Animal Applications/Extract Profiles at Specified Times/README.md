# Extract Profiles at Specified Times (Animal)

This example uses the **SimcypR** package with the **V26 Animal** engine. It loads `Default_Animal.wksz`, runs a simulation, plots the full systemic concentration–time profile, then extracts and overlays concentrations at chosen sampling times.

## This script shows

1. How to initialise SimcypR for V26 Animal (`InitialiseAnimal`).
2. How to load a workspace, run a simulation, and open the `.db`.
3. How to plot the full systemic CT profile with `PlotProfile_DB` (no manual `GetProfile_DB` extract for the line).
4. How to extract Csys at specified times with `GetProfileAtSpecifiedTimes_DB`.
5. How to overlay those points on the `PlotProfile_DB` ggplot object.

## Files

| File | Role |
|---|---|
| `Extract Profiles at Specified Times Animal.R` | Example script |
| `Default_Animal.wksz` | Dog workspace (species is defined in the workspace) |

## How to run this example

1. Ensure the Simcyp Simulator **V26** and matching **SimcypR** package are installed.
2. Keep `Default_Animal.wksz` in the same folder as the script.
3. Open `Extract Profiles at Specified Times Animal.R` and run the script.

## Notes

- On V26, use `InitialiseAnimal(requestedVersion = 26)`. Do not pass `SpeciesID$Monkey` (or other species) into `Initialise()` — the workspace sets the animal species.
- `GetProfileAtSpecifiedTimes_DB` returns a **numeric vector** of concentrations (nearest simulated time if an exact time is missing).
