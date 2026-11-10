# Run a bespoke sensitivity analysis through SimcypR

This example uses the **SimcypR** package with the **V25 Human** engine (the workspace is a V25 workspace). 
The example is split over two files to demonstrate a bespoke sensitivity analysis.

The first file `Run Simulations.R` runs all the necessary simulations using **SimcypR** to interface with the **Simcyp Simulator** and 
the `Large molecule GSA_V25.wksz` workspace. and writes the simulation results and other ancillary files to csv. 

The second file `Run GSA.R` reads in saved data and produces results from the various uncertainty and sensitivity analysis routines. 
As all the simulations have already been run and stored in the folders of this directory,  this script can be run independently, calling these saved simulations.

## These scripts show

1. How to load a workspace, run a simulation, and open the `.db`. 
2. How to check that parameters are being successfully overwritten using `SetParameter` and `GetParameter`.
3. How to extract a profile using `GetProfile_DB`.
4. How to run looping operations for to generate simulations according to a given tabular design 
5. How to use the `Sensitivity` package to conduct Morris Test and Variance-based sensitivity analysis within R.
6. How to use `spacefillr' package to generate efficient designs for exploring parameter space.
7. How to use decoupling, which allows third party software products to be used in conjunction with the **Simcyp Simulator**.

## Files

| File | Role |
|---|---|
| `Run simulations.R` | Sampling R script |
| `Run GSA.R` | Analysis R script |
| `Large molecule GSA_V25.wksz` | Human large molecule workspace (V25) |

## How to run this example

1. Ensure the Simcyp Simulator **V26** and matching **SimcypR** package are installed. The V25 Simulator is also needed, because the workspace is a V25 workspace.
2. Keep `Large molecule GSA_V25.wksz` in the same folder as the script.
3  Keep the three sub folders and their contents in the same folder as the script.
3. Open `Run simulations.R` and run the script.
3. Open `Run GSA.R` and run the script.


## Notes

- The simulations for the variance-based sensitivity analysis are time consuming. The script `Run Simulations.R` shows how to generate designs, and the looping
operations that are required in practical use. This file may be easily adapted. 
- For the purposes of rapidly generating results, `Run GSA.R` will work with the simulation data that is supplied in the sub folders.   

## Running with a V26 workspace

The workspace supplied is a **V25** workspace, so `Run simulations.R` starts the V25 engine with
`SimcypR::Initialise(species = SpeciesID$Human, requestedVersion = 25, verbose = FALSE)`
and loads `Large molecule GSA_V25`, even though the V26 **SimcypR** package is used. 
`InitialiseHuman()` cannot be used for this, which is why `Initialise()` is called.

To run the example with a V26 workspace instead:

1. Open `Large molecule GSA_V25.wksz` in the **V26** Simcyp Simulator and save it again. This upgrades the workspace to V26.
   Save it with a new name, for example `Large molecule GSA_V26.wksz`, in the same folder as the script.
2. In `Run simulations.R`, change the workspace name to `name.file = c("Large molecule GSA_V26")`.
3. In `Run simulations.R`, replace the `SimcypR::Initialise(species = SpeciesID$Human, requestedVersion = 25, verbose = FALSE)` call with
   `SimcypR::InitialiseHuman(requestedVersion = 26, verbose = FALSE)`.
4. Run the script as normal.

Simulation results from a V26 workspace may differ slightly from the supplied V25 results in the sub folders.

