# Run a bespoke sensitivity analysis through SimcypR

This example uses the **SimcypR** package with the **V25 Human** engine. 
The example is split over two files to demonstrate a bespoke sensitivity analysis.

The first file `Run Simulations.R` runs all the necessary simulations using **SimcypR** to interface with the **Simcyp Simulator** and 
the `Large molecule GSA_V25.wksz` workspace. and writes the simulation results and other ancillary files to csv. 

The second file `Run GSA.R` reads in saved data and produces results from the various uncertainty and sensitivity analysis routines. 
As all the simulations have already been run and stored in the folders of this directory,  this script can be run independently, calling these saved simulations.

## These scripts shows

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
| `Large molecule GSA_V25.wksz` | Human large molecule workspace |

## How to run this example

1. Ensure the Simcyp Simulator **V26** and matching **SimcypR** package are installed.
2. Keep `Large molecule GSA.wksz` in the same folder as the script.
3  Keep the three sub folders and their contents in the same folder as the script.
3. Open `Run simulations.R` and run the script.
3. Open `Run GSA.R` and run the script.


## Notes

- The simulations for the variance-based sensitivity analysis are time consuming. The script `Run Simulations.R` shows how to generate designs, and the looping
operations that are required in practical use. This file may be easily adapted. 
- For the purposes of rapidly generating results, `Run GSA.R` will work with the simulation data that is supplied in the sub folders.   

