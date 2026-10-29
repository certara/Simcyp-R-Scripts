# Run efficient global sensitivity analysis (GSA) using a meta-modelling approach

This sub-directory contains the files described below. The entire directory must be downloaded. The scripts should be placed in the same directory as the various data files.


| File | Role |
|---|---|
| `Run simulations.R` | R script for running Simcyp simulations |
| `Meta model fitting.R` | R script for fitting meta models and GSA |
| `Theophylline_V25.wksz` | Human small molecule workspace |
| `Theo Design.csv` | Design points |
| `Theo Scales.csv` | Parameter scalings |
| `Theo Results.csv` | Results file |

The first script `Run Simulations.R` runs all the necessary simulations using **SimcypR** to interface with the **Simcyp Simulator**. As input the script takes `Theo Design.csv`. For each point in the design, the parameter values in the `Theophylline_V25.wksz` workspace are overwritten and a simulation is run. Three PK summary measures are extracted from the simulation and stored. The output file from this script is a dataset of three PK summary measures corresponding to the design. This script *requires a valid Simcyp license*.

The second script `Meta model fitting.R` builds meta models and runs GSA using the *Sensitivity* package. As input the script takes `Theo Design.csv`, `Theo Scales.csv` and `Theo Results.csv`. The script is self contained - all the necessary data for replicating the tutorial paper example is provided. The script can be run *without a Simcyp license*. 


## How to run this example

To run simulations in Simcyp using `Run simulations.R`

1. Ensure the Simcyp Simulator **V25** and matching **SimcypR** package are installed.
2. Keep `Theophylline_V25.wksz` and `Theo Design.csv` in the same folder as the script.
3. Open `Run simulations.R` and run the script.

To run GSA in R using `Meta model fitting.R`

1. Keep `Theo Design.csv`, `Theo Scales.csv` and `Theo Results.csv` in the same folder as the script.
2. Open `Meta model fitting.R` and run the script.


## Notes

- The `Run Simulations.R` can be easly adapted to take a different Simcyp workspace and different parameter tags. A design specific to the workspace is required, which requies trivial adaption of the `Build a design.R` script (located in parent directory).
- The `Meta model fitting.R` script can be used in conjunction with other PBPK modelling platforms. It only requires design, results and scaling files. Some additional edits to the script will be required for outputs other than the three PK measures studied in this tutorial.

