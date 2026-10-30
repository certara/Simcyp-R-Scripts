# Simcyp R Example Scripts

> [!IMPORTANT]
> **You are viewing the `main` branch, which is the maintenance branch.** It is updated continuously, so the example scripts here may be ahead of the Simcyp Simulator version you are using.
>
> To get the example scripts that match a specific Simcyp release, download that release from the [Releases page](https://github.com/certara/Simcyp-R-Scripts/releases) instead:
>
> | Release | Download |
> |---------|----------|
> | V25 | [Simcyp-R-Scripts v25](https://github.com/certara/Simcyp-R-Scripts/releases/tag/v25) |
> | V24 | [Simcyp-R-Scripts v24](https://github.com/certara/Simcyp-R-Scripts/releases/tag/v24) |
> | V23 | [Simcyp-R-Scripts v23](https://github.com/certara/Simcyp-R-Scripts/releases/tag/v23) |

## Description

**SimcypR** is an R package that supports R script interactions with the Simcyp Simulator (V25-V26). It provides functions to initialise Simcyp, load, modify and save workspaces, run simulations and interrogate results.

New for Version 26, **SimcypResults** is a companion R package for post-simulation reporting and analysis. It provides tools for PBPK model QC, concentration-time plots, PK tables, forest plots, observed-versus-predicted assessments and Word-ready table formatting. It does not drive the Simulator. Use SimcypR for that.

> **Note:** From Version 26 the package is named `SimcypR` (previously `Simcyp`). Update `library(Simcyp)` and `Simcyp::` calls in existing scripts to `library(SimcypR)` and `SimcypR::`.

## System Requirements

* Windows 10 or Windows 11 (64-bit)
* A licensed Simcyp Simulator, Version 25 or Version 26, installed on the same computer. The Simcyp Animal Simulator is needed only for animal workflows.
* R 4.6 or later, 64-bit 
* RStudio (recommended)
* Internet access to CRAN, or a site-approved CRAN mirror, during installation

## Installation

SimcypR and SimcypResults are **not on CRAN**. They are distributed by Simcyp support and can be requested from simcyp.support@certara.com. Once the request is processed, a folder named *Simcyp Simulator V26 R Packages* is added to your Simcyp Members Area (https://members.certara.co.uk) under *Software Downloads*. It contains:

| File | Description |
|------|-------------|
| `SimcypR_26.0...zip` | SimcypR Windows binary package |
| `SimcypResults_...tar.gz` | SimcypResults package |
| `Install_Simcyp_Packages.R` | Interactive installer for both packages and their CRAN dependencies |

**Do not unzip the package archives.**

### Recommended: use the installation script

1. Save `Install_Simcyp_Packages.R` in the same folder as the two package archives.
2. Open it in RStudio and click **Source**.

The script checks your R version, installs the required CRAN dependencies, installs SimcypR and then SimcypResults, and confirms the installed versions. It can also offer to run the optional unit tests, which are off by default and can take a long time.

### Alternative: manual installation

Install the CRAN packages listed under `Imports` in each package's `DESCRIPTION` file, then install both archives without unzipping them (RStudio: Tools > Install Packages > Package Archive File), or run:

```r
install.packages("path/to/SimcypR_26.0...zip", repos = NULL, type = "win.binary")
install.packages("path/to/SimcypResults_...tar.gz", repos = NULL, type = "source")
```

The current released versions are SimcypR `V26.0...` and SimcypResults `...`.

## Getting started

```r
library(SimcypR)
InitialiseHuman(verbose = FALSE)      # Human Simulator, V26 and later
SetWorkspace("YourWorkspaceName.wksz")
```

* Use `InitialiseAnimal()` for the Animal Simulator (V26 and later).
* Use `Initialise()` for V25 Human and Animal engines.
* `SaveWorkspace()` saves a modified workspace and requires a V26 engine.
* Open the vignettes with `browseVignettes(package = "SimcypR")`.
* For SimcypResults, run `library(SimcypResults)` and then `browse_SimcypResults_help()`.

> **Limitation:** the Workspace API in SimcypR is for single-session use. Concurrent sessions are not supported. Use Simcyp Cloud for large-scale or parallel simulations.

## What's new in Version 26

* Connects to Simcyp Human and Animal Simulator V25 and V26
* Dedicated `InitialiseHuman()` and `InitialiseAnimal()` functions
* `SaveWorkspace()` saves workspaces modified from R (V26 only)
* `CalcStats()` for PK and numeric summary statistics, and `PlotRibbon_DB()` for ribbon plots
* Observed-versus-predicted goodness-of-fit functions (`GetObsvsPredProfileData_DB()`, `GetObsvsPredPKData_DB()`, `PlotObsvsPred()`)
* Multi-individual profile extraction (`GetProfiles_DB()`, `GetProfilesAtSpecifiedTimes_DB()`)
* `GetPKParameters_DB()` replaces `GetAUCFrom_DB()`
* Experimental Shiny applications for the Human Simulator, launched with `RunShinyAppDashboard()`
* The new SimcypResults reporting package

See the Release Notes supplied with the packages for the full list of new and updated functions.

> The Shiny applications are experimental research tools. Check their outputs independently before use.

## Examples

Navigate to one of the example folders within the repository to view usage details.

Example scripts for SimcypResults are contained within the package itself and can be viewed by running `browse_SimcypResults_help()`. If you would like to run those examples exactly as written with the example simulation files listed, please email simcyp.support@certara.com and they will share the example Simcyp workspaces and results files in your Members Area account.

## Support

Please email simcyp.support@certara.com for any queries.

## Learning material:
If you are a Simcyp license members, E-Learning videos can be viewed via your members area.

Other learning videos can be found on YouTube:

* [Simcyp R package overview](https://youtu.be/9KuWlyrPMko?si=nuXK41YlvI4IWx28)
* [Simcyp R package overview v21](https://youtu.be/htlC2EvDWTk?si=NEKdw5i6u3oCI43S)
* [Summary statistics and ribbon plots](https://youtu.be/qn5YVLGPis4?si=7N85NKH8S229tITR)
* [Forest plots](https://youtu.be/bi6M6IpCYdc?si=gzMcdbEolDqc7hpl)
* [SimRFlow](https://www.youtube.com/watch?v=TJxSOk37EXI)

## Publications:
* [SimRFlow: An R-based workflow for automated high-throughput PBPK simulation with the Simcyp® simulator](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2022.929200/full)
* [Using the Simcyp R Package for PBPK Simulation Workflows With the Simcyp Simulator](https://pubmed.ncbi.nlm.nih.gov/40179011/)
