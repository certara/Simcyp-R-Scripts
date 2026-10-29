# Run efficient global sensitivity analysis (GSA) using a meta-modelling approach

This directory contains all the materials required to replicate the example workflow from the recent tutorial paper published in AAPS. The full workflow generates an Owen scrambled Sobol sequence using the *Spacefillr* package and scales it to model parameter space. The design points are run using a **Simcyp V25** workspace for Theophylline with the **SimcypR** package used to interface with the **V25 Human** engine. The final step builds meta-models and performs GSA using the saved input and output data. All analysis is condicted using *R* packages, with **Simcyp V25** utilised as an ODE solver.

For readers who are interested in testing the meta-model methodology, the data files supplied in the two sub-directories are sufficient to replicate the example in the tutorial paper. The *read me* files in subdirectories provide details on the workflow.



