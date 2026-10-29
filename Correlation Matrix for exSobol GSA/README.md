READ ME - GitHub 

# Correlation matrix for exSobol GSA 

In this example, we create a correlation matrix between Sensitivity Analysis (SA) parameters, an input necessary to run the extended Sobol GSA method. We infer SA parameter correlations by running a population simulation. This way we can circumvent large numbers of experiments in a wet lab. 

We demonstrate this method by using an Alfentanil Workspace as a simple example. 

 

## This script shows 

01. How to run a Simcyp simulation and store the results to a database file.  

02. How to extract the simulated parameter values from the database file. 

03. How to visually examine correlation between them. 

04. How to create a correlation matrix, including how to make sure that the matrix is positive definite, a mathematical property that a genuine correlation matrix needs to satisfy. 

 

## How to run this example  

Download two files - “making_corMat.R” and “Alfentanil_1x100.wksz”, then open “making_corMat.R” and run the script line by line.   

 

## References 

* Dan Liu, Linzhong Li, Amin Rostami-Hodjegan, Frederic Y Bois, Masoud Jamei (2020) Considerations and Caveats when Applying Global Sensitivity Analysis Methods to Physiologically Based Pharmacokinetic Models; AAPS Journal 22, 93. 

* Nick Higham (2002) Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22, 329--343. 
