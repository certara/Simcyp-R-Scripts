## In this exercise we demonstrate the application of two methods for global sensitivity analysis (GSA): Morris test and variance-based.
## We additionally demonstrate the use of space-filling sampling to examine the bounding behavior of the model using an efficient design.
## In GSA we focus on several outputs from a large molecule workspace that cannot be analysed using GSA in the simulator itself.
## We use the SimcypR package to interface with the Simcyp Simulator. 
## All analysis is conducted using R packages: the Simcyp Simulator is used solely as an ODE-solver within this exercise.

#This is the second of 2 scripts and deals entirely with UA and SA, reading in simulation data for analysis.
#####################################################################################################################################

## 1. Initialization
## Here we read in the R packages for analysis

## Clear memory and remove all objects: 
rm(list=ls())

## Load Required Packages for the analysis

library("ggplot2")
library("sensitivity")
library("reshape")
library("spacefillr")
library("SimcypR")
library("flux")




## 2. Compute results
## Here we read in the design and outputs, compute summary measures for raw profiles as required, and calculate sensitivity measures.
## We initially run through reading in the results from space filling sampling and plot out the concentration-time profiles.
## We move on to processing the simulations from Morris Test and Sobol to extract summary measures to feed into these routines.

## Set script to source file location
path_user = SimcypR::ScriptLocation()
setwd(path_user)

## Here we assume all results have been organised into sub directories.
## Output files are read in from the sub directory 'Space Filling'
profiles_Vascular = read.csv('Space Filling/Vascular_Kid_SF.csv')
profiles_Endothelial = read.csv('Space Filling/Endothelial_Kid_SF.csv')
profiles_Interstitial = read.csv('Space Filling/Interstitial_Kid_SF.csv')
profiles_Intracellular = read.csv('Space Filling/Intracellular_Kid_SF.csv')

Time = read.csv('Space Filling/Time.csv')
points = dim(profiles_Vascular)[2] - 1

##The saved files include a leading column of labels, that we remove before analysis.
Time = Time[,2]
profiles_Vascular = profiles_Vascular[,2:(points+1)]
profiles_Endothelial = profiles_Endothelial[,2:(points+1)]
profiles_Interstitial = profiles_Interstitial[,2:(points+1)]
profiles_Intracellular = profiles_Intracellular[,2:(points+1)]


## Plotting profiles and PK summary outputs
profiles_melt = melt(profiles_Vascular)
profiles_melt$time = rep(Time, points)

Drug_profiles_Vascular = ggplot(data=profiles_melt, aes(x=time, y=value, group=variable)) +
  geom_line(data=profiles_melt, aes(x=time, y=value, group=variable), colour='darkblue', alpha=0.4) +
  ylab("Plasma concentration") +
  xlab("Time (hours)")


## Plotting profiles and PK summary outputs
profiles_melt = melt(profiles_Endothelial)
profiles_melt$time = rep(Time, points)

Drug_profiles_Endothelial = ggplot(data=profiles_melt, aes(x=time, y=value, group=variable)) +
  geom_line(data=profiles_melt, aes(x=time, y=value, group=variable), colour='darkblue', alpha=0.4) +
  ylab("Plasma concentration") +
  xlab("Time (hours)")


## Plotting profiles and PK summary outputs
profiles_melt = melt(profiles_Interstitial)
profiles_melt$time = rep(Time, points)

Drug_profiles_Interstitial = ggplot(data=profiles_melt, aes(x=time, y=value, group=variable)) +
  geom_line(data=profiles_melt, aes(x=time, y=value, group=variable), colour='darkblue', alpha=0.4) +
  ylab("Plasma concentration") +
  xlab("Time (hours)")


## Plotting profiles and PK summary outputs
profiles_melt = melt(profiles_Intracellular)
profiles_melt$time = rep(Time, points)

Drug_profiles_Intracellular = ggplot(data=profiles_melt, aes(x=time, y=value, group=variable)) +
  geom_line(data=profiles_melt, aes(x=time, y=value, group=variable), colour='darkblue', alpha=0.4) +
  ylab("Plasma concentration") +
  xlab("Time (hours)")


## For combining we use the multiplot function reproduced below from the R cookbook

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(Drug_profiles_Vascular, Drug_profiles_Endothelial, Drug_profiles_Interstitial, Drug_profiles_Intracellular, cols=2)


## We move onto the Morris Test results
## The results have been written in a rectangular format as before and require a little pre-processing before results are extracted.
## The columns of the data correspond to the data from each run.  Each row contains the data from all runs of the model at a given time point.


## Here we assume all results have been organised into sub directories.
## Output files are read in from the sub directory 'Morris'
profiles_Vascular = read.csv('Morris/Vascular_Kid_Morris.csv')
profiles_Endothelial = read.csv('Morris/Endothelial_Kid_Morris.csv')
profiles_Interstitial = read.csv('Morris/Interstitial_Kid_Morris.csv')
profiles_Intracellular = read.csv('Morris/Intracellular_Kid_Morris.csv')

Time = read.csv('Morris/Time.csv')
points = dim(profiles_Vascular)[2] - 1



##The saved files include a leading column of labels, that we remove before analysis.
Time = Time[,2]
profiles_Vascular = profiles_Vascular[,2:(points+1)]
profiles_Endothelial = profiles_Endothelial[,2:(points+1)]
profiles_Interstitial = profiles_Interstitial[,2:(points+1)]
profiles_Intracellular = profiles_Intracellular[,2:(points+1)]

## Now we loop through each column of output, and extract PK summary measures of interest
## Here we use Cmax, Tmax and AUC but there is flexibility to customize as required here

#Initialize empty vectors for results
AUC = c()
Cmax = c()
Tmax = c()

## We consider just Morris Test applied to profiles_Interstitial here, but the other outputs are trivial adaptions.
## We loop through and pull out AUC using the flux package and calculate a Cmax and Tmax for each point using simple calculations.
## We append each new point to our existing vector of results.
for(i in 1:points)
{
  run = profiles_Interstitial[,i]
  AUC = c(AUC, flux::auc(Time, run))
  Cmax = c(Cmax, max(run))
  Tmax = c(Tmax, Time[which.max(run)])
}

## We now create a Morris test object so that we can compute indices.
## We use the exact same options as we used to originally create the design. 
## We will need to overwrite the design to ensure it corresponds exactly to our model output.

## The design on the unit hypercube is read in here
MorrisDesignUnit = read.csv('Morris/MorrisDesignUnit.csv')

## Format to remove the leading column of labels
n.param = dim(MorrisDesignUnit)[2] - 1
MorrisDesignUnit = MorrisDesignUnit[,2:(n.param+1)]

## We now convert the dataframe to a matrix and un-name before passing to the Morris Object
MorrisDesignUnit = as.matrix(MorrisDesignUnit)
MorrisDesignUnit = unname(MorrisDesignUnit)

## Create Morris Test object then overwrite the design of the Morris test object
MorrisObject <- sensitivity::morris(model = NULL, factors = n.param , r = 20, design = list(type = "oat", levels = 8, grid.jump = 4))
MorrisObject$X = MorrisDesignUnit 

## We now loop through here and calculate and plot the mu* and sigma measures for AUC, Cmax and Tmax
par(mfrow=c(3,1))
#Now pass each column of summary data 
tell(MorrisObject,AUC)
print(MorrisObject)
plot(MorrisObject)
title(main='AUC')

tell(MorrisObject,Cmax)
print(MorrisObject)
plot(MorrisObject)
title(main='Cmax')

tell(MorrisObject,Tmax)
print(MorrisObject)
plot(MorrisObject)
title(main='Tmax')

## The initial plotting is a little basic. In the below we redefine our original design dataframe and append the Morris test results.
## In this dataframe format we can easily use the features of ggplot and add parameter labels.

names = c("fu", 
          "Kup_adip", 
          "Kre_adip", 
          "CLnu",
          "KD")
id    = c("idFu1", 
          "BioKupAdipose",
          "BioKreAdipose",
          "OligosNucleaseClearance",
          "RiscComplexBindingKd")

central = c(0.25,  0.447,  0.00023386,  2.91,  0.000435)
lowB  = c(0.1, 0.0447, 0.000023386, 1.41, 0.0000435)
upB   = c(0.4,   4.47,   0.0023386,  3.41,   0.00435)
mu = c(0.2, NA, NA, 2.91, NA)
sigma = c(0.05, NA, NA, 0.5, NA)   
distribution = c("normal", "uniform", "uniform", "normal", "uniform")

## Results Table, initially populated by just info on our parameters.
para.results = data.frame(names, id, central, mu, sigma, lowB, upB, distribution)

## Here we begin to loop through, calculate indices and pin them to our results dataframe.
tell(MorrisObject,AUC)
mu.star <- apply(MorrisObject$ee, 2, function(x) mean(abs(x)))
sigma <- apply(MorrisObject$ee, 2, sd)

para.results$AUC_mu.star = mu.star
para.results$AUC_sigma = sigma

tell(MorrisObject,Cmax)
mu.star <- apply(MorrisObject$ee, 2, function(x) mean(abs(x)))
sigma <- apply(MorrisObject$ee, 2, sd)

para.results$Cmax_mu.star = mu.star
para.results$Cmax_sigma = sigma

tell(MorrisObject,Tmax)
mu.star <- apply(MorrisObject$ee, 2, function(x) mean(abs(x)))
sigma <- apply(MorrisObject$ee, 2, sd)

para.results$Tmax_mu.star = mu.star
para.results$Tmax_sigma = sigma

## Here we see our full table of results from Morris Test 
para.results

## We show a plot of combined results for the three measures below using multiplot.
## We plot parameter labels in lieu of points to display results.
Morris_AUC = ggplot(para.results, aes(x=AUC_mu.star,y=AUC_sigma, label=names)) +
  geom_label(position = "nudge")

Morris_Cmax = ggplot(para.results, aes(x=Cmax_mu.star,y=Cmax_sigma, label=names)) +
  geom_label(position = "nudge") +
  scale_x_continuous(limits=c(0,0.0015))

Morris_Tmax = ggplot(para.results, aes(x=Tmax_mu.star,y=Tmax_sigma, label=names)) +
  geom_label(position = "nudge")


multiplot(Morris_AUC, Morris_Cmax, Morris_Tmax, cols=1)



## We finish the script with the Sobol results
## The results have been already been post-processed and need reading in.
## We assume that results are in the sub directory Sobol


Time = read.csv('Sobol/Time.csv')
Time = Time[,2]

Results_Vascular = read.csv('Sobol/Results_Vascular.csv')
Results_Endothelial = read.csv('Sobol/Results_Endothelial.csv')
Results_Interstitial = read.csv('Sobol/Results_Interstitial.csv')
Results_Intracellular = read.csv('Sobol/Results_Intracellular.csv')

## We now create a Sobol test object so that we can compute indices.
## We use the exact same options as we used to originally crate the design. 
## We require the exact same Sobol design that corresponds to results.
## We generate a Sobol object, and overwrite the design using the design we saved in the Run simulations script. 
## This step would not be required if we were running in the same R session and the object was saved.

## The design on the unit hypercube is read in here
## We also read in and replace the values of X1 and X2 that have already been written.
SobolDesignUnit = read.csv('Sobol/SobolDesignUnit.csv')
SobolDX1 = read.csv('Sobol/SobolD_X1.csv')
SobolDX2 = read.csv('Sobol/SobolD_X2.csv')

## Format to remove the leading column of labels
n.param = dim(SobolDesignUnit)[2] - 1
SobolDesignUnit = SobolDesignUnit[,2:(n.param+1)]
SobolDX1 = SobolDX1[,2:(n.param+1)]
SobolDX2 = SobolDX2[,2:(n.param+1)]

## We now convert the dataframes to a matrices and un-name before passing to the Sobol Object
SobolDesignUnit = as.matrix(SobolDesignUnit)
SobolDesignUnit = unname(SobolDesignUnit)

SobolDX1 = as.matrix(SobolDX1)
SobolDX1 = unname(SobolDX1)

SobolDX2 = as.matrix(SobolDX2)
SobolDX2 = unname(SobolDX2)

##We now have all the elements of the design we used and can overwrite the Sobol object.

## Here we specify a Sobol object with the same options as before, and we overwrite X1, X2 and X before supplying results 
m = 1000
X1 <- data.frame(matrix(runif(n.param * m), nrow = m))
X2 <- data.frame(matrix(runif(n.param * m), nrow = m))

## Create Sobol object
Sobol_object = sensitivity::sobolSalt(model = NULL, X1, X2, scheme="A", nboot = 10)

## Overwrite the design of the Sobol object
Sobol_object$X1 = SobolDX1
Sobol_object$X2 = SobolDX2
Sobol_object$X = SobolDesignUnit


## We now loop through here and compute results
## Pass each column of summary data using the tell command. We display the full table of results in each case.
## We also extract the main effect and total effect components of results.
## As bootstrapping has been used to derive CI's the main and total effects results are dataframes.

##We write a script here for plotting main effects against total effects with a CI for indices shown in each dimension.
##Only results for the Endothelial concentrations are given below, but the other compartments require trivial editing.


##In the code below the basic results use tell() and the main and total effects dataframes may be extracted.
##The script below is a little more elaborate and formats the results for a plot of main and total effects (with CI) 
##Only the non zero parameters are explicitly plotted in the code below.

library(ggplot2)
library(dplyr)
library(tibble)

prep <- function(df, prefix) {
  df |>
    rownames_to_column("param") |>
    transmute(param,
              "{prefix}"        := original,
              "{prefix}min"     := `min. c.i.`,
              "{prefix}max"     := `max. c.i.`)
}

#Pass AUC for Endothelial
tell(Sobol_object,Results_Endothelial$AUC)
print(Sobol_object)
S = Sobol_object$S
T = Sobol_object$T

id    = c("idFu1", 
          "BioKupAdipose",
          "BioKreAdipose",
          "OligosNucleaseClearance",
          "RiscComplexBindingKd")

plot_df <- inner_join(prep(S, "x"), prep(T, "y"), by = "param")
plot_df$param = id

thr <- 0.05   # points with both effects below this are treated as the "near-zero" cluster
plot_df$cluster <- plot_df$x < thr & plot_df$y < thr

cluster_lab <- paste(plot_df$param[plot_df$cluster], collapse = "\n")

AUC_plot <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, colour = "grey40", alpha = 0.6) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0, colour = "grey40", alpha = 0.6,
                orientation = "y") +
  geom_point(size = 2.5, colour = "#119DD9", alpha = 0.8) +
  # individual labels for the distinct points
  geom_text_repel(data = subset(plot_df, !cluster), aes(label = param),
                  size = 3, nudge_y = 0.05, min.segment.length = 0) +
  # one combined label for the near-zero cluster, offset to open space
  annotate("label", x = 0.2, y = 0.01, label = cluster_lab, size = 2.5,
           hjust = 0, fill = "white") +
  annotate("segment", x = 0.2, y = 0.01, xend = 0.01, yend = 0.01,
           colour = "grey40", linewidth = 0.4,
           arrow = arrow(length = unit(1.5, "mm"))) +
  #coord_equal() +
  labs(x = "Main effect", y = "Total effect", , title = "AUC") +
  scale_x_continuous(limits=c(-0.05,1)) +
  scale_y_continuous(limits=c(-0.05,1)) +
  theme_classic(base_size = 10)

#Pass Cmax for Endothelial
tell(Sobol_object,Results_Endothelial$Cmax)
print(Sobol_object)
S = Sobol_object$S
T = Sobol_object$T

plot_df <- inner_join(prep(S, "x"), prep(T, "y"), by = "param")
plot_df$param = id

thr <- 0.05   # points with both effects below this are treated as the "near-zero" cluster
plot_df$cluster <- plot_df$x < thr & plot_df$y < thr

cluster_lab <- paste(plot_df$param[plot_df$cluster], collapse = "\n")

Cmax_plot <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, colour = "grey40", alpha = 0.6) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0, colour = "grey40", alpha = 0.6,
                orientation = "y") +
  geom_point(size = 2.5, colour = "#119DD9", alpha = 0.8) +
  # individual labels for the distinct points
  geom_text_repel(data = subset(plot_df, !cluster), aes(label = param),
                  size = 3, nudge_y = 0.05, min.segment.length = 0) +
  # one combined label for the near-zero cluster, offset to open space
  annotate("label", x = 0.2, y = 0.01, label = cluster_lab, size = 2.5,
           hjust = 0, fill = "white") +
  annotate("segment", x = 0.2, y = 0.01, xend = 0.01, yend = 0.01,
           colour = "grey40", linewidth = 0.4,
           arrow = arrow(length = unit(1.5, "mm"))) +
  #coord_equal() +
  labs(x = "Main effect", y = "Total effect", title = "Cmax") +
  scale_x_continuous(limits=c(-0.05,1)) +
  scale_y_continuous(limits=c(-0.05,1)) +
  theme_classic(base_size = 10)

#Pass Tmax for Endothelial
tell(Sobol_object,Results_Endothelial$Tmax)
print(Sobol_object)
S = Sobol_object$S
T = Sobol_object$T

plot_df <- inner_join(prep(S, "x"), prep(T, "y"), by = "param")
plot_df$param = id

thr <- 0.05   # points with both effects below this are treated as the "near-zero" cluster
plot_df$cluster <- plot_df$x < thr & plot_df$y < thr

cluster_lab <- paste(plot_df$param[plot_df$cluster], collapse = "\n")

Tmax_plot <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, colour = "grey40", alpha = 0.6) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0, colour = "grey40", alpha = 0.6,
                orientation = "y") +
  geom_point(size = 2.5, colour = "#119DD9", alpha = 0.8) +
  # individual labels for the distinct points
  geom_text_repel(data = subset(plot_df, !cluster), aes(label = param),
                  size = 3, nudge_y = 0.05, min.segment.length = 0) +
  # one combined label for the near-zero cluster, offset to open space
  annotate("label", x = 0.2, y = 0.01, label = cluster_lab, size = 2.5,
           hjust = 0, fill = "white") +
  annotate("segment", x = 0.2, y = 0.01, xend = 0.01, yend = 0.01,
           colour = "grey40", linewidth = 0.4,
           arrow = arrow(length = unit(1.5, "mm"))) +
  #coord_equal() +
  scale_x_continuous(limits=c(-0.05,1)) +
  scale_y_continuous(limits=c(-0.05,1)) +
  labs(x = "Main effect", y = "Total effect", title = "Tmax") +
  theme_classic(base_size = 10)


multiplot(AUC_plot, Cmax_plot, Tmax_plot, cols=2)

