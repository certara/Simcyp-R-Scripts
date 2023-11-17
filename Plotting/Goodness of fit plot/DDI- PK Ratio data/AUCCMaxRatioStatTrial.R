############################################################################################################

# Clear memory: 
# rm(list=ls())

# Load Required Packages 
library("Simcyp")
library("RSQLite")
library(PKNCA)  # Compute standard Non-Compartmental Analysis (NCA) parameters for typical pharmacokinetic analyses and summarize them.

Skewness <- function (Values){
  data<- Values
  n<- length(data)
  mean<- mean(data)
  sd<- sd(data)
  Skewness<- (n*sum((data-mean)^3))/( ( (n-1)*(n-2) )*sd^3 )
  return(Skewness)
}


AUCCMaxRatioStat <- function ( compound, Alpha, Upper, Lower,  conn, Last_Dose, Trials) 
  # automatically extracts and calculates AUC/CMAx ratios
{
  
  nPop <- GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate)
  nTrials <- GetParameter(SimulationParameterID$Group , CategoryID$SimulationData, CompoundID$Substrate)
  
  quantiles_Lower <- Lower * 0.01  
  quantiles_Upper <- Upper * 0.01 
  level<- 1-(Alpha/2) 
  
  cmax_trialvalues <- array(data.frame(), dim = nTrials)
  auc_trialvalues <- array(data.frame(), dim = nTrials)
  
  # Vector to store outputs for each trial 
  PopValuesAUC<- vector()
  PopValuesCmax<- vector()
  
  for (i in 1:nPop){
     nTrial <- GetIndividualValue_DB(i, "idGroupNo", conn)
     AUC_deets<- Simcyp::GetAUCFrom_DB(ProfileID$CsysFull ,CompoundID$Substrate,individual = i,conn, allDoses = TRUE)
     #AUC_deets<- Simcyp::GetAUCFrom_DB(ProfileID$CsysFull,CompoundID$Substrate,individual = 1,conn, allDoses = TRUE)
     
     # SPLIT by inhibition group 
     auc <- split( AUC_deets , f = AUC_deets$Inhibition)
     # substrate pk WITHOUT the presence of inhibitor
     no_inh <- data.frame(auc[["0"]])
     # substrate pk WITH the presence of inhibitor
     inhibition <- data.frame(auc[["1"]])
     
    if (Last_Dose==FALSE){  # Overall dose 
      no_inh <- head(no_inh, 1)
      inhibition <- head(inhibition, 1)
    }
     
     if (Last_Dose==TRUE){  # Last dose 
       no_inh <- tail(no_inh, 1)
       inhibition <- tail(inhibition, 1)
       
     }
     
     # Ratio calculation
     cmax_ratio <- inhibition[1,"Cmax"]/no_inh[1,"Cmax"]
     auc_ratio <- inhibition[1,"AUC"]/no_inh[1,"AUC"]
     
     #Result for each trial
     cmax_trialvalues[[nTrial]] <- c(cmax_trialvalues[[nTrial]], cmax_ratio)
     auc_trialvalues[[nTrial]] <- c(auc_trialvalues[[nTrial]], auc_ratio)
     
     # Full Pop result
     PopValuesCmax<- c(PopValuesCmax, cmax_ratio)  
     PopValuesAUC<- c(PopValuesAUC, auc_ratio)  
     

  }
  
  if (Trials==TRUE){ 
    message("Returning Trial Geomean stats")
  # Trial stats summary 
  cmax <- sapply(cmax_trialvalues, geomean)
  AUC <- sapply(auc_trialvalues, geomean)

  #  Geomean Stats summary
  # cmax
  Mean <- geomean(cmax)
  Median <- median(cmax)
  SD <- sd(cmax)
  SD_Lower<- Mean-SD
  SD_Upper<- Mean+SD
  CV<- SD/ Mean/100  # Expressed as a probability
  GeoMean <- geomean(cmax)
  GeoSD <- geosd(cmax)
  GeoCV <- geocv(cmax)/100   # Expressed as a probability
  ASD<-log(GeoSD) 
  error <- qt(level,df=nPop-1)*ASD/sqrt(nPop)
  CI_Lower<- GeoMean - error  # Check calculation for Geometric mean!!
  CI_Upper<- GeoMean + error
  Min <- min(cmax)
  Max <- max(cmax)
  percentiles <- quantile(cmax, c(quantiles_Lower, quantiles_Upper))
  Centile_Lower <- percentiles[[1]]
  Centile_Upper <- percentiles[[2]]
  Fold <- Max/Min
  #Skewness<-skewness(cmax)
  Skew<-Skewness(cmax) 
  CmaxRValues<- data.frame( Mean, Median, GeoMean, GeoCV,
                            CI_Lower, CI_Upper, Centile_Lower ,Centile_Upper, CV, Min, Max, Fold,SD, SD_Lower, SD_Upper, Skew)    

  
  # AUC
  Mean <- geomean(AUC)
  Median <- median(AUC)
  SD <- sd(AUC)
  SD_Lower<- Mean-SD
  SD_Upper<- Mean+SD
  CV<- SD/ Mean/100  # Expressed as a probability
  GeoMean <- geomean(AUC)
  GeoSD <- geosd(AUC)
  GeoCV <- geocv(AUC)/100   # Expressed as a probability
  ASD<-log(GeoSD) 
  error <- qt(level,df=nPop-1)*ASD/sqrt(nPop)
  CI_Lower<- GeoMean - error  # Check calculation for Geometric mean!!
  CI_Upper<- GeoMean + error
  Min <- min(AUC)
  Max <- max(AUC)
  percentiles <- quantile(AUC, c(quantiles_Lower, quantiles_Upper))
  Centile_Lower <- percentiles[[1]]
  Centile_Upper <- percentiles[[2]]
  Fold <- Max/Min
  #Skewness<-skewness(AUC)
  Skew<-Skewness(AUC) 
  
  
  AUCRValues<- data.frame( Mean, Median, GeoMean, GeoCV,
                           CI_Lower, CI_Upper, Centile_Lower ,Centile_Upper, CV, Min, Max, Fold,SD, SD_Lower, SD_Upper, Skew)    
   
  DDIRatio<-rbind(AUCRValues,CmaxRValues)
  row.names(DDIRatio)<- c("AUC_ratio_dose1", "Cmax_ratio_dose1")
  }
  
  
  
  # Population stats
  if (Trials==FALSE){ 
    message("Returning Population stats")
  #  AUC Ratio Stats summary
  Mean <- mean(PopValuesAUC)
  Median <- median(PopValuesAUC)
  SD <- sd(PopValuesAUC)
  SD_Lower<- Mean-SD
  SD_Upper<- Mean+SD
  CV<- SD/ Mean/100  # Expressed as a probability
  GeoMean <- geomean(PopValuesAUC)
  GeoSD <- geosd(PopValuesAUC)
  GeoCV <- geocv(PopValuesAUC)/100   # Expressed as a probability
  ASD<-log(GeoSD) 
  error <- qt(level,df=nPop-1)*ASD/sqrt(nPop)
  CI_Lower<- GeoMean - error
  CI_Upper<- GeoMean + error
  Min <- min(PopValuesAUC)
  Max <- max(PopValuesAUC)
  percentiles <- quantile(PopValuesAUC, c(quantiles_Lower, quantiles_Upper))
  Centile_Lower <- percentiles[[1]]
  Centile_Upper <- percentiles[[2]]
  Fold <- Max/Min
  #Skewness<-skewness(PopValuesAUC)
  Skew<-Skewness(PopValuesAUC) 
  
  AUCRValues<- data.frame( Mean, Median, GeoMean, GeoCV,
                            CI_Lower, CI_Upper, Centile_Lower ,Centile_Upper, CV, Min, Max, Fold,SD, SD_Lower, SD_Upper, Skew)    
  
  # Cmax Ratio Stats summary
  Mean <- mean(PopValuesCmax)
  Median <- median(PopValuesCmax)
  SD <- sd(PopValuesCmax)
  SD_Lower<- Mean-SD
  SD_Upper<- Mean+SD
  CV<- SD/ Mean/100  # Expressed as a probability
  GeoMean <- geomean(PopValuesCmax)
  GeoSD <- geosd(PopValuesCmax)
  GeoCV <- geocv(PopValuesCmax)/100   # Expressed as a probability
  ASD<-log(GeoSD) 
  error <- qt(level,df=nPop-1)*ASD/sqrt(nPop)
  CI_Lower<- GeoMean - error
  CI_Upper<- GeoMean + error
  Min <- min(PopValuesCmax)
  Max <- max(PopValuesCmax)
  percentiles <- quantile(PopValuesCmax, c(quantiles_Lower, quantiles_Upper))
  Centile_Lower <- percentiles[[1]]
  Centile_Upper <- percentiles[[2]]
  Fold <- Max/Min
  #Skewness<-skewness(PopValuesCmax)
  Skew<-Skewness(PopValuesCmax) 
  
  CmaxRValues<- data.frame( Mean, Median, GeoMean, GeoCV,
                            CI_Lower, CI_Upper, Centile_Lower ,Centile_Upper, CV, Min, Max, Fold,SD, SD_Lower, SD_Upper, Skew)    
  
  
  DDIRatio<-rbind(AUCRValues,CmaxRValues)
  row.names(DDIRatio)<- c("AUC_ratio_dose1", "Cmax_ratio_dose1")
  #DDIRatio<-data.frame(t(DDIRatio))
  #names(DDIRatio)[names(DDIRatio) == "X1"] <- "AUC_Ratio"
  #names(DDIRatio)[names(DDIRatio) == "X2"] <- "Cmax_Ratio"
  }
  
  return(DDIRatio)    
}