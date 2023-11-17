####################################
# Simulate and Extract DDI results
####################################

source('AUCCMaxRatioStatTrial.R')  # can output both trial and population stats
# source('AUCCMaxRatioStatPop.R')

GetForestData_DB<- function( Alpha, Upper, Lower, conn, Last_Dose, Trials ){

DDI_stats<- AUCCMaxRatioStat(  CompoundID$Substrate, Alpha, Upper, Lower, 
                               conn, Last_Dose, Trials )
 
# Update to extract from database files instead.  
Substrate<-  GetCompoundParameter(CompoundParameterID$Drug_Name, CompoundID$Substrate)  # Obtained from Simcyp Workspace
Inhibitor1<-  GetCompoundParameter(CompoundParameterID$Drug_Name, CompoundID$Inhibitor1) # Obtained from Simcyp Workspace


 
ValueRatioPop<-  data.frame(  File= fs::path_ext_remove(basename(conn@dbname)) , 
                              VictimCompound=Substrate ,  
                              PerpCompound=Inhibitor1, 
                              PKparameter= paste(c("AUCt_ratio", "Cmax_ratio")) ,
                              Mean= DDI_stats[[1]] , Median=DDI_stats[[2]] , GeoMean= DDI_stats[[3]], 
                              GeoCV= DDI_stats[[4]] , 
                              CI90_Lower= DDI_stats[[5]] , CI90_Upper=DDI_stats[[6]] , Centile5th_Lower=DDI_stats[[7]] , Centile95th_Upper=DDI_stats[[8]] ,
                              Skewness=DDI_stats[[16]] , ArithCV=DDI_stats[[9]], Min=DDI_stats[[10]] , Max= DDI_stats[[11]], Fold= DDI_stats[[12]], StdDev= DDI_stats[[13]] ) 

return(ValueRatioPop)  
    

}

