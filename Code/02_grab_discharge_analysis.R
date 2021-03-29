##############################################################################################################
#' @title Water chemistry and discharge

#' @author Bobby Hensley email: hensley@battelleecology.org

#' @description R script which pulls water chemistry data and discharge measurements from NEON sites 

##############################################################################################################
#### load libraries ####
library(neonUtilities)
library(plyr)
library(lubridate)

#### Concentration discharge for ARIK #### 
#' Set site and date range
  siteName="ARIK"
  startDate="2016-01"
  endDate="2020-04"
#' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
#' Creates data table for each filetered sample analyte. 
#' Non-detects are replaced with half-detection limit.
#' Outliers (>2 stdev away from mean) are set to NA. 
#' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
    grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
    Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabNa$analyteConcentration)
    for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
    grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
    grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
    Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabK$analyteConcentration)
    for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
    grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
    grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
    Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabCa$analyteConcentration)
    for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
    grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
    grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
    Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabMg$analyteConcentration)
    for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
    grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
    grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
    Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabSi$analyteConcentration)
    for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
    grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
    grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
    Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabTDS$analyteConcentration)
    for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
    grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
    grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
    Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabCl$analyteConcentration)
    for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
    grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
    grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
    Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabF$analyteConcentration)
    for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
    grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
    grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
    for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
    Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabBr$analyteConcentration)
    for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
    grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
    grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
    Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabDIC$analyteConcentration)
    for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
    grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
    grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
    Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabSO4$analyteConcentration)
    for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
    grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
    grabpH<-grabpH[,c("collectDate","initialSamplepH")]
    #' pH should never be a non-detect
    Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabpH$initialSamplepH)
    for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
    grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
    grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
    Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabFe$analyteConcentration)
    for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
    grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
    grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
    Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabMn$analyteConcentration)
    for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
    grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
    grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
    Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabNO3$analyteConcentration)
    for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
    grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
    grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
    Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabNH4$analyteConcentration)
    for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
    grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
    grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
    Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabDOC$analyteConcentration)
    for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
    grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
    grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
    for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
    Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(grabTDP$analyteConcentration)
    for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
    grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
#' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
#' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
#' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-siteStats2
#Pulls L1 discharge data 
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
#' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                               h=mean(streamStage),Q=mean(totalDischarge))  
#' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  allSiteMeans<-siteStats
  #allSiteMeans<-rbind(allSiteMeans,siteStats)
#' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
#' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
#' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  
# C-Q plots and regressions #
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
#' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Na~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[1,2]<-"Na"
    regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(K~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[2,2]<-"K"
    regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Ca~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[3,2]<-"Ca"
    regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Mg~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[4,2]<-"Mg"
    regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Si~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[5,2]<-"Si"
    regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(TDS~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[6,2]<-"TDS"
    regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Cl~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[7,2]<-"Cl"
    regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(F~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[8,2]<-"F"
    regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Br~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[9,2]<-"Br"
    regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(DIC~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[10,2]<-"DIC"
    regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(SO4~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[11,2]<-"SO4"
    regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
    fit<-lm(pH~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[12,2]<-"pH"
    regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Fe~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[13,2]<-"Fe"
    regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(Mn~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[14,2]<-"Mn"
    regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(NO3~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[15,2]<-"NO3"
    regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(NH4~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[16,2]<-"NH4"
    regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(DOC~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[17,2]<-"DOC"
    regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
    fit<-lm(TDP~Q, data=logData)
    cf<-round(coef(fit),digits=2)
    rsq<-round(summary(fit)$r.squared,digits=2)
    eq<-paste0("y = ",cf[2]," x + ",cf[1])
    r2<-paste0("R2 = ",rsq)
    abline(coef(fit),lty=2)
    mtext(eq,3,line=-2)
    mtext(r2,3,line = -3)
    regValues[18,2]<-"TDP"
    regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
    regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
    regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  allRegressionData<-regValues    
  #allRegressionData<-rbind(allRegressionData,regValues)

### Concentration discharge for BIGC #### 
  #' Set site and date range
  siteName="BIGC"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filtered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for BLDE ####
  #' Set site and date range
  siteName="BLDE"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for BLUE ####
  #' Set site and date range
  siteName="BLUE"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for BLWA #### 
  #' Set site and date range
  siteName="BLWA"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldDataADCP$collectDate<-as.POSIXct(dsc_fieldDataADCP$startDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldDataADCP[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for CARI ####
  #' Set site and date range
  siteName="CARI"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for COMO ####  
  #' Set site and date range
  siteName="COMO"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for CUPE ####  
  #' Set site and date range
  siteName="CUPE"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for FLNT ####  
  #' Set site and date range
  siteName="FLNT"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldDataADCP$collectDate<-as.POSIXct(dsc_fieldDataADCP$startDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldDataADCP[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for GUIL ####  
  #' Set site and date range
  siteName="GUIL"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for HOPB ####  
  #' Set site and date range
  siteName="HOPB"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for KING ####  
  #' Set site and date range
  siteName="KING"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for LECO ####  
  #' Set site and date range
  siteName="LECO"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for LEWI ####  
  #' Set site and date range
  siteName="LEWI"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for MART ####  
  #' Set site and date range
  siteName="MART"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for MAYF ####    
  #' Set site and date range
  siteName="MAYF"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for MCDI ####   
  #' Set site and date range
  siteName="MCDI"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for MCRA ####    
  #' Set site and date range
  siteName="MCRA"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for OKSR ####    
  #' Set site and date range
  siteName="OKSR"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for POSE ####    
  #' Set site and date range
  siteName="POSE"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for PRIN ####   
  #' Set site and date range
  siteName="PRIN"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for REDB ####    
  #' Set site and date range
  siteName="REDB"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for SYCA ####    
  #' Set site and date range
  siteName="SYCA"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)

#### Concentration discharge for TECR ####    
  #' Set site and date range
  siteName="TECR"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for TOMB ####   
  #' Set site and date range
  siteName="TOMB"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldDataADCP$collectDate<-as.POSIXct(dsc_fieldDataADCP$startDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldDataADCP[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for WALK ####    
  #' Set site and date range
  siteName="WALK"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
#### Concentration discharge for WLOU ####    
  #' Set site and date range
  siteName="WLOU"
  startDate="2016-01"
  endDate="2020-04"
  #' Pulls L1 grab sample data
  grabData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteName, startdate=startDate, 
                                         enddate=endDate, package="expanded", check.size = F)
  for(i in 1:length(grabData)) {assign(names(grabData)[i], grabData[[i]])}
  swc_externalLabDataByAnalyte$startDateTime<-as.POSIXct(swc_externalLabDataByAnalyte$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  swc_externalLabDataByAnalyte<-swc_externalLabDataByAnalyte[,c("collectDate","sampleID","analyte","analyteConcentration")]
  swc_externalLabDataByAnalyte<-na.omit(swc_externalLabDataByAnalyte)
  #' Creates data table for each filetered sample analyte. 
  #' Non-detects are replaced with half-detection limit.
  #' Outliers (>2 stdev away from mean) are set to NA. 
  #' Replicate samples are averaged.
  grabNa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Na"),]
  grabNa<-grabNa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<=0){grabNa[i,3]=0.0005}}
  Q <- quantile(grabNa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNa$analyteConcentration)
  for(i in 1:nrow(grabNa)){if(grabNa[i,3]<(Q[1]-1.5*iqr)|grabNa[i,3]>(Q[2]+1.5*iqr)){grabNa[i,3]=NA}}
  grabNa<-plyr::ddply(grabNa,c("collectDate"),summarise,Na=mean(analyteConcentration)) 
  grabK<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="K"),]
  grabK<-grabK[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabK)){if(grabK[i,3]<=0){grabK[i,3]=0.0005}}
  Q <- quantile(grabK$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabK$analyteConcentration)
  for(i in 1:nrow(grabK)){if(grabK[i,3]<(Q[1]-1.5*iqr)|grabK[i,3]>(Q[2]+1.5*iqr)){grabK[i,3]=NA}}
  grabK<-plyr::ddply(grabK,c("collectDate"),summarise,K=mean(analyteConcentration))   
  grabCa<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Ca"),]
  grabCa<-grabCa[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<=0){grabCa[i,3]=0.0005}}
  Q <- quantile(grabCa$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCa$analyteConcentration)
  for(i in 1:nrow(grabCa)){if(grabCa[i,3]<(Q[1]-1.5*iqr)|grabCa[i,3]>(Q[2]+1.5*iqr)){grabCa[i,3]=NA}}
  grabCa<-plyr::ddply(grabCa,c("collectDate"),summarise,Ca=mean(analyteConcentration)) 
  grabMg<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mg"),]
  grabMg<-grabMg[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<=0){grabMg[i,3]=0.005}}
  Q <- quantile(grabMg$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMg$analyteConcentration)
  for(i in 1:nrow(grabMg)){if(grabMg[i,3]<(Q[1]-1.5*iqr)|grabMg[i,3]>(Q[2]+1.5*iqr)){grabMg[i,3]=NA}}
  grabMg<-plyr::ddply(grabMg,c("collectDate"),summarise,Mg=mean(analyteConcentration))   
  grabSi<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Si"),]
  grabSi<-grabSi[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<=0){grabSi[i,3]=0.005}}
  Q <- quantile(grabSi$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSi$analyteConcentration)
  for(i in 1:nrow(grabSi)){if(grabSi[i,3]<(Q[1]-1.5*iqr)|grabSi[i,3]>(Q[2]+1.5*iqr)){grabSi[i,3]=NA}}
  grabSi<-plyr::ddply(grabSi,c("collectDate"),summarise,Si=mean(analyteConcentration)) 
  grabTDS<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDS"),]
  grabTDS<-grabTDS[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<=0){grabTDS[i,3]=0.05}}
  Q <- quantile(grabTDS$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDS$analyteConcentration)
  for(i in 1:nrow(grabTDS)){if(grabTDS[i,3]<(Q[1]-1.5*iqr)|grabTDS[i,3]>(Q[2]+1.5*iqr)){grabTDS[i,3]=NA}}
  grabTDS<-plyr::ddply(grabTDS,c("collectDate"),summarise,TDS=mean(analyteConcentration)) 
  grabCl<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Cl"),]
  grabCl<-grabCl[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<=0){grabCl[i,3]=0.005}}
  Q <- quantile(grabCl$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabCl$analyteConcentration)
  for(i in 1:nrow(grabCl)){if(grabCl[i,3]<(Q[1]-1.5*iqr)|grabCl[i,3]>(Q[2]+1.5*iqr)){grabCl[i,3]=NA}}
  grabCl<-plyr::ddply(grabCl,c("collectDate"),summarise,Cl=mean(analyteConcentration)) 
  grabF<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="F"),]
  grabF<-grabF[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabF)){if(grabF[i,3]<=0){grabF[i,3]=0.005}}
  Q <- quantile(grabF$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabF$analyteConcentration)
  for(i in 1:nrow(grabF)){if(grabF[i,3]<(Q[1]-1.5*iqr)|grabF[i,3]>(Q[2]+1.5*iqr)){grabF[i,3]=NA}}
  grabF<-plyr::ddply(grabF,c("collectDate"),summarise,F=mean(analyteConcentration)) 
  grabBr<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Br"),]
  grabBr<-grabBr[,c("collectDate","sampleID","analyteConcentration")] 
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<=0){grabBr[i,3]=0.005}}
  Q <- quantile(grabBr$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabBr$analyteConcentration)
  for(i in 1:nrow(grabBr)){if(grabBr[i,3]<(Q[1]-1.5*iqr)|grabBr[i,3]>(Q[2]+1.5*iqr)){grabBr[i,3]=NA}}
  grabBr<-plyr::ddply(grabBr,c("collectDate"),summarise,Br=mean(analyteConcentration)) 
  grabDIC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DIC"),]
  grabDIC<-grabDIC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<=0){grabDIC[i,3]=0.0125}}
  Q <- quantile(grabDIC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDIC$analyteConcentration)
  for(i in 1:nrow(grabDIC)){if(grabDIC[i,3]<(Q[1]-1.5*iqr)|grabDIC[i,3]>(Q[2]+1.5*iqr)){grabDIC[i,3]=NA}}
  grabDIC<-plyr::ddply(grabDIC,c("collectDate"),summarise,DIC=mean(analyteConcentration))   
  grabSO4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="SO4"),]
  grabSO4<-grabSO4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<=0){grabSO4[i,3]=0.005}}
  Q <- quantile(grabSO4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabSO4$analyteConcentration)
  for(i in 1:nrow(grabSO4)){if(grabSO4[i,3]<(Q[1]-1.5*iqr)|grabSO4[i,3]>(Q[2]+1.5*iqr)){grabSO4[i,3]=NA}}
  grabSO4<-plyr::ddply(grabSO4,c("collectDate"),summarise,SO4=mean(analyteConcentration))   
  grabpH<-swc_domainLabData[(swc_domainLabData$sampleType=="ALK"),]
  grabpH<-grabpH[,c("collectDate","initialSamplepH")]
  #' pH should never be a non-detect
  Q <- quantile(grabpH$initialSamplepH, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabpH$initialSamplepH)
  for(i in 1:nrow(grabpH)){if(grabpH[i,2]<(Q[1]-1.5*iqr)|grabpH[i,2]>(Q[2]+1.5*iqr)){grabpH[i,2]=NA}}
  grabpH<-plyr::ddply(grabpH,c("collectDate"),summarise,pH=mean(initialSamplepH))  
  grabFe<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Fe"),]
  grabFe<-grabFe[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<=0){grabFe[i,3]=0.0005}}
  Q <- quantile(grabFe$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabFe$analyteConcentration)
  for(i in 1:nrow(grabFe)){if(grabFe[i,3]<(Q[1]-1.5*iqr)|grabFe[i,3]>(Q[2]+1.5*iqr)){grabFe[i,3]=NA}}
  grabFe<-plyr::ddply(grabFe,c("collectDate"),summarise,Fe=mean(analyteConcentration)) 
  grabMn<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="Mn"),]
  grabMn<-grabMn[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<=0){grabMn[i,3]=0.0005}}
  Q <- quantile(grabMn$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabMn$analyteConcentration)
  for(i in 1:nrow(grabMn)){if(grabMn[i,3]<(Q[1]-1.5*iqr)|grabMn[i,3]>(Q[2]+1.5*iqr)){grabMn[i,3]=NA}}
  grabMn<-plyr::ddply(grabMn,c("collectDate"),summarise,Mn=mean(analyteConcentration)) 
  grabNO3<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NO3+NO2 - N"),]
  grabNO3<-grabNO3[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<=0){grabNO3[i,3]=0.0135}}
  Q <- quantile(grabNO3$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNO3$analyteConcentration)
  for(i in 1:nrow(grabNO3)){if(grabNO3[i,3]<(Q[1]-1.5*iqr)|grabNO3[i,3]>(Q[2]+1.5*iqr)){grabNO3[i,3]=NA}}
  grabNO3<-plyr::ddply(grabNO3,c("collectDate"),summarise,NO3=mean(analyteConcentration)) 
  grabNH4<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="NH4 - N"),]
  grabNH4<-grabNH4[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<=0){grabNH4[i,3]=0.002}}
  Q <- quantile(grabNH4$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabNH4$analyteConcentration)
  for(i in 1:nrow(grabNH4)){if(grabNH4[i,3]<(Q[1]-1.5*iqr)|grabNH4[i,3]>(Q[2]+1.5*iqr)){grabNH4[i,3]=NA}}
  grabNH4<-plyr::ddply(grabNH4,c("collectDate"),summarise,NH4=mean(analyteConcentration)) 
  grabDOC<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="DOC"),]
  grabDOC<-grabDOC[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<=0){grabDOC[i,3]=0.05}}
  Q <- quantile(grabDOC$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabDOC$analyteConcentration)
  for(i in 1:nrow(grabDOC)){if(grabDOC[i,3]<(Q[1]-1.5*iqr)|grabDOC[i,3]>(Q[2]+1.5*iqr)){grabDOC[i,3]=NA}}
  grabDOC<-plyr::ddply(grabDOC,c("collectDate"),summarise,DOC=mean(analyteConcentration)) 
  grabTDP<-swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte=="TDP"),]
  grabTDP<-grabTDP[,c("collectDate","sampleID","analyteConcentration")]
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<=0){grabTDP[i,3]=0.0005}}
  Q <- quantile(grabTDP$analyteConcentration, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(grabTDP$analyteConcentration)
  for(i in 1:nrow(grabTDP)){if(grabTDP[i,3]<(Q[1]-1.5*iqr)|grabTDP[i,3]>(Q[2]+1.5*iqr)){grabTDP[i,3]=NA}}
  grabTDP<-plyr::ddply(grabTDP,c("collectDate"),summarise,TDP=mean(analyteConcentration)) 
  #' Remerges individual dataframes to create one wide format table
  grabAll<-merge(grabNa,grabK,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCa,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMg,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSi,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDS,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabCl,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabF,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabBr,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDIC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabSO4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabpH,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabFe,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabMn,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabNO3,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabNH4,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll<-merge(grabAll,grabDOC,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T) 
  grabAll<-merge(grabAll,grabTDP,by.x="collectDate",by.y="collectDate",all.x=T,all.y=T)
  grabAll$siteID=siteName
  #' Caclulates mean for each solutes  
  siteStats<-data.frame(matrix(ncol=20,nrow=1))
  colnames(siteStats)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP","Q")
  siteStats[1,1]=siteName
  siteStats[1,2]=mean(grabAll$Na,na.rm=T)
  siteStats[1,3]=mean(grabAll$K,na.rm=T)
  siteStats[1,4]=mean(grabAll$Ca,na.rm=T)
  siteStats[1,5]=mean(grabAll$Mg,na.rm=T)
  siteStats[1,6]=mean(grabAll$Si,na.rm=T)
  siteStats[1,7]=mean(grabAll$TDS,na.rm=T)
  siteStats[1,8]=mean(grabAll$Cl,na.rm=T)
  siteStats[1,9]=mean(grabAll$F,na.rm=T)
  siteStats[1,10]=mean(grabAll$Br,na.rm=T)
  siteStats[1,11]=mean(grabAll$DIC,na.rm=T)
  siteStats[1,12]=mean(grabAll$SO4,na.rm=T)
  siteStats[1,13]=mean(grabAll$pH,na.rm=T)
  siteStats[1,14]=mean(grabAll$Fe,na.rm=T)
  siteStats[1,15]=mean(grabAll$Mn,na.rm=T)
  siteStats[1,16]=mean(grabAll$NO3,na.rm=T)
  siteStats[1,17]=mean(grabAll$NH4,na.rm=T)
  siteStats[1,18]=mean(grabAll$DOC,na.rm=T)
  siteStats[1,19]=mean(grabAll$TDP,na.rm=T)
  #' Caclulates stdev for each solutes  
  siteStats2<-data.frame(matrix(ncol=19,nrow=1))
  colnames(siteStats2)<-c("siteID","Na","K","Ca","Mg","Si","TDS","Cl","F","Br","DIC","SO4","pH","Fe","Mn","NO3","NH4","DOC","TDP")
  siteStats2[1,1]=siteName
  siteStats2[1,2]=sd(grabAll$Na,na.rm=T)
  siteStats2[1,3]=sd(grabAll$K,na.rm=T)
  siteStats2[1,4]=sd(grabAll$Ca,na.rm=T)
  siteStats2[1,5]=sd(grabAll$Mg,na.rm=T)
  siteStats2[1,6]=sd(grabAll$Si,na.rm=T)
  siteStats2[1,7]=sd(grabAll$TDS,na.rm=T)
  siteStats2[1,8]=sd(grabAll$Cl,na.rm=T)
  siteStats2[1,9]=sd(grabAll$F,na.rm=T)
  siteStats2[1,10]=sd(grabAll$Br,na.rm=T)
  siteStats2[1,11]=sd(grabAll$DIC,na.rm=T)
  siteStats2[1,12]=sd(grabAll$SO4,na.rm=T)
  siteStats2[1,13]=sd(grabAll$pH,na.rm=T)
  siteStats2[1,14]=sd(grabAll$Fe,na.rm=T)
  siteStats2[1,15]=sd(grabAll$Mn,na.rm=T)
  siteStats2[1,16]=sd(grabAll$NO3,na.rm=T)
  siteStats2[1,17]=sd(grabAll$NH4,na.rm=T)
  siteStats2[1,18]=sd(grabAll$DOC,na.rm=T)
  siteStats2[1,19]=sd(grabAll$TDP,na.rm=T)
  allSiteStdevs<-rbind(allSiteStdevs,siteStats2)
  #' Pulls L1 discharge data
  dischargeData<-neonUtilities::loadByProduct(dpID="DP1.20048.001", site=siteName, startdate=startDate, 
                                              enddate=endDate, package="expanded", check.size = F) 
  for(i in 1:length(dischargeData)) {assign(names(dischargeData)[i], dischargeData[[i]])}
  dsc_fieldData$startDateTime<-as.POSIXct(dsc_fieldData$collectDate,format="%Y-%m-%dT%H:%M", tz="UTC")
  dischargeData<-dsc_fieldData[,c("collectDate","streamStage","totalDischarge","totalDischargeUnits")]
  for(i in 1:nrow(dischargeData)){if(dischargeData[i,4]=="cubicMetersPerSecond"){dischargeData[i,3]=dischargeData[i,3]*1000}}
  dischargeData<-dischargeData[,c("collectDate","streamStage","totalDischarge")]
  #' Averages any replicate discharge measurements
  dischargeData<-plyr::ddply(dischargeData,c("collectDate"),summarise,
                             h=mean(streamStage),Q=mean(totalDischarge))  
  #' Calculates average discharge
  siteStats[1,20]=mean(dischargeData$Q,na.rm=T)
  #allSiteMeans<-siteStats
  allSiteMeans<-rbind(allSiteMeans,siteStats)
  #' Rounds date to make grab and discharge timestamps match  
  grabAll$collectDate<-lubridate::floor_date(grabAll$collectDate,unit="day")
  dischargeData$collectDate<-lubridate::floor_date(dischargeData$collectDate,unit="day")
  #' Matches values collected on the same day
  mergedData<-merge(grabAll,dischargeData,by.x="collectDate",by.y="collectDate",all.x=T,all.y=F)
  #' Creates a new dataframe of Log transformed data for fitting linear regerssions (C-Q relations typically power functions).
  logData<-mergedData
  logData$Na<-log10(logData$Na)  
  logData$K<-log10(logData$K)
  logData$Ca<-log10(logData$Ca)  
  logData$Mg<-log10(logData$Mg)
  logData$Si<-log10(logData$Si)  
  logData$TDS<-log10(logData$TDS)
  logData$Cl<-log10(logData$Cl)  
  logData$F<-log10(logData$F)
  logData$Br<-log10(logData$Br)  
  logData$DIC<-log10(logData$DIC)
  logData$SO4<-log10(logData$SO4)  
  #` pH already a Log scale and not transformed`
  logData$Fe<-log10(logData$Fe)  
  logData$Mn<-log10(logData$Mn)
  logData$NO3<-log10(logData$NO3)  
  logData$NH4<-log10(logData$NH4)
  logData$DOC<-log10(logData$DOC)  
  logData$TDP<-log10(logData$TDP)
  logData$Q<-log10(logData$Q) 
  #' Creates an empty dataframe to be populated with fitted regression values
  regValues<-data.frame(matrix(ncol=5,nrow=18))
  colnames(regValues)<-c("siteID","solute","slope","p-value","R-squared")
  regValues$siteID=siteName
  #' Plots and fits reressions for each solute
  plot(Na~Q, data=logData, col="blue",pch=18, ylab="Log Na (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Na~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[1,2]<-"Na"
  regValues[1,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[1,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[1,5]<-round(summary(fit)$r.squared,digits=2)
  plot(K~Q, data=logData, col="blue",pch=18, ylab="Log K (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(K~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[2,2]<-"K"
  regValues[2,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[2,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[2,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Ca~Q, data=logData, col="blue",pch=18, ylab="Log Ca (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Ca~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[3,2]<-"Ca"
  regValues[3,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[3,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[3,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mg~Q, data=logData, col="blue",pch=18, ylab="Log Mg (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mg~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[4,2]<-"Mg"
  regValues[4,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[4,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[4,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Si~Q, data=logData, col="blue",pch=18, ylab="Log Si (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Si~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[5,2]<-"Si"
  regValues[5,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[5,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[5,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDS~Q, data=logData, col="blue",pch=18, ylab="Log TDS (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDS~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[6,2]<-"TDS"
  regValues[6,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[6,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[6,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Cl~Q, data=logData, col="blue",pch=18, ylab="Log Cl (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Cl~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[7,2]<-"Cl"
  regValues[7,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[7,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[7,5]<-round(summary(fit)$r.squared,digits=2)
  plot(F~Q, data=logData, col="blue",pch=18, ylab="Log F (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(F~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[8,2]<-"F"
  regValues[8,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[8,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[8,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Br~Q, data=logData, col="blue",pch=18, ylab="Log Br (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Br~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[9,2]<-"Br"
  regValues[9,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[9,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[9,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DIC~Q, data=logData, col="blue",pch=18, ylab="Log DIC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DIC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[10,2]<-"DIC"
  regValues[10,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[10,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[10,5]<-round(summary(fit)$r.squared,digits=2)
  plot(SO4~Q, data=logData, col="blue",pch=18, ylab="Log SO4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(SO4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[11,2]<-"SO4"
  regValues[11,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[11,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[11,5]<-round(summary(fit)$r.squared,digits=2)
  plot(pH~Q, data=logData, col="blue",pch=18, ylab="pH", xlab="Log Q (L/s)")
  fit<-lm(pH~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[12,2]<-"pH"
  regValues[12,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[12,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[12,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Fe~Q, data=logData, col="blue",pch=18, ylab="Log Fe (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Fe~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[13,2]<-"Fe"
  regValues[13,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[13,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[13,5]<-round(summary(fit)$r.squared,digits=2)
  plot(Mn~Q, data=logData, col="blue",pch=18, ylab="Log Mn (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(Mn~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[14,2]<-"Mn"
  regValues[14,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[14,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[14,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NO3~Q, data=logData, col="blue",pch=18, ylab="Log NO3 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NO3~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[15,2]<-"NO3"
  regValues[15,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[15,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[15,5]<-round(summary(fit)$r.squared,digits=2)
  plot(NH4~Q, data=logData, col="blue",pch=18, ylab="Log NH4 (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(NH4~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[16,2]<-"NH4"
  regValues[16,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[16,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[16,5]<-round(summary(fit)$r.squared,digits=2)
  plot(DOC~Q, data=logData, col="blue",pch=18, ylab="Log DOC (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(DOC~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[17,2]<-"DOC"
  regValues[17,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[17,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[17,5]<-round(summary(fit)$r.squared,digits=2)
  plot(TDP~Q, data=logData, col="blue",pch=18, ylab="Log TDP (ug/L)", xlab="Log Q (L/s)")
  fit<-lm(TDP~Q, data=logData)
  cf<-round(coef(fit),digits=2)
  rsq<-round(summary(fit)$r.squared,digits=2)
  eq<-paste0("y = ",cf[2]," x + ",cf[1])
  r2<-paste0("R2 = ",rsq)
  abline(coef(fit),lty=2)
  mtext(eq,3,line=-2)
  mtext(r2,3,line = -3)
  regValues[18,2]<-"TDP"
  regValues[18,3]<-round(summary(fit)$coefficients[2,1],digits=2)
  regValues[18,4]<-round(summary(fit)$coefficients[2,4],digits=2)
  regValues[18,5]<-round(summary(fit)$r.squared,digits=2)
  #allRegressionData<-regValues    
  allRegressionData<-rbind(allRegressionData,regValues)
  
 