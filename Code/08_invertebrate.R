########################################################################################################
#' @title D08invert

#' @author Bobby Hensley email: hensley@battelleecology.org

#' @description R script which downloads and examines NEON D08 invertebrate data
#'
#' Created 8/24/2021
########################################################################################################
#' Loads required packages from library 
library(neonUtilities)
# library(tidyr)
# library(ggplot2)
# library(lubridate)
# library(plyr)
# library(streamMetabolizer)


#' Downloads invertebrate data from NEON API into R environment
invert<-loadByProduct(dpID="DP1.20120.001", site="MAYF", startdate="2014-01", enddate="2021-08", check.size = F)
list2env(invert, .GlobalEnv)
invertSummaryMAYF<-plyr::ddply(inv_taxonomyProcessed,c("order"),summarise,totalCountMAYF=sum(estimatedTotalCount)) 
invertSummaryMAYF<-invertSummaryMAYF[(invertSummaryMAYF$order!="NA"),]
boutsMAYF<-as.numeric(nrow(as.data.frame(unique(inv_taxonomyProcessed$collectDate))))
totalCollectMAYF<-sum(invertSummaryMAYF$totalCount,na.rm=T)
invertSummaryMAYF$avgCountMAYF=invertSummaryMAYF$totalCount/boutsMAYF
invertSummaryMAYF$MAYFpercent=invertSummaryMAYF$totalCount/totalCollectMAYF

invert<-loadByProduct(dpID="DP1.20120.001", site="BLWA", startdate="2014-01", enddate="2021-08", check.size = F)
list2env(invert, .GlobalEnv)
invertSummaryBLWA<-plyr::ddply(inv_taxonomyProcessed,c("order"),summarise,totalCountBLWA=sum(estimatedTotalCount)) 
invertSummaryBLWA<-invertSummaryBLWA[(invertSummaryBLWA$order!="NA"),]
boutsBLWA<-as.numeric(nrow(as.data.frame(unique(inv_taxonomyProcessed$collectDate))))
totalCollectBLWA<-sum(invertSummaryBLWA$totalCount,na.rm=T)
invertSummaryBLWA$avgCountBLWA=invertSummaryBLWA$totalCount/boutsBLWA
invertSummaryBLWA$BLWApercent=invertSummaryBLWA$totalCount/totalCollectBLWA

invert<-loadByProduct(dpID="DP1.20120.001", site="TOMB", startdate="2014-01", enddate="2021-08", check.size = F)
list2env(invert, .GlobalEnv)
invertSummaryTOMB<-plyr::ddply(inv_taxonomyProcessed,c("order"),summarise,totalCountTOMB=sum(estimatedTotalCount)) 
invertSummaryTOMB<-invertSummaryTOMB[(invertSummaryTOMB$order!="NA"),]
boutsTOMB<-as.numeric(nrow(as.data.frame(unique(inv_taxonomyProcessed$collectDate))))
totalCollectTOMB<-sum(invertSummaryTOMB$totalCount,na.rm=T)
invertSummaryTOMB$avgCountTOMB=invertSummaryTOMB$totalCount/boutsTOMB
invertSummaryTOMB$TOMBpercent=invertSummaryTOMB$totalCount/totalCollectTOMB

invertSummaryAll<-merge(invertSummaryMAYF,invertSummaryBLWA,,by.x="order",by.y="order",all.x=T,all.y=T)
invertSummaryAll<-merge(invertSummaryAll,invertSummaryTOMB,,by.x="order",by.y="order",all.x=T,all.y=T)