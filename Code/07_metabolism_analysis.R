########################################################################################################
#' @title metabolism analysis

#' @author Bobby Hensley email: hensley@battelleecology.org

#' @description R script which downloads, cleans, and fits the streamMetabolizer model to NEON D08 data
#'
#' Created 8/18/2021
########################################################################################################
#' Loads required packages from library 
library(neonUtilities)
library(tidyr)
library(ggplot2)
library(lubridate)
library(plyr)
library(streamMetabolizer)


#' Downloads water quality data from NEON API into R environment
waq<-loadByProduct(dpID="DP1.20288.001", site="MAYF", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(waq, .GlobalEnv)
waqMAYF<-waq_instantaneous[(waq_instantaneous$horizontalPosition=="102"),]
waq<-loadByProduct(dpID="DP1.20288.001", site="BLWA", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(waq, .GlobalEnv)
waqBLWA<-waq_instantaneous[(waq_instantaneous$horizontalPosition=="103"),]
waq<-loadByProduct(dpID="DP1.20288.001", site="TOMB", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(waq, .GlobalEnv)
waqTOMB<-waq_instantaneous[(waq_instantaneous$horizontalPosition=="103"),]


#' Creates dataframe of just DO measurements and removes any quality flagged measurements
doMAYF<-waqMAYF[,c("startDateTime","dissolvedOxygen","dissolvedOxygenFinalQF","dissolvedOxygenSaturation","dissolvedOxygenSatFinalQF")]
doMAYF<-doMAYF[(doMAYF$dissolvedOxygenFinalQF=="0"),]
doBLWA<-waqBLWA[,c("startDateTime","dissolvedOxygen","dissolvedOxygenFinalQF","dissolvedOxygenSaturation","dissolvedOxygenSatFinalQF")]
doBLWA<-doBLWA[(doBLWA$dissolvedOxygenFinalQF=="0"),]
doTOMB<-waqTOMB[,c("startDateTime","dissolvedOxygen","dissolvedOxygenFinalQF","dissolvedOxygenSaturation","dissolvedOxygenSatFinalQF")]
doTOMB<-doTOMB[(doTOMB$dissolvedOxygenFinalQF=="0"),]


#' Calculates an hourly running average and fills any data gaps (up to 6 hrs) using a spline curve 
doMAYF$startDateTime<-lubridate::round_date(doMAYF$startDateTime,unit="60 minute")
doMAYF<-plyr::ddply(doMAYF,c("startDateTime"),summarise,DO.obs=mean(dissolvedOxygen),dissolvedOxygenSaturation=mean(dissolvedOxygenSaturation))
doMAYF$DO.obs<-zoo::na.spline(doMAYF$DO.obs,maxgap=6)
doMAYF$dissolvedOxygenSaturation<-zoo::na.spline(doMAYF$dissolvedOxygenSaturation,maxgap=6)
doBLWA$startDateTime<-lubridate::round_date(doBLWA$startDateTime,unit="60 minute")
doBLWA<-plyr::ddply(doBLWA,c("startDateTime"),summarise,DO.obs=mean(dissolvedOxygen),dissolvedOxygenSaturation=mean(dissolvedOxygenSaturation))
doBLWA$DO.obs<-zoo::na.spline(doBLWA$DO.obs,maxgap=6)
doBLWA$dissolvedOxygenSaturation<-zoo::na.spline(doBLWA$dissolvedOxygenSaturation,maxgap=6)
doTOMB$startDateTime<-lubridate::round_date(doTOMB$startDateTime,unit="60 minute")
doTOMB<-plyr::ddply(doTOMB,c("startDateTime"),summarise,DO.obs=mean(dissolvedOxygen),dissolvedOxygenSaturation=mean(dissolvedOxygenSaturation))
doTOMB$DO.obs<-zoo::na.spline(doTOMB$DO.obs,maxgap=6)
doTOMB$dissolvedOxygenSaturation<-zoo::na.spline(doTOMB$dissolvedOxygenSaturation,maxgap=6)


#' Downloads elevation of surface water data from NEON API into R environment
eos<-loadByProduct(dpID="DP1.20016.001", site="MAYF", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(eos, .GlobalEnv)
eosMAYF<-EOS_30_min[(EOS_30_min$horizontalPosition=="132"),]
eos<-loadByProduct(dpID="DP1.20016.001", site="BLWA", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(eos, .GlobalEnv)
eosBLWA<-EOS_30_min
eos<-loadByProduct(dpID="DP1.20016.001", site="TOMB", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(eos, .GlobalEnv)
eosTOMB<-EOS_30_min


#' Converts elevation to estimated mean depth. For MAYF this comes from salt releases (Q/u/w=d). For BLWA and TOMB it is estimated from ADCP profiles (at similar stage 7.6 for BLWA and 9.2 for TOMB)
eosMAYF<-eosMAYF[,c("startDateTime","surfacewaterElevMean")]
eosMAYF$startDateTime<-lubridate::round_date(eosMAYF$startDateTime,unit="60 minute")
eosMAYF<-plyr::ddply(eosMAYF,c("startDateTime"),summarise,elevation=mean(surfacewaterElevMean))
eosMAYF$elevation<-zoo::na.spline(eosMAYF$elevation,maxgap=6)
for(i in 1:nrow(eosMAYF)){if(is.na(eosMAYF[i,2])){eosMAYF[i,2]=70.35}}
eosMAYF$depth<-eosMAYF$elevation-70
eosBLWA<-eosBLWA[,c("startDateTime","surfacewaterElevMean")]
eosBLWA$startDateTime<-lubridate::round_date(eosBLWA$startDateTime,unit="60 minute")
eosBLWA<-plyr::ddply(eosBLWA,c("startDateTime"),summarise,elevation=mean(surfacewaterElevMean))
eosBLWA$elevation<-zoo::na.spline(eosBLWA$elevation,maxgap=6)
for(i in 1:nrow(eosBLWA)){if(is.na(eosBLWA[i,2])){eosBLWA[i,2]=22.55}}
eosBLWA$depth<-eosBLWA$elevation-15
eosTOMB<-eosTOMB[,c("startDateTime","surfacewaterElevMean")]
eosTOMB$startDateTime<-lubridate::round_date(eosTOMB$startDateTime,unit="60 minute")
eosTOMB<-plyr::ddply(eosTOMB,c("startDateTime"),summarise,elevation=mean(surfacewaterElevMean))
eosTOMB$elevation<-zoo::na.spline(eosTOMB$elevation,maxgap=6)
for(i in 1:nrow(eosTOMB)){if(is.na(eosTOMB[i,2])){eosTOMB[i,2]=10.51}}
eosTOMB$depth<-eosTOMB$elevation-1.3


#' Downloads MAYF surface water temperature data from NEON API into R environment
tsw<-loadByProduct(dpID="DP1.20053.001", site="MAYF", startdate="2019-08", enddate="2019-08", check.size = F)
list2env(tsw, .GlobalEnv)
tswMAYF<-TSW_30min[(TSW_30min$horizontalPosition=="102"),]
tswMAYF<-tswMAYF[,c("startDateTime","surfWaterTempMean")]
tswMAYF$startDateTime<-lubridate::round_date(tswMAYF$startDateTime,unit="60 minute")
tswMAYF<-plyr::ddply(tswMAYF,c("startDateTime"),summarise,temp.water=mean(surfWaterTempMean))
tswMAYF$temp.water<-zoo::na.spline(tswMAYF$temp.water,maxgap=6)
#' For BLWA and TOMB no temp chain was deployed at this time. Need to read in csv files of L0 data created using restR.
tswBLWA<-read.csv(file="tswBLWA.csv")
tswBLWA$startDateTime<-as.POSIXct(tswBLWA$startDateTime,format="%Y-%m-%d %H:%M", tz="UTC")
tswTOMB<-read.csv(file="tswTOMB.csv")
tswTOMB$startDateTime<-as.POSIXct(tswTOMB$startDateTime,format="%Y-%m-%d %H:%M", tz="UTC")


#' Merges dataframes and formats for streamMetabolizer model
metabDataMAYF<-merge(doMAYF,eosMAYF,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
metabDataMAYF<-merge(metabDataMAYF,tswMAYF,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
metabDataBLWA<-merge(doBLWA,eosBLWA,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
metabDataBLWA<-merge(metabDataBLWA,tswBLWA,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
metabDataTOMB<-merge(doTOMB,eosTOMB,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
metabDataTOMB<-merge(metabDataTOMB,tswTOMB,by.x="startDateTime",by.y="startDateTime",all.x=T,all.y=F)
siteLong=-87.41
siteLat=32.96
metabDataMAYF$DO.sat<-metabDataMAYF$DO.obs*(100/metabDataMAYF$dissolvedOxygenSaturation)
metabDataMAYF$solar.time<-streamMetabolizer::calc_solar_time(metabDataMAYF$startDateTime,longitude=siteLong)
metabDataMAYF$light<-streamMetabolizer::calc_light(metabDataMAYF$solar.time,siteLat,siteLong,max.PAR=2000) 
metabDataMAYF<-metabDataMAYF[,c("solar.time","DO.obs","DO.sat","depth","temp.water","light")]
siteLong=-87.80
siteLat=32.54
metabDataBLWA$DO.sat<-metabDataBLWA$DO.obs*(100/metabDataBLWA$dissolvedOxygenSaturation)
metabDataBLWA$solar.time<-streamMetabolizer::calc_solar_time(metabDataBLWA$startDateTime,longitude=siteLong)
metabDataBLWA$light<-streamMetabolizer::calc_light(metabDataBLWA$solar.time,siteLat,siteLong,max.PAR=2000) 
metabDataBLWA<-metabDataBLWA[,c("solar.time","DO.obs","DO.sat","depth","temp.water","light")]
siteLong=-88.16
siteLat=31.85
metabDataTOMB$DO.sat<-metabDataTOMB$DO.obs*(100/metabDataTOMB$dissolvedOxygenSaturation)
metabDataTOMB$solar.time<-streamMetabolizer::calc_solar_time(metabDataTOMB$startDateTime,longitude=siteLong)
metabDataTOMB$light<-streamMetabolizer::calc_light(metabDataTOMB$solar.time,siteLat,siteLong,max.PAR=2000) 
metabDataTOMB<-metabDataTOMB[,c("solar.time","DO.obs","DO.sat","depth","temp.water","light")]


#############################################################################################################################################################################################################

#' Plots MAYF data for visual inspection  
metabDataMAYF %>% unitted::v() %>%
        mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
        select(solar.time, starts_with('DO')) %>%
        gather(type, DO.value, starts_with('DO')) %>%
        mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
        ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +facet_grid(units ~ ., scale='free_y') + theme_bw() + scale_color_discrete('variable')
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
metabDataMAYF %>% unitted::v() %>%
        select(solar.time, depth, temp.water, light) %>%
        gather(type, value, depth, temp.water, light) %>%
        mutate(
                type=ordered(type, levels=c('depth','temp.water','light')),
                units=ordered(labels[type], unname(labels))) %>%
        ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
        facet_grid(units ~ ., scale='free_y') + theme_bw() +
        scale_color_discrete('variable')


#' Runs mle model for MAYF
mle_name <- mm_name(type='mle')
mle_specs<-streamMetabolizer::specs(mle_name)
mle_specs
metabModel<-streamMetabolizer::metab(mle_specs,data=metabDataMAYF)
metabModel
plot_metab_preds(metabModel,y_lim=list(GPP = c(-10, 20), ER = c(-25, 35))) 
plot_DO_preds(metabModel)
metabResultsMAYF<-get_params(metabModel)
metabResultsMAYF$PR<-metabResultsMAYF$GPP.daily/metabResultsMAYF$ER.daily*-1


#' Plots BLWA data for visual inspection  
metabDataBLWA %>% unitted::v() %>%
        mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
        select(solar.time, starts_with('DO')) %>%
        gather(type, DO.value, starts_with('DO')) %>%
        mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
        ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +facet_grid(units ~ ., scale='free_y') + theme_bw() + scale_color_discrete('variable')
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
metabDataBLWA %>% unitted::v() %>%
        select(solar.time, depth, temp.water, light) %>%
        gather(type, value, depth, temp.water, light) %>%
        mutate(
                type=ordered(type, levels=c('depth','temp.water','light')),
                units=ordered(labels[type], unname(labels))) %>%
        ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
        facet_grid(units ~ ., scale='free_y') + theme_bw() +
        scale_color_discrete('variable')


#' Runs mle model for BLWA
mle_name <- mm_name(type='mle')
mle_specs<-streamMetabolizer::specs(mle_name)
mle_specs
metabModel<-streamMetabolizer::metab(mle_specs,data=metabDataBLWA)
metabModel
plot_metab_preds(metabModel,y_lim=list(GPP = c(-10, 20), ER = c(-25, 35))) 
plot_DO_preds(metabModel)
metabResultsBLWA<-get_params(metabModel)
metabResultsBLWA$PR<-metabResultsBLWA$GPP.daily/metabResultsBLWA$ER.daily*-1


#' Plots TOMB data for visual inspection  
metabDataTOMB %>% unitted::v() %>%
        mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
        select(solar.time, starts_with('DO')) %>%
        gather(type, DO.value, starts_with('DO')) %>%
        mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
        ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +facet_grid(units ~ ., scale='free_y') + theme_bw() + scale_color_discrete('variable')
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
metabDataTOMB %>% unitted::v() %>%
        select(solar.time, depth, temp.water, light) %>%
        gather(type, value, depth, temp.water, light) %>%
        mutate(
                type=ordered(type, levels=c('depth','temp.water','light')),
                units=ordered(labels[type], unname(labels))) %>%
        ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
        facet_grid(units ~ ., scale='free_y') + theme_bw() +
        scale_color_discrete('variable')


#' Runs mle model for TOMB
mle_name <- mm_name(type='mle')
mle_specs<-streamMetabolizer::specs(mle_name)
mle_specs
metabModel<-streamMetabolizer::metab(mle_specs,data=metabDataTOMB)
metabModel
plot_metab_preds(metabModel,y_lim=list(GPP = c(-10, 20), ER = c(-25, 35))) 
plot_DO_preds(metabModel)
metabResultsTOMB<-get_params(metabModel)
metabResultsTOMB$PR<-metabResultsTOMB$GPP.daily/metabResultsTOMB$ER.daily*-1

