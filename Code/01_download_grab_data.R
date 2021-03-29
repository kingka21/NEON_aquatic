#### Script for downloading and organizing NEON datasets
### code created by Katelyn King 16-OCT-2019, modified by KK on Sep28 
### code adapted from Bobby Hensley 
###Modified by Jennifer Edmonds 21Sep2020


#### load libraries ####
# install neonUtilities if you have not already 
# C:\Users\INBRE\Documents\NEON_aquatic
install.packages("neonUtilities")
#load neonUtilities
library(stringr)
library(neonUtilities)
library(ggplot2)
library(dplyr)
library(tidyr)

#### load datasets ####
## avg: "all" to download all data (default), or number of minutes in the averaging interval. 
##      only applicable to IS data.

## site: defaults to "all", meaning all sites with available data; 
##      ARIK =  Arikaree River CO       BARC = Barco Lake FL          BIGC = Upper Big Creek CA   
##      BLDE = Black Deer Creek WY      BLUE = Blue River OK          BLWA = Black Warrior River AL
##      CARI = Caribou Creek AK         COMO = Como Creek CO          CRAM = Crampton Lake WI     
##      CUPE = Rio Cupeyes PR           FLNT = Flint River GA         GUIL = Rio Guilarte PR
##      HOPB = Lower Hop Brook MA       KING = Kings Creek KS  L      ECO = LeConte Creek TN     
##      LEWI = Lewis Run VA             LIRO = Little Rock Lake WI    MART = Marta Creek WA
##      MAYF = Mayfield Creek AL        MCDI = McDiffett Creek KS     MCRA = McRae Creek OR       
##      OKSR = Oksrukuyik Creek AK      POSE = Posey Creek VA         PRIN = Pringle Creek TX       
##      PRLA = Prairie Lake ND          PRPO = Prairie Pothole ND     REDB = Red Butte Creek UT
##      SUGG = Suggs Lake FL            SYCA = Sycamore Creek AZ      TECR = Teakettle Creek CA        
##      TOMB = Lower Tombigbee River AL TOOK = Toolik Lake AK         WALK = Walker Branch TN
##      WLOU = West St Louis Creek CO       

## Need the data product ID (dpID)
##      20016.001 = elevation of surface water
##      20288.001 = water quality
##      20033.001 = nitrate in surface water
##      20092.001 = groundwater 
##      20163.001 = chl
##      20206.001 = isotopes in surface water 
##      20276.001 = isotopes in groundwater
##      20048.001 = manual discharge measures (about 26/year/site) (correct)


#####All 27 sites surface water grab ####
#Grab samples of surface water chemistry including general chemistry (DOC), anions, cations, and nutrients. streams 26 times per year
nutrients_SWgrab_allsites <- loadByProduct(dpID="DP1.20093.001", 
                           site=c(
                             'HOPB','POSE','KING', 'WALK','LECO','MAYF','PRIN',
                             'BLDE','COMO','MART', 'BIGC','CARI', 'WLOU', "FLNT", 
                             "MCDI", 'LEWI', "BLUE", "TECR", "REDB", "SYCA", 
                             "MCRA", "OKSR", "ARIK", "GUIL", "CUPE", "TOMB", "BLWA"
                           ),
                           startdate="2017-01", 
                           enddate="2020-09", 
                           package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                           check.size = F)  ### check the size of the file before you download it 

for(i in 1:length(nutrients_SWgrab_allsites)) {assign(names(nutrients_SWgrab_allsites)[i], nutrients_SWgrab_allsites[[i]])}   #calls the table waq_instances

SWgrab_chem_dat_allsites<-as.data.frame(swc_externalLabDataByAnalyte) #table that has water chem #contains 23,340 observations
SWgrab_info_dat_allsites<-as.data.frame(swc_fieldSuperParent) #table that has lat/lon and elevation, DO, waterTemp, maxDepth
SWgrab_lab_dat_allsites<-as.data.frame(swc_domainLabData)#table that has pH values that we can trust

#first step is to remove lines where sampleID ends with .2 or .3#
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites, !grepl("FIL.3$", sampleID))
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, !grepl("FIL.2$", sampleID))
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, !grepl("RAW.3$", sampleID))
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, !grepl("RAW.2$", sampleID))
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, !grepl("PCN.3$", sampleID))
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, !grepl("PCN.2$", sampleID))
#contains 60,945 observations

SWgrab_chem_dat_allsites_QC$sampleID<-gsub(".FIL" , "" , SWgrab_chem_dat_allsites_QC$sampleID)
SWgrab_chem_dat_allsites_QC$sampleID<-gsub(".RAW" , "" , SWgrab_chem_dat_allsites_QC$sampleID)
SWgrab_chem_dat_allsites_QC$sampleID<-gsub(".PCN" , "" , SWgrab_chem_dat_allsites_QC$sampleID)

#Remove flagged data, shipmentWarmQF == 0 indicates no flags # 
SWgrab_chem_dat_allsites_QC<-filter(SWgrab_chem_dat_allsites_QC, shipmentWarmQF == 0)

#Setting negative values for nutrients to 0# 53269 observations 
SWgrab_chem_dat_allsites_QC[SWgrab_chem_dat_allsites_QC <0] <- 0 #change negative values to 0 

#generate new table using pivot data for water chemistry parameters# #1652 rows 

SWgrab_chem_dat_allsites_PIVOT<-pivot_wider(SWgrab_chem_dat_allsites_QC, 
                                            id_cols= c(siteID, sampleID, collectDate), 
                                            names_from=analyte, 
                                            values_from=analyteConcentration,
                                            values_fn = list(analyteConcentration = mean)) ## TDN had one instance where it was repeated, this takes the average of just those

#select only dissolved variables and otheres discussed in meeting on Sept 21,2020
#left with 18 variables 
swc_all<-select(SWgrab_chem_dat_allsites_PIVOT, siteID, sampleID, collectDate, Br, Ca, Cl, DIC, DOC, F, Fe, K, Mg, Mn, Na, 'NH4 - N',
                'NO3+NO2 - N', pH, Si, SO4, TDP, TDS)

swc_all<-na.omit(swc_all) # get rid of rows with NA, 1504

#replaces reported values below the detection limit with a value equal to half the detection limit
#Bobby Hensley (9/28/2020) 
swc_all$collectDate<-as.POSIXct(swc_all$collectDate,format="%m/%d/%Y %H:%M", tz="UTC")

#' Replace Br non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,4]<=0){swc_all[i,4]=0.005}}
#' Replace Ca non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,5]<=0){swc_all[i,5]=0.0005}}
#' Replace Cl non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,6]<=0){swc_all[i,6]=0.005}}
#' Replace DIC non-dectects with half detection limit (0.100 prior to 2017-02-09, 0.025 after)
changeDate<-as.POSIXct("02/09/2017 00:00",format="%m/%d/%Y %H:%M", tz="UTC")
for(i in 1:nrow(swc_all)){if(swc_all[i,7]<=0){if(swc_all[i,3]<=changeDate){swc_all[i,7]=0.050}else{swc_all[i,7]=0.0125}}}
#' Replace DIC non-dectects with half detection limit (0.100 prior to 2019-05-28, 0.970 after)
changeDate<-as.POSIXct("05/28/2019 00:00",format="%m/%d/%Y %H:%M", tz="UTC")
for(i in 1:nrow(swc_all)){if(swc_all[i,8]<=0){if(swc_all[i,3]<=changeDate){swc_all[i,8]=0.050}else{swc_all[i,8]=0.0485}}}
#' Replace F non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,9]<=0){swc_all[i,9]=0.005}}
#' Replace Fe non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,10]<=0){swc_all[i,10]=0.0005}}
#' Replace K non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,11]<=0){swc_all[i,11]=0.0005}}
#' Replace Mg non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,12]<=0){swc_all[i,12]=0.005}}
#' Replace Mn non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,13]<=0){swc_all[i,13]=0.0005}}
#' Replace Na non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,14]<=0){swc_all[i,14]=0.0005}}
#' Replace NH4N non-dectects with half detection limit (0.020 prior to 2019-07-03, 0.004 after)
changeDate<-as.POSIXct("05/28/2019 00:00",format="%m/%d/%Y %H:%M", tz="UTC")
for(i in 1:nrow(swc_all)){if(swc_all[i,15]<=0){if(swc_all[i,3]<=changeDate){swc_all[i,15]=0.010}else{swc_all[i,15]=0.002}}}
#' Replace NO3NO2 non-dectects with half detection limit (0.027)
for(i in 1:nrow(swc_all)){if(swc_all[i,16]<=0){swc_all[i,16]=0.0135}}
#' Replace Si non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,18]<=0){swc_all[i,18]=0.005}}
#' Replace SO4 non-dectects with half detection limit (0.010)
for(i in 1:nrow(swc_all)){if(swc_all[i,19]<=0){swc_all[i,19]=0.005}}
#' Replace TDP non-dectects with half detection limit (0.001)
for(i in 1:nrow(swc_all)){if(swc_all[i,20]<=0){swc_all[i,20]=0.0005}}
#' Replace TDS non-dectects with half detection limit (0.100)
for(i in 1:nrow(swc_all)){if(swc_all[i,21]<=0){swc_all[i,21]=0.050}}


##Join tables containing chemistry data and pH domain data##
#first remove all sample lines in domain data containing ANC
SWgrab_lab_dat_allsites_noANC<-filter(SWgrab_lab_dat_allsites, !grepl(".ANC", domainSampleID))

#next remove all sample lines in domain data containing ALK.2 and REP.2 and REP.3 and .1
SWgrab_lab_dat_allsites_QC<-filter(SWgrab_lab_dat_allsites_noANC, !grepl(".ALK.2", domainSampleID))
SWgrab_lab_dat_allsites_QC<-filter(SWgrab_lab_dat_allsites_QC, !grepl(".REP2", domainSampleID))
SWgrab_lab_dat_allsites_QC<-filter(SWgrab_lab_dat_allsites_QC, !grepl(".REP3", domainSampleID))
SWgrab_lab_dat_allsites_QC<-filter(SWgrab_lab_dat_allsites_QC, !grepl("ALK.1", domainSampleID))

#Remove all columns except those needed to meld tables, leaving only pH from domain dataset
SWgrab_lab_pH_allsites<-select(SWgrab_lab_dat_allsites_QC, parentSampleID, initialSamplepH) 

#join table with external lab data with table from domain lab that has pH in it
SWgrab_chem_dat_allsites_pH <- left_join(swc_all, SWgrab_lab_pH_allsites,  by =c("sampleID" = "parentSampleID"))

write.csv(SWgrab_chem_dat_allsites_pH, 'Data/surface_water_grab_QCQC_ph.csv', row.names = FALSE)



