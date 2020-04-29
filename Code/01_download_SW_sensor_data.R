#### Script for downloading and organizing NEON datasets ##### 
### code written by Katelyn King 16-OCT-2019 
### code adapted from Bobby Hensley ### 

#### load libraries ####
# install neonUtilities if you have not already 
install.packages("neonUtilities")
# load neonUtilities
library(neonUtilities)
library(ggplot2)
library(dplyr)

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
##      20033.001 = nitrate in surface water sensor
##      20092.001 = grab samples -groundwater 
##      20163.001 = chl in grab
##      20206.001 = isotopes in surface water grab 
##      20276.001 = isotopes in groundwater
##      20048.001 = manual discharge measures (about 26/year/site)
##      20042.011 = PAR
####Water quality SENSOR DATA ####
#In situ sensor-based specific conductivity, concentration of chlorophyll a, dissolved oxygen content, fDOM concentration (fluorescent dissolved material), pH, and turbidity 
#available as one-, five-, and thirty-minute averages in surface water 
## pull out just one of the sensors horizontal position
##  Upstream = 101              Downstream = 102
##  US overhanging = 111        DS overhanging = 112
##  Lake/river buoy = 103
    ##  NOTE: SUNA and fDOM only at downstream stream station or lake/river buoy
#specificConductance, dissolvedOxygen, pH, chlorophyll, turbidity,fDOM, nitrate? water level? isotopes?  temp? 
#10 sites: ARIK, COMO, KING, MAYF, TOMB, BLUE, CUPE, MART, HOPB, WALK
#FDOM, chla, conductivity 
####Chlorophyll a ####
##ARIK site 
ARIK_sensor <- loadByProduct(dpID="DP1.20288.001", 
                           site=c("ARIK"),
                           startdate="2017-01",  #year and month
                           enddate="2019-12" ,
                           package="expanded",
                           check.size = T , 
                           avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

for(i in 1:length(ARIK_sensor)) {assign(names(ARIK_sensor)[i], ARIK_sensor[[i]])}   #calls the table 

ARIK<-waq_instantaneous #table with data is waq_instantaneous
ARIK_QA<-select(ARIK, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                  chlorophyllNullQF, chlorophyllGapQF,
                  chlorophyllSpikeQF, chlorophyllValidCalQF,
                  chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                  chlorophyllAlphaQF, chlorophyllBetaQF,
                  chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
ARIK_QA[is.na(ARIK_QA)] <- 0 #change NAs to 0 
ARIK_QA <-filter(ARIK_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                   chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                  chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                  chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2226260 obs to 1064315 after flags removed 
ARIK_QA<-filter(ARIK_QA, horizontalPosition == "102" ) #102= downstream sensor 
#ARIK_QA$startDateTime<-as.POSIXct(ARIK_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

## COMO and KING sites 
comoking_sensor <- loadByProduct(dpID="DP1.20288.001", 
                                  site=c('COMO', 'KING'),
                                  startdate="2017-01",  #year and month
                                  enddate="2019-12" ,
                                  package="expanded",
                                  check.size = T , 
                                  avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.
for(i in 1:length(comoking_sensor)) {assign(names(comoking_sensor)[i], comoking_sensor[[i]])}   #calls the table 
#table with data is waq_instantaneous
COMO_KING<-waq_instantaneous
COMO_KING_QA<-select(COMO_KING, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                chlorophyllNullQF, chlorophyllGapQF,
                chlorophyllSpikeQF, chlorophyllValidCalQF,
                chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                chlorophyllAlphaQF, chlorophyllBetaQF,
                chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
COMO_KING_QA[is.na(COMO_KING_QA)] <- 0 #change NAs to 0 
COMO_KING_QA <-filter(COMO_KING_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                   chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                   chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                   chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 4691649 obs to 2405831 after flags removed 
COMO_KING_QA<-filter(COMO_KING_QA, horizontalPosition == "102" ) #102= downstream sensor 
#COMO_KING_QA$startDateTime<-as.POSIXct(COMO_KING_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

#MAYF and TOMB sites 
mayftomb_sensor <- loadByProduct(dpID="DP1.20288.001", 
                                  site=c( 'MAYF', "TOMB"
                                          ),
                                  startdate="2017-01",  #year and month
                                  enddate="2019-12" ,
                                  package="expanded",
                                  check.size = T , 
                                  avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

for(i in 1:length(mayftomb_sensor)) {assign(names(mayftomb_sensor)[i], mayftomb_sensor[[i]])}   #calls the table 
#table with data is waq_instantaneous
MAYF_TOMB<-waq_instantaneous
MAYF_TOMB_QA<-select(MAYF_TOMB, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                     chlorophyllNullQF, chlorophyllGapQF,
                     chlorophyllSpikeQF, chlorophyllValidCalQF,
                     chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                     chlorophyllAlphaQF, chlorophyllBetaQF,
                     chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
MAYF_TOMB_QA[is.na(MAYF_TOMB_QA)] <- 0 #change NAs to 0 
MAYF_TOMB_QA <-filter(MAYF_TOMB_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                        chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                        chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                        chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2264882 obs to 1491737 after flags removed 
MAYF_TOMB_QA<-filter(MAYF_TOMB_QA, horizontalPosition == "102" ) #102= downstream sensor 
#MAYF_TOMB_QA$startDateTime<-as.POSIXct(MAYF_TOMB_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

#BLUE and CUPE 
bluecupe_sensor <- loadByProduct(dpID="DP1.20288.001", 
                                  site=c( "BLUE", "CUPE"),
                                  startdate="2017-01",  #year and month
                                  enddate="2019-12" ,
                                  package="expanded",
                                  check.size = T , 
                                  avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.
for(i in 1:length(bluecupe_sensor)) {assign(names(bluecupe_sensor)[i], bluecupe_sensor[[i]])}   #calls the table 
#table with data is waq_instantaneous
BLUE_CUPE<-waq_instantaneous
BLUE_CUPE_QA<-select(BLUE_CUPE, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                     chlorophyllNullQF, chlorophyllGapQF,
                     chlorophyllSpikeQF, chlorophyllValidCalQF,
                     chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                     chlorophyllAlphaQF, chlorophyllBetaQF,
                     chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
BLUE_CUPE_QA[is.na(BLUE_CUPE_QA)] <- 0 #change NAs to 0 
BLUE_CUPE_QA <-filter(BLUE_CUPE_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                        chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                        chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                        chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2730975 obs to 1197141 after flags removed 
BLUE_CUPE_QA<-filter(BLUE_CUPE_QA, horizontalPosition == "102" ) #102= downstream sensor 
#BLUE_CUPE_QA$startDateTime<-as.POSIXct(BLUE_CUPE_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 


## Mart 
mart_sensor <- loadByProduct(dpID="DP1.20288.001", 
                                  site=c('MART'),
                                  startdate="2017-01",  #year and month
                                  enddate="2019-12" ,
                                  package="expanded",
                                  check.size = T , 
                                  avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

for(i in 1:length(mart_sensor)) {assign(names(mart_sensor)[i], mart_sensor[[i]])}   #calls the table 
#table with data is waq_instantaneous
MART<-waq_instantaneous
MART_QA<-select(MART, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                     chlorophyllNullQF, chlorophyllGapQF,
                     chlorophyllSpikeQF, chlorophyllValidCalQF,
                     chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                     chlorophyllAlphaQF, chlorophyllBetaQF,
                     chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
MART_QA[is.na(MART_QA)] <- 0 #change NAs to 0 
MART_QA <-filter(MART_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                        chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                        chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                        chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2044465 obs to 953432 after flags removed 
MART_QA<-filter(MART_QA, horizontalPosition == "102" ) #102= downstream sensor 
#MART_QA$startDateTime<-as.POSIXct(MART_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

## HOPB site first half (too big to do all at once)
hopb_sensor1 <- loadByProduct(dpID="DP1.20288.001", 
                             site=c( 'HOPB'),
                             startdate="2017-01",  #year and month
                             enddate="2018-12" ,
                             package="expanded",
                             check.size = T , 
                             avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.
for(i in 1:length(hopb_sensor1)) {assign(names(hopb_sensor1)[i], hopb_sensor1[[i]])}   #calls the table 
#table with data is waq_instantaneous
HOPB1<-waq_instantaneous
HOPB1_QA<-select(HOPB1, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                     chlorophyllNullQF, chlorophyllGapQF,
                     chlorophyllSpikeQF, chlorophyllValidCalQF,
                     chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                     chlorophyllAlphaQF, chlorophyllBetaQF,
                     chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
HOPB1_QA[is.na(HOPB1_QA)] <- 0 #change NAs to 0 
HOPB1_QA <-filter(HOPB1_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                        chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                        chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                        chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2706902 obs to 410405 after flags removed 
HOPB1_QA<-filter(HOPB1_QA, horizontalPosition == "102" ) #102= downstream sensor 
#HOPB1_QA$startDateTime<-as.POSIXct(HOPB1_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

#hopb sites half 2 
hopb_sensor2 <- loadByProduct(dpID="DP1.20288.001", 
                              site=c( 'HOPB'),
                              startdate="2019-01",  #year and month
                              enddate="2019-12" ,
                              package="expanded",
                              check.size = T , 
                              avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

for(i in 1:length(hopb_sensor2)) {assign(names(hopb_sensor2)[i], hopb_sensor2[[i]])}   #calls the table 
#table with data is waq_instantaneous
HOPB2<-waq_instantaneous
HOPB2_QA<-select(HOPB2, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                     chlorophyllNullQF, chlorophyllGapQF,
                     chlorophyllSpikeQF, chlorophyllValidCalQF,
                     chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                     chlorophyllAlphaQF, chlorophyllBetaQF,
                     chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
HOPB2_QA[is.na(HOPB2_QA)] <- 0 #change NAs to 0 
HOPB2_QA <-filter(HOPB2_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                        chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                        chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                        chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2045867 obs to 558361 after flags removed 
HOPB2_QA<-filter(HOPB2_QA, horizontalPosition == "102" ) #102= downstream sensor 
#HOPB2_QA$startDateTime<-as.POSIXct(HOPB2_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

#WALK site
walk_sensor <- loadByProduct(dpID="DP1.20288.001", 
                                 site=c( "WALK"),
                                 startdate="2017-01",  #year and month
                                 enddate="2019-12" ,
                                 package="expanded",
                                 check.size = T , 
                                 avg= 'all') # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

for(i in 1:length(walk_sensor)) {assign(names(walk_sensor)[i], walk_sensor[[i]])}   #calls the table 
#table with data is waq_instantaneous
WALK<-waq_instantaneous
WALK_QA<-select(WALK, siteID, horizontalPosition, startDateTime, chlorophyll, chlorophyllRangeQF, chlorophyllStepQF,
                 chlorophyllNullQF, chlorophyllGapQF,
                 chlorophyllSpikeQF, chlorophyllValidCalQF,
                 chlorophyllSuspectCalQF, chlorophyllPersistenceQF,
                 chlorophyllAlphaQF, chlorophyllBetaQF,
                 chlorophyllFinalQF, chlorophyllFinalQFSciRvw) 
WALK_QA[is.na(WALK_QA)] <- 0 #change NAs to 0 
WALK_QA <-filter(WALK_QA, chlorophyllRangeQF == 0 & chlorophyllStepQF == 0 &
                    chlorophyllNullQF == 0 & chlorophyllGapQF == 0 & chlorophyllSpikeQF == 0 &
                    chlorophyllValidCalQF == 0 & chlorophyllSuspectCalQF == 0 & chlorophyllPersistenceQF == 0 & chlorophyllAlphaQF == 0 &
                    chlorophyllBetaQF == 0 & chlorophyllFinalQF == 0 & chlorophyllFinalQFSciRvw == 0) 
# 2131068 obs to 1202093 after flags removed 
WALK_QA<-filter(WALK_QA, horizontalPosition == "102" ) #102= downstream sensor 
#WALK_QA$startDateTime<-as.POSIXct(WALK_QA$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 

#merge all the tables into one 
chl_sensor_10sites<- gtools::smartbind(ARIK_QA, BLUE_CUPE_QA, 
                               COMO_KING_QA, HOPB1_QA, HOPB2_QA, 
                               MART_QA, MAYF_TOMB_QA, WALK_QA)

#put year, day, month and time in a different column?! 
chl_sensor_10sites$startDateTime<-as.POSIXct(chl_sensor_10sites$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 
chl_sensor_10sites$DATE<-as.Date(chl_sensor_10sites$startDateTime,format="%Y-%m-%d")
chl_sensor_10sites$YEAR<-format(chl_sensor_10sites$DATE,format="%Y")
chl_sensor_10sites$MONTH<-format(chl_sensor_10sites$DATE,format="%m")
chl_sensor_10sites$DAY<-format(chl_sensor_10sites$DATE,format="%d")

chl_sensor<-select(chl_sensor_10sites, siteID, chlorophyll, startDateTime, DATE, YEAR, MONTH, DAY)

write.csv(chl_sensor, 'Data/sw_chl_sensor.csv', row.names = FALSE)

#read in data 
#sensor1<-read.csv('Data/surface_water_sensor.csv')
#waq_instantaneous$monthYr = format(as.Date(waq_instantaneous$DATE), "%m.%d")

#single site investigation check out diernal patterns
BLUE<-filter(chl_sensor, siteID == "BLUE" )
TOMB<-filter(chl_sensor, siteID == "TOMB" ) #no info for this site
ARIK<-filter(chl_sensor, siteID == "ARIK" )
COMO<-filter(chl_sensor, siteID == "COMO" )
KING<-filter(chl_sensor, siteID == "KING" )
MAYF<-filter(chl_sensor, siteID == "MAYF" )
CUPE<-filter(chl_sensor, siteID == "CUPE" )
MART<-filter(chl_sensor, siteID == "MART" )
HOPB<-filter(chl_sensor, siteID == "HOPB" )
WALK<-filter(chl_sensor, siteID == "WALK" )
  
#chose some days to look at winter, spring, summer, fall
nona<-WALK[!with(WALK,is.na(chlorophyll)),]
winter<-filter(nona, DATE == "2018-12-22")
spring<-filter(nona, DATE == "2019-04-29")
summer<-filter(nona, DATE == "2019-07-26")
fall<-filter(nona, DATE == "2019-10-10")

ggplot(data = summer, aes(x=startDateTime, y=chlorophyll))  +
  geom_point(aes(colour=siteID)) + 
  theme(legend.position = "none")


####Conductivity ####
ARIK_QA<-select(ARIK, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                specificConductanceNullQF, specificConductanceGapQF,
                specificConductanceSpikeQF, specificConductanceValidCalQF,
                specificCondSuspectCalQF, specificConductancePersistQF,
                specificConductanceAlphaQF, specificConductanceBetaQF,
                specificCondFinalQF, specificCondFinalQFSciRvw) 
ARIK_QA[is.na(ARIK_QA)] <- 0 #change NAs to 0 
ARIK_QA <-filter(ARIK_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                   specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                   specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                   specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
ARIK_QA<-filter(ARIK_QA, horizontalPosition == "102" ) #102= downstream sensor 

## COMO and KING sites 
COMO_KING_QA<-select(COMO_KING, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                     specificConductanceNullQF, specificConductanceGapQF,
                     specificConductanceSpikeQF, specificConductanceValidCalQF,
                     specificCondSuspectCalQF, specificConductancePersistQF,
                     specificConductanceAlphaQF, specificConductanceBetaQF,
                     specificCondFinalQF, specificCondFinalQFSciRvw) 
COMO_KING_QA[is.na(COMO_KING_QA)] <- 0 #change NAs to 0 
COMO_KING_QA <-filter(COMO_KING_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                        specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                        specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                        specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
COMO_KING_QA<-filter(COMO_KING_QA, horizontalPosition == "102" ) #102= downstream sensor 

#MAYF and TOMB sites 
MAYF_TOMB_QA<-select(MAYF_TOMB, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                     specificConductanceNullQF, specificConductanceGapQF,
                     specificConductanceSpikeQF, specificConductanceValidCalQF,
                     specificCondSuspectCalQF, specificConductancePersistQF,
                     specificConductanceAlphaQF, specificConductanceBetaQF,
                     specificCondFinalQF, specificCondFinalQFSciRvw) 
MAYF_TOMB_QA[is.na(MAYF_TOMB_QA)] <- 0 #change NAs to 0 
MAYF_TOMB_QA <-filter(MAYF_TOMB_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                        specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                        specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                        specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
MAYF_TOMB_QA<-filter(MAYF_TOMB_QA, horizontalPosition == "102" ) #102= downstream sensor 

#BLUE and CUPE 
BLUE_CUPE_QA<-select(BLUE_CUPE, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                     specificConductanceNullQF, specificConductanceGapQF,
                     specificConductanceSpikeQF, specificConductanceValidCalQF,
                     specificCondSuspectCalQF, specificConductancePersistQF,
                     specificConductanceAlphaQF, specificConductanceBetaQF,
                     specificCondFinalQF, specificCondFinalQFSciRvw) 
BLUE_CUPE_QA[is.na(BLUE_CUPE_QA)] <- 0 #change NAs to 0 
BLUE_CUPE_QA <-filter(BLUE_CUPE_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                        specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                        specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                        specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
BLUE_CUPE_QA<-filter(BLUE_CUPE_QA, horizontalPosition == "102" ) #102= downstream sensor 

## Mart 
MART_QA<-select(MART, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                specificConductanceNullQF, specificConductanceGapQF,
                specificConductanceSpikeQF, specificConductanceValidCalQF,
                specificCondSuspectCalQF, specificConductancePersistQF,
                specificConductanceAlphaQF, specificConductanceBetaQF,
                specificCondFinalQF, specificCondFinalQFSciRvw) 
MART_QA[is.na(MART_QA)] <- 0 #change NAs to 0 
MART_QA <-filter(MART_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                   specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                   specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                   specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
MART_QA<-filter(MART_QA, horizontalPosition == "102" ) #102= downstream sensor 

## HOPB site first half (too big to do all at once)
HOPB1_QA<-select(HOPB1, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                 specificConductanceNullQF, specificConductanceGapQF,
                 specificConductanceSpikeQF, specificConductanceValidCalQF,
                 specificCondSuspectCalQF, specificConductancePersistQF,
                 specificConductanceAlphaQF, specificConductanceBetaQF,
                 specificCondFinalQF, specificCondFinalQFSciRvw) 
HOPB1_QA[is.na(HOPB1_QA)] <- 0 #change NAs to 0 
HOPB1_QA <-filter(HOPB1_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                    specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                    specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                    specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
HOPB1_QA<-filter(HOPB1_QA, horizontalPosition == "102" ) #102= downstream sensor 

#hopb sites half 2 

HOPB2_QA<-select(HOPB2, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                 specificConductanceNullQF, specificConductanceGapQF,
                 specificConductanceSpikeQF, specificConductanceValidCalQF,
                 specificCondSuspectCalQF, specificConductancePersistQF,
                 specificConductanceAlphaQF, specificConductanceBetaQF,
                 specificCondFinalQF, specificCondFinalQFSciRvw) 
HOPB2_QA[is.na(HOPB2_QA)] <- 0 #change NAs to 0 
HOPB2_QA <-filter(HOPB2_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                    specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                    specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                    specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
HOPB2_QA<-filter(HOPB2_QA, horizontalPosition == "102" ) #102= downstream sensor 

#WALK site
WALK_QA<-select(WALK, siteID, horizontalPosition, startDateTime, specificConductance, specificConductanceRangeQF, specificConductanceStepQF,
                specificConductanceNullQF, specificConductanceGapQF,
                specificConductanceSpikeQF, specificConductanceValidCalQF,
                specificCondSuspectCalQF, specificConductancePersistQF,
                specificConductanceAlphaQF, specificConductanceBetaQF,
                specificCondFinalQF, specificCondFinalQFSciRvw) 
WALK_QA[is.na(WALK_QA)] <- 0 #change NAs to 0 
WALK_QA <-filter(WALK_QA, specificConductanceRangeQF == 0 & specificConductanceStepQF == 0 &
                   specificConductanceNullQF == 0 & specificConductanceGapQF == 0 & specificConductanceSpikeQF == 0 &
                   specificConductanceValidCalQF == 0 & specificCondSuspectCalQF == 0 & specificConductancePersistQF == 0 & specificConductanceAlphaQF == 0 &
                   specificConductanceBetaQF == 0 & specificCondFinalQF == 0 & specificCondFinalQFSciRvw == 0) 
WALK_QA<-filter(WALK_QA, horizontalPosition == "102" ) #102= downstream sensor 

#merge all the tables into one 
cond_sensor_10sites<- gtools::smartbind(ARIK_QA, BLUE_CUPE_QA, 
                                       COMO_KING_QA, HOPB1_QA, HOPB2_QA, 
                                       MART_QA, MAYF_TOMB_QA, WALK_QA)

#put year, day, month and time in a different column?! 
cond_sensor_10sites$startDateTime<-as.POSIXct(cond_sensor_10sites$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 
cond_sensor_10sites$DATE<-as.Date(cond_sensor_10sites$startDateTime,format="%Y-%m-%d")
cond_sensor_10sites$YEAR<-format(cond_sensor_10sites$DATE,format="%Y")
cond_sensor_10sites$MONTH<-format(cond_sensor_10sites$DATE,format="%m")
cond_sensor_10sites$DAY<-format(cond_sensor_10sites$DATE,format="%d")

cond_sensor<-select(cond_sensor_10sites, siteID, specificConductance, startDateTime, DATE, YEAR, MONTH, DAY)

write.csv(cond_sensor, 'Data/sw_cond_sensor.csv', row.names = FALSE)

#### FDOM ####
ARIK_QA<-select(ARIK, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                fDOMNullQF, fDOMGapQF,
                fDOMSpikeQF, fDOMValidCalQF,
                fDOMSuspectCalQF, fDOMPersistenceQF,
                fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                fDOMFinalQF, fDOMFinalQFSciRvw) 
ARIK_QA[is.na(ARIK_QA)] <- 0 #change NAs to 0 
ARIK_QA <-filter(ARIK_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                   fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                   fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                   fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2226260 to 1169704
ARIK_QA<-filter(ARIK_QA, horizontalPosition == "102" ) #102= downstream sensor 

## COMO and KING sites 
COMO_KING_QA<-select(COMO_KING, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                     fDOMNullQF, fDOMGapQF,
                     fDOMSpikeQF, fDOMValidCalQF,
                     fDOMSuspectCalQF, fDOMPersistenceQF, 
                     fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                     fDOMFinalQF, fDOMFinalQFSciRvw) 
COMO_KING_QA[is.na(COMO_KING_QA)] <- 0 #change NAs to 0 
COMO_KING_QA <-filter(COMO_KING_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                        fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                        fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                        fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#4691649 to 2568603
COMO_KING_QA<-filter(COMO_KING_QA, horizontalPosition == "102" ) #102= downstream sensor 

#MAYF and TOMB sites 
MAYF_TOMB_QA<-select(MAYF_TOMB, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                     fDOMNullQF, fDOMGapQF,
                     fDOMSpikeQF, fDOMValidCalQF,
                     fDOMSuspectCalQF, fDOMPersistenceQF,
                     fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                     fDOMFinalQF, fDOMFinalQFSciRvw) 
MAYF_TOMB_QA[is.na(MAYF_TOMB_QA)] <- 0 #change NAs to 0 
MAYF_TOMB_QA <-filter(MAYF_TOMB_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                        fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                        fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                        fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2264882 to 1051365
# there is only 101 (Upstream) for MAYF_TOMB for horizontalPosition even though note said only at downstream

#BLUE and CUPE 
BLUE_CUPE_QA<-select(BLUE_CUPE, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                     fDOMNullQF, fDOMGapQF,
                     fDOMSpikeQF, fDOMValidCalQF,
                     fDOMSuspectCalQF, fDOMPersistenceQF,
                     fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                     fDOMFinalQF, fDOMFinalQFSciRvw) 
BLUE_CUPE_QA[is.na(BLUE_CUPE_QA)] <- 0 #change NAs to 0 
BLUE_CUPE_QA <-filter(BLUE_CUPE_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                        fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                        fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                        fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2730975 to 1126373
BLUE_CUPE_QA<-filter(BLUE_CUPE_QA, horizontalPosition == "101" ) 

## Mart 
MART_QA<-select(MART, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                fDOMNullQF, fDOMGapQF,
                fDOMSpikeQF, fDOMValidCalQF,
                fDOMSuspectCalQF, fDOMPersistenceQF,
                fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                fDOMFinalQF, fDOMFinalQFSciRvw) 
MART_QA[is.na(MART_QA)] <- 0 #change NAs to 0 
MART_QA <-filter(MART_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                   fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                   fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                   fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2044464 to 1053916
MART_QA<-filter(MART_QA, horizontalPosition == "102" ) #102= downstream sensor 

## HOPB site first half (too big to do all at once)
HOPB1_QA<-select(HOPB1, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                 fDOMNullQF, fDOMGapQF,
                 fDOMSpikeQF, fDOMValidCalQF,
                 fDOMSuspectCalQF, fDOMPersistenceQF,
                 fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                 fDOMFinalQF, fDOMFinalQFSciRvw) 
HOPB1_QA[is.na(HOPB1_QA)] <- 0 #change NAs to 0 
HOPB1_QA <-filter(HOPB1_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                    fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                    fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                    fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2706902 to 649142
HOPB1_QA<-filter(HOPB1_QA, horizontalPosition == "101" )  

#hopb sites half 2 
HOPB2_QA<-select(HOPB2, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                 fDOMNullQF, fDOMGapQF,
                 fDOMSpikeQF, fDOMValidCalQF,
                 fDOMSuspectCalQF, fDOMPersistenceQF,
                 fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                 fDOMFinalQF, fDOMFinalQFSciRvw) 
HOPB2_QA[is.na(HOPB2_QA)] <- 0 #change NAs to 0 
HOPB2_QA <-filter(HOPB2_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                    fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                    fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                    fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2045867 to 514312
HOPB2_QA<-filter(HOPB2_QA, horizontalPosition == "102" ) #102= downstream sensor 

#WALK site
WALK_QA<-select(WALK, siteID, horizontalPosition, startDateTime, fDOM, fDOMRangeQF, fDOMStepQF,
                fDOMNullQF, fDOMGapQF,
                fDOMSpikeQF, fDOMValidCalQF,
                fDOMSuspectCalQF, fDOMPersistenceQF,
                fDOMAlphaQF, fDOMBetaQF, fDOMTempQF, fDOMAbsQF,
                fDOMFinalQF, fDOMFinalQFSciRvw) 
WALK_QA[is.na(WALK_QA)] <- 0 #change NAs to 0 
WALK_QA <-filter(WALK_QA, fDOMRangeQF == 0 & fDOMStepQF == 0 &
                   fDOMNullQF == 0 & fDOMGapQF == 0 & fDOMSpikeQF == 0 & fDOMTempQF == 0 & fDOMAbsQF == 0 &
                   fDOMValidCalQF == 0 & fDOMSuspectCalQF == 0 & fDOMPersistenceQF == 0 & fDOMAlphaQF == 0 &
                   fDOMBetaQF == 0 & fDOMFinalQF == 0 & fDOMFinalQFSciRvw == 0) 
#2131068 to 1065968
WALK_QA<-filter(WALK_QA, horizontalPosition == "101" ) #101= upstream sensor 

#merge all the tables into one 
fDOM_sensor_10sites<- gtools::smartbind(ARIK_QA, BLUE_CUPE_QA, 
                                        COMO_KING_QA, HOPB1_QA, HOPB2_QA, 
                                        MART_QA, MAYF_TOMB_QA, WALK_QA)

#put year, day, month and time in a different column?! 
fDOM_sensor_10sites$startDateTime<-as.POSIXct(fDOM_sensor_10sites$startDateTime,format="%Y-%m-%dT%H:%M:%OS") #extract date and time to a different format 
fDOM_sensor_10sites$DATE<-as.Date(fDOM_sensor_10sites$startDateTime,format="%Y-%m-%d")
fDOM_sensor_10sites$YEAR<-format(fDOM_sensor_10sites$DATE,format="%Y")
fDOM_sensor_10sites$MONTH<-format(fDOM_sensor_10sites$DATE,format="%m")
fDOM_sensor_10sites$DAY<-format(fDOM_sensor_10sites$DATE,format="%d")

fDOM_sensor<-select(fDOM_sensor_10sites, siteID, fDOM, startDateTime, DATE, YEAR, MONTH, DAY)

write.csv(fDOM_sensor, 'Data/sw_fDOM_sensor.csv', row.names = FALSE)

####daily averages of chlorophyll a ####
sensor_chl <- sensor1 %>%
  group_by(siteID, DATE) %>%
  summarise( 
    n=n(),
    mean=mean(chlorophyll, na.rm=TRUE)
  )

# sensor data that isn't NA
nona<-sensor_chl[!with(sensor_chl,is.na(mean)),]

#all sites 
ggplot(data = nona, aes(x=monthYr, y=mean)) +
                 geom_point(aes(colour=siteID)) 

BLUE<-filter(nona, siteID == "BLUE" )
ggplot(data = BLUE, aes(x=DATE, y=mean)) +
        geom_point(aes(colour=siteID)) 

ARIK<-filter(nona, siteID == "ARIK" )
ggplot(data = ARIK, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

COMO<-filter(nona, siteID == "COMO" )
ggplot(data = COMO, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

KING<-filter(nona, siteID == "KING" )
ggplot(data = KING, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

MAYF<-filter(nona, siteID == "MAYF" )
ggplot(data = MAYF, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

TOMB<-filter(nona, siteID == "TOMB" )
ggplot(data = TOMB, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

FLNT<-filter(nona, siteID == "FLNT" )
ggplot(data = FLNT, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

CUPE<-filter(nona, siteID == "CUPE" )
ggplot(data = CUPE, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 

GUIL<-filter(nona, siteID == "GUIL" )
ggplot(data = GUIL, aes(x=monthYr, y=mean)) +
  geom_point(aes(colour=siteID)) 


#### try to merge grab with sensor ####

# having a full join by patient ID
full_df <- full_join(df1, df2, by = "ID")

# select obs within 5 days
result <- full_df %>%
  # including the day difference for your reference
  mutate(Day.diff = abs(as.Date(Start.Date, "%m/%d/%Y") - as.Date(Lab.date, "%m/%d/%Y"))) %>%
  # filtering the data frame to keep the difference within 5 days
  filter(Day.diff <= 5)


chl_sensor<-select(sites, siteID, DATE, MONTH, DAY, monthYr, chlorophyll)
#daily averages 
chl_avg <- chl_sensor %>%
  group_by(siteID, monthYr ) %>%
  summarise( 
    n=n(),
    mean=mean(chlorophyll, na.rm=TRUE)
  )

sw_data$YEAR<-format(sw_data$DATE,format="%Y")
sw_data$MONTH<-format(sw_data$DATE,format="%m")
sw_data$DAY<-format(sw_data$DATE,format="%d")
surface_grab_2019<-filter(sw_data, YEAR == "2019")
surface_grab_2019$monthYr = format(as.Date(surface_grab_2019$DATE), "%m.%d") 

combined<-left_join(chl_avg, surface_grab_2019,  by = c('siteID', 'monthYr'))
#write.csv(combined, 'Data/ten_sites.csv', row.names = FALSE)

#select out 12 variables + chlorophyll ### removed Conductivity!! 
variables<-select(combined, siteID, monthYr, mean,
                  dissolvedOrganicCarbon, waterIron, 
                  waterNitrateAndNitriteN, waterManganese, waterChlorine, 
                  waterFluorine, waterPotassium, waterSodium, 
                  waterCalcium, dissolvedInorganicCarbon, waterMagnesium)

#change the name of "mean" to cholophyll 
names(variables)[3]<-paste("chlorophyll_mean")

#select only rows with all of the data points 
nona<-variables[!with(variables,is.na(externalConductance)),]


#chla
my_chla <- variables %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(chlorophyll_mean, na.rm=TRUE),
    sd=sd(chlorophyll_mean, na.rm=TRUE)
  )

my_chla$siteID <-ordered(my_chla$siteID, 
                         levels=c("MAYF", "COMO", "MART", 
                                  "HOPB", "TOMB", "WALK", "CUPE",
                                  "ARIK", "KING", "BLUE"))

print(ggplot(my_chla) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("chla") + theme(axis.text.x = element_text(angle = 90))


