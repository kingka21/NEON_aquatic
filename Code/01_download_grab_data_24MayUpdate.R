#### Script for downloading and organizing NEON datasets ##### 
### code written by Katelyn King 16-OCT-2019 
### code adapted from Bobby Hensley ### 
###Modified by Jennifer Edmonds 24May2020####

#### load libraries ####
# install neonUtilities if you have not already 
# C:\Users\INBRE\Documents\NEON_aquatic
install.packages("neonUtilities")
# load neonUtilities
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


###############################################################################################################
#All 27 sites: Grab samples of surface water chemistry including general chemistry (DOC), anions, cations, and nutrients. streams 26 times per year
nutrients_SWgrab_allsites <- loadByProduct(dpID="DP1.20093.001", 
                           site=c(
                             'HOPB','POSE','KING', 'WALK','LECO','MAYF','PRIN',
                             'BLDE','COMO','MART', 'BIGC','CARI', 'WLOU', "FLNT", 
                             "MCDI", 'LEWI', "BLUE", "TECR", "REDB", "SYCA", 
                             "MCRA", "OKSR", "ARIK", "GUIL", "CUPE", "TOMB", "BLWA"
                           ),
                           startdate="2017-01", 
                           enddate="2019-12", 
                           package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                           check.size = F)  ### check the size of the file before you download it 

for(i in 1:length(nutrients_SWgrab_allsites)) {assign(names(nutrients_SWgrab_allsites)[i], nutrients_SWgrab_allsites[[i]])}   #calls the table waq_instances
SWgrab_chem_dat_allsites<-as.data.frame(swc_externalLabDataByAnalyte) #table that has water chem
SWgrab_info_dat_allsites<-as.data.frame(swc_fieldSuperParent) #table that has lat/lon and elevation, DO, waterTemp, maxDepth
write.csv(SWgrab_chem_dat_allsites, 'Data/SW_grab_allsites_unaltered.csv', row.names = FALSE)#contains 29,641 observations

#first step is to remove lines where sampleID ends with .2 or .3#
SWgrab_chem_dat_allsites_REMOVE<-filter(SWgrab_chem_dat_allsites, !grepl("FIL.3$", sampleID))
SWgrab_chem_dat_allsites_REMOVE<-filter(SWgrab_chem_dat_allsites_REMOVE, !grepl("FIL.2$", sampleID))
SWgrab_chem_dat_allsites_REMOVE<-filter(SWgrab_chem_dat_allsites_REMOVE, !grepl("RAW.3$", sampleID))
SWgrab_chem_dat_allsitesL_REMOVE<-filter(SWgrab_chem_dat_allsites_REMOVE, !grepl("RAW.2$", sampleID))
SWgrab_chem_dat_allsites_REMOVE<-filter(SWgrab_chem_dat_allsites_REMOVE, !grepl("PCN.3$", sampleID))
SWgrab_chem_dat_allsites_REMOVE<-filter(SWgrab_chem_dat_allsites_REMOVE, !grepl("PCN.2$", sampleID))
write.csv(SWgrab_chem_dat_allsites_REMOVE, 'Data/SW_grab_allsites_REMOVE.csv', row.names = FALSE)#contains 26,499 observations

SWgrab_chem_dat_allsites_REMOVE$sampleID<-gsub(".FIL" , "" , SWgrab_chem_dat_allsites_REMOVE$sampleID)
SWgrab_chem_dat_allsites_REMOVE$sampleID<-gsub(".RAW" , "" , SWgrab_chem_dat_allsites_REMOVE$sampleID)
SWgrab_chem_dat_allsites_REMOVE$sampleID<-gsub(".PCN" , "" , SWgrab_chem_dat_allsites_REMOVE$sampleID)
#check file to make sure removal of extensions was completed#
write.csv(SWgrab_chem_dat_allsites_REMOVE, 'Data/surface_water_grab_allsites.csv', row.names = FALSE)

#Remove flagged data#
SWgrab_chem_dat_allsites_QF<-filter(SWgrab_chem_dat_allsites_REMOVE, shipmentWarmQF == 0)

#Setting negative values for nutrients to 0#
SWgrab_chem_dat_allsites_QF[SWgrab_chem_dat_allsites_QF <0] <- 0 #change negative values to 0 
write.csv(SWgrab_chem_dat_allsites_QF, 'Data/surface_water_grab_allsites_QF.csv', row.names = FALSE)#23,121 observations

#generate new table using pivot data for water chemistry parameters#
SWgrab_chem_dat_allsites_PIVOT<-pivot_wider(SWgrab_chem_dat_allsites_QF, id_cols= c(siteID, sampleID, collectDate), names_from=analyte, values_from=analyteConcentration)

#check file
write.csv(SWgrab_chem_dat_allsites_PIVOT, 'Data/surface_water_grab_allsites_PIVOT.csv', row.names = FALSE)



########################################################################################################

#Focus sites for surfacewater, final 10 sites (+3 replacements), 3 years of data 2017-2019
SWgrab_FINALsubset <- loadByProduct(dpID="DP1.20093.001", 
                                              site=c('KING','MAYF','COMO',"HOPB","WALK", "MART", "BLUE", "ARIK", "TOMB", "CUPE",'GUIL','POSE','BLWA'),
                                              startdate="2017-01",  #year and month
                                              enddate="2019-12", 
                                              package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                                              check.size = F)  ### check the size of the file before you download it 


for(i in 1:length(SWgrab_FINALsubset)) {assign(names(SWgrab_FINALsubset)[i], SWgrab_FINALsubset[[i]])}   #calls the table waq_instances
SWgrab_chem_dat_subsetFINAL<-as.data.frame(swc_externalLabDataByAnalyte) #table that has water chem
SWgrab_info_dat_subsetFINAL<-as.data.frame(swc_fieldSuperParent) #table that has lat/lon and elevation, DO, waterTemp, maxDepth
write.csv(SWgrab_chem_dat_subsetFINAL, 'Data/surface_water_grab_subsetFINAL_NoModification.csv', row.names = FALSE)

#32,521 observations (13 sites) in initial download of surface water chemistry, grab samples
#first step is to remove lines where sampleID ends with .2 or .3, remaining observations=18,693#
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL, !grepl("FIL.3$", sampleID))
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, !grepl("FIL.2$", sampleID))
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, !grepl("RAW.3$", sampleID))
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, !grepl("RAW.2$", sampleID))
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, !grepl("PCN.3$", sampleID))
SWgrab_chem_dat_subsetFINAL_REMOVE<-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, !grepl("PCN.2$", sampleID))

#28,158 observations, check file to make sure removal of lines with FIL.3, FIL.2, RAW.3, RAW.2, PCN.3, and PCN.2 was completed#
write.csv(SWgrab_chem_dat_subsetFINAL_REMOVE, 'Data/surface_water_grab_subsetFINAL_REMOVED.csv', row.names = FALSE)

SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID<-gsub(".FIL" , "" , SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID)
SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID<-gsub(".RAW" , "" , SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID)
SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID<-gsub(".PCN" , "" , SWgrab_chem_dat_subsetFINAL_REMOVE$sampleID)
#check file to make sure removal of extensions was completed#
write.csv(SWgrab_chem_dat_subsetFINAL_REMOVE, 'Data/surface_water_grab_subsetFINAL.csv', row.names = FALSE)

#remove all values with shipmentWarmQF=1, reduced to 16,721#
SWgrab_chem_dat_subsetFINAL_QF <-filter(SWgrab_chem_dat_subsetFINAL_REMOVE, shipmentWarmQF == 0)

#Setting negative values for nutrients to 0#
SWgrab_chem_dat_subsetFINAL_QF[SWgrab_chem_dat_subsetFINAL_QF <0] <- 0 #change negative values to 0 

#25,336 observations, check file to make sure removal of flagged data and zeros was completed#
write.csv(SWgrab_chem_dat_subsetFINAL_QF, 'Data/surface_water_grab_subsetFINAL_QF.csv', row.names = FALSE)

#generate new table using pivot data for water chemistry parameters#
SWgrab_chem_dat_subsetFINAL_PIVOT<-pivot_wider(SWgrab_chem_dat_subsetFINAL_QF, id_cols= c(siteID, sampleID, collectDate), names_from=analyte, values_from=analyteConcentration)

#check file to make sure pivot was executed properly"
write.csv(SWgrab_chem_dat_subsetFINAL_PIVOT, 'Data/surface_water_grab_QAQC.csv', row.names = FALSE)

################################################################################################################

#Focus sites for groundwater, final 10 sites (+3 replacements), 3 years of data 2017-2019
GWgrab_FINALsubset <- loadByProduct(dpID="DP1.20092.001", 
                                               site=c('KING','MAYF','COMO',"HOPB","WALK", "MART", "BLUE", "ARIK", "TOMB", "CUPE", "GUIL", "POSE", "BLWA"),
                                               startdate="2017-01", 
                                               enddate="2019-12", 
                                               package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                                               check.size = F)  ### check the size of the file before you download it 

for(i in 1:length(GWgrab_FINALsubset)) {assign(names(GWgrab_FINALsubset)[i], GWgrab_FINALsubset[[i]])}   #calls the table waq_instances
GWgrab_chem_dat_subsetFINAL<-as.data.frame(gwc_externalLabDataByAnalyte) #table that has water chem in groundwater
GWgrab_info_dat_subsetFINAL<-as.data.frame(gwc_fieldSuperParent) #table that has lat/lon and elevation, DO, waterTemp, maxDepth

#6,390 observations in the initial downloaded data set for groundwater grab samples, water chemistry#
write.csv(GWgrab_chem_dat_subsetFINAL, 'Data/groundwater_grab_subsetFINAL_NOMODIFICATION.csv', row.names = FALSE)
#There are no replicate observations in the groundwater dataset
#Remove extensions on samplesID#
GWgrab_chem_dat_subsetFINAL$sampleID<-gsub(".FIL" , "" , GWgrab_chem_dat_subsetFINAL$sampleID)
GWgrab_chem_dat_subsetFINAL$sampleID<-gsub(".RAW" , "" , GWgrab_chem_dat_subsetFINAL$sampleID)
GWgrab_chem_dat_subsetFINAL$sampleID<-gsub(".PCN" , "" , GWgrab_chem_dat_subsetFINAL$sampleID)
#check file to make sure removal of extensions was completed#
write.csv(GWgrab_chem_dat_subsetFINAL, 'Data/groundwater_grab_subsetFINAL_REMOVED.csv', row.names = FALSE)

# following flagging removal#
GWgrab_chem_dat_subsetFINAL_QF <-filter(GWgrab_chem_dat_subsetFINAL, shipmentWarmQF == 0)

#Setting negative values for nutrients to 0#
GWgrab_chem_dat_subsetFINAL_QF[GWgrab_chem_dat_subsetFINAL_QF <0] <- 0 #change negative values to 0 

#5,778 observations following removal of flagged data#
write.csv(GWgrab_chem_dat_subsetFINAL_QF_noNA, 'Data/groundwater_grab_subsetFINAL_QF.csv', row.names = FALSE)

#generate new table using pivot data for groundwater chemistry parameters#
GWgrab_chem_dat_subsetFINAL_PIVOT<-pivot_wider(GWgrab_chem_dat_subsetFINAL_QF, id_cols= c(siteID, collectDate, sampleID), names_from=analyte, values_from=analyteConcentration)

#write datafrom from memory to csv file"
write.csv(GWgrab_chem_dat_subsetFINAL_PIVOT, 'Data/groundwater_grab_QAQC.csv', row.names = FALSE)

############################################################################################################################

#get discharge data for focus 10 focus sites (+3 replacement sites), taken manually 
discharge_grab_subsetFINAL <- loadByProduct(dpID="DP1.20048.001", 
                                site=c( 'KING','MAYF','COMO',"HOPB","WALK", "MART", "BLUE", "ARIK", "TOMB", "CUPE","GUIL","POSE","BLWA"),
                                startdate="2017-01", 
                                enddate="2019-12", 
                                package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                                check.size = F)  ### check the size of the file before you download it 

for(i in 1:length(discharge_grab_subsetFINAL)) {assign(names(discharge_grab_subsetFINAL)[i], discharge_grab_subsetFINAL[[i]])}   #calls the table waq_instances
discharge_subsetFINAL<-as.data.frame(dsc_individualFieldData) #table that has all measurements for discharge
stage_subsetFINAL<-as.data.frame(dsc_fieldData)#table tht has stage readings from field

#These dataframes are Level 1 data, need to convert to Level 2 data using conv.calc.Q function before using in further analyses# 

##############     END OF CODE EDITED 24May2020   ###########################################################
#############################################################################################################

#get chl-a data from surface water grab smaples 
chl_grab <- loadByProduct(dpID="DP1.20163.001", 
                          site=c(
                            'HOPB','POSE','KING', 'WALK','LECO','MAYF','PRIN',
                            'BLDE','COMO','MART', 'BIGC','CARI', 'WLOU', "FLNT", 
                            "MCDI", 'LEWI', "BLUE", "TECR", "REDB", "SYCA", 
                            "MCRA", "OKSR", "ARIK", "GUIL", "CUPE", "TOMB", "BLWA"
                          ),
                          startdate="2012-01", 
                          enddate="2019-09", 
                          package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                          check.size = F)  ### check the size of the file before you download it 

#get isotope data from surface water grab smaples 
iso_grab <- loadByProduct(dpID="DP1.20206.001", 
                          site=c(
                            'HOPB','POSE','KING', 'WALK','LECO','MAYF','PRIN',
                            'BLDE','COMO','MART', 'BIGC','CARI', 'WLOU', "FLNT", 
                            "MCDI", 'LEWI', "BLUE", "TECR", "REDB", "SYCA", 
                            "MCRA", "OKSR", "ARIK", "GUIL", "CUPE", "TOMB", "BLWA"
                          ),
                          startdate="2012-01", 
                          enddate="2019-09", 
                          package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                          check.size = F)  ### check the size of the file before you download it 

write.csv(grab_chl_dat, 'Data/surface_water_chla_grab.csv', row.names = FALSE)

#############################################################
#### read in saved data and do some exploratory analysis #### 
#############################################################

#read in chem data and info for each site to do analysis or graph 
#note sampleID in grab_dat and parentSampleID in grab_info are the identifiers to join
grab_dat<-read.csv('Data/surface_water_grab.csv', header = TRUE)
grab_info<-read.csv('Data/grab_info.csv', header = TRUE)
grab_chl<-read.csv('Data/surface_water_chla_grab.csv', header = TRUE) ##maybe named location can help join?
grab_discharge<-read.csv('Data/grab_discharge.csv', header = TRUE)

### Extracts date only to a new column 
grab_dat$DATE<-as.Date(grab_dat$collectDate,format="%Y-%m-%d")
grab_discharge$DATE<-as.Date(grab_discharge$collectDate,format="%Y-%m-%d")

### select out only needed columns from these tables and join all tables into one 
grab_dat<-subset(grab_dat, select = -c(namedLocation, collectDate, laboratoryName, coolerTemp, #this function you list columns to remove
                                       sampleCondition, remarks, shipmentWarmQF, externalLabDataQF, 
                                       receivedBy, shipmentCondition, shipmentLateQF)) 
grab_info<-select(grab_info, decimalLatitude, decimalLongitude, elevation, parentSampleID, dissolvedOxygen, specificConductance, waterTemp)
grab_discharge<-select(grab_discharge, siteID, DATE, totalDischarge)
grab_chl<-select(grab_chl, siteID, collectDate, analyte, analyteConcentration) # function lists columns that I want to save
## need to rotate table with chl data and select out chla
grab_chl2<-tidyr::pivot_wider(grab_chl, id_cols = c(siteID,collectDate), names_from = analyte,
            values_from = analyteConcentration)  ### need to figure out which site along the reach to choose

#join tables into one big table 
sw_data<-left_join(grab_dat, grab_info, by = c("sampleID" = "parentSampleID")) 
test <-left_join(grab_dat, grab_discharge, by = c("siteID", "DATE")) #want to match observances by site and date


#### exploring variation within a site across years (temporal) #### 
HOPB<- filter(grab_dat, siteID == "HOPB")
plot(HOPB$DATE, HOPB$waterTotalOrganicCarbon, type='l')  

#all sites 
ggplot(data = grab_dat, aes(x=DATE, y=waterTotalOrganicCarbon)) + geom_line(aes(colour=siteID))

#select out only 2016-current #seems to messy to actually use 
grab_dat$YEAR<-lubridate::year(grab_dat$DATE)
recent_years<-filter(grab_dat, YEAR >= 2016)
ggplot(data = recent_years, aes(x=DATE, y=waterTotalOrganicCarbon)) + geom_line(aes(colour=siteID))


#### Graphing 10 sites of interest for exploring variation across sites (spatial) #### 
#graph bar plot with error bars 
grab_dat<-read.csv('Data/surface_water_grab_subsetCOMBINED.csv', header = TRUE)
grab_info<-read.csv('Data/grab_info_subsetCOMBINED.csv', header = TRUE)

#pH 
my_pH <- grab_dat %>%
  group_by (siteID) %>%
  summarise( 
    n=n(),
    mean=mean(pH, na.rm=TRUE),
    sd=sd(pH, na.rm=TRUE)
  )
my_pH$siteID<-ordered(my_pH$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))



print(ggplot(my_pH) +
  geom_bar( aes(x=siteID, y=mean), stat= "identity", fill="skyblue") +
  geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("pH") + theme(axis.text.x = element_text(angle = 90)) 

# conductivity
my_con <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(externalConductance, na.rm=TRUE),
    sd=sd(externalConductance, na.rm=TRUE)
  )
my_con$siteID<-ordered(my_con$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


ggplot(my_con) +
  geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
  geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("conductivity")

# "waterChlorine"
my_Cl <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterChlorine, na.rm=TRUE),
    sd=sd(waterChlorine, na.rm=TRUE)
  )
my_Cl$siteID<-ordered(my_Cl$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Cl) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Chlorine") + theme(axis.text.x = element_text(angle = 90))

#"waterIron"  
my_iron <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterIron, na.rm=TRUE),
    sd=sd(waterIron, na.rm=TRUE)
  )
my_iron$siteID<-ordered(my_iron$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_iron) + 
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") + 
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
        ylab("iron") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.07, 0.72)

#"dissolvedOrganicCarbon" 
my_DOC <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(dissolvedOrganicCarbon, na.rm=TRUE),
    sd=sd(dissolvedOrganicCarbon, na.rm=TRUE)
  )
my_DOC$siteID<-ordered(my_DOC$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))

print(ggplot(my_DOC) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("DOC") + theme(axis.text.x = element_text(angle = 90))


#"waterManganese"  
my_Mn <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterManganese, na.rm=TRUE),
    sd=sd(waterManganese, na.rm=TRUE)
  )
my_Mn$siteID<-ordered(my_Mn$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Mn) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Manganese") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.01, 0.04)


# "waterPotassium"      
my_K <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterPotassium, na.rm=TRUE),
    sd=sd(waterPotassium, na.rm=TRUE)
  )
my_K$siteID<-ordered(my_K$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_K) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Potassium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.3, 3)

# "waterSodium"      
my_Na <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterSodium, na.rm=TRUE),
    sd=sd(waterSodium, na.rm=TRUE)
  )
my_Na$siteID<-ordered(my_Na$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Na) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Sodium") + theme(axis.text.x = element_text(angle = 90)) + ylim(0, 18)


# "waterFluorine"      
my_F <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterFluorine, na.rm=TRUE),
    sd=sd(waterFluorine, na.rm=TRUE)
  )
my_F$siteID<-ordered(my_F$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_F) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Fluorine") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.07, 0.75)

# "waterCalcium"      
my_Ca <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterCalcium, na.rm=TRUE),
    sd=sd(waterCalcium, na.rm=TRUE)
  )
my_Ca$siteID<-ordered(my_Ca$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Ca) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Calcium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.7, 116)

# "dissolvedInorganicCarbon"      
my_DIC <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(dissolvedInorganicCarbon, na.rm=TRUE),
    sd=sd(dissolvedInorganicCarbon, na.rm=TRUE)
  )
my_DIC$siteID<-ordered(my_DIC$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_DIC) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("DIC") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.7, 100)

# "waterMagnesium"      
my_Mg <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterMagnesium, na.rm=TRUE),
    sd=sd(waterMagnesium, na.rm=TRUE)
  )
my_Mg$siteID<-ordered(my_Mg$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Mg) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Magnesium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.7, 41)

# "waterNitrateAndNitriteN"      
my_NitrateNitrite <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterNitrateAndNitriteN, na.rm=TRUE),
    sd=sd(waterNitrateAndNitriteN, na.rm=TRUE)
  )
my_NitrateNitrite$siteID<-ordered(my_NitrateNitrite$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' , 'WALK' , 'CUPE', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_NitrateNitrite) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Nitrate+Nitrite") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.04, 0.46)

##### other tips from Bobby 
#### Seperates measurements by location 
##  Upstream = 101              Downstream = 102
##  US overhanging = 111        DS overhanging = 112
##  Lake/river buoy = 103
##  Lake inlet = 130            Lake outlet = 140
##  NOTE: SUNA and fDOM only at downstream stream station or lake/river buoy
##  NOTE: function removes values, so enter term for locations you DON'T want
wqvalues101<-wqvalues[(wqvalues$horizontalPosition=="101"),]
wqvalues102<-wqvalues[(wqvalues$horizontalPosition=="102"),]

### Plots data ###
plot(wqvalues101$startDateTime,wqvalues101$dissolvedOxygen,type="l",col="blue",main="ARIK DO",xlab="Date",ylab="DO (mg/L)")
lines(wqvalues102$startDateTime,wqvalues102$dissolvedOxygen,col="red")
grid(nx=NULL,ny=NULL, col="lightgray",lty="dotted",lwd=par("lwd"), equilogs=TRUE)
legend("bottomleft",legend=c("Upstream", "Downstream"),
       col=c("blue", "red"), lty=1:1, cex=0.8)
