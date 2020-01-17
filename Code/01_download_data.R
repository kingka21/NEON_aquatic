#### Script for downloading and organizing NEON datasets ##### 
### code written by Katelyn King 16-OCT-2019 
### code adapted from Bobby Hensley ### 

#### load libraries ####
# install neonUtilities if you have not already 
install.packages("neonUtilities")
# load neonUtilities
library(neonUtilities)
?loadByProduct
?getDatatable
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
##      20033.001 = nitrate in surface water
##      20092.001 = groundwater 
##      20163.001 = chl
##      20206.001 = isotopes in surface water 
##      20276.001 = isotopes in groundwater
##      20048.001 = manual discharge measures (about 26/year/site)

#Water quality Sensor data EXAMPLE! 
#In situ sensor-based specific conductivity, concentration of chlorophyll a, dissolved oxygen content, fDOM concentration (fluorescent dissolved material), pH, and turbidity 
#available as one-, five-, and thirty-minute averages in surface water 
waterchem_sensor <- loadByProduct(dpID="DP1.20288.001", 
                           site=c("MOAB","ONAQ"),
                           startdate="2018-05",  #year and month
                           enddate="2018-08" ,
                           check.size = T , 
                           avg= "30") # 'all', or the averaging interval to download, in minutes. Only applicable to sensor data.

#All 27 sites: Grab samples of surface water chemistry including general chemistry (DOC), anions, cations, and nutrients. streams 26 times per year
nutrients_grab <- loadByProduct(dpID="DP1.20093.001", 
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

#get discharge data which was taken manually  
discharge_grab <- loadByProduct(dpID="DP1.20048.001", 
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


# Turn data into a dataframe (can I use the get datatable function?)
for(i in 1:length(nutrients_grab)) {assign(names(nutrients_grab)[i], nutrients_grab[[i]])}   #calls the table waq_instances
grab_chem_dat<-as.data.frame(swc_externalLabData) #table that has water chem
grab_info_dat<-as.data.frame(swc_fieldSuperParent) #table that has lat/lon and elevation, DO, waterTemp, maxDepth

for(i in 1:length(chl_grab)) {assign(names(chl_grab)[i], chl_grab[[i]])}   #calls the table 
grab_chl_dat<-as.data.frame(alg_algaeExternalLabDataPerSample) #table that has total chlorophyll a and chlorophyll a (dont know the diff right now)

for(i in 1:length(discharge_grab)) {assign(names(discharge_grab)[i], discharge_grab[[i]])}   #calls the table waq_instances
grab_discharge_dat<-as.data.frame(dsc_fieldData) #table that has total discharge in cubic meters per sec 

write.csv(grab_chem_dat, 'Data/surface_water_grab.csv', row.names = FALSE)
write.csv(grab_info_dat, 'Data/grab_info.csv', row.names = FALSE)
write.csv(grab_chl_dat, 'Data/surface_water_chla_grab.csv', row.names = FALSE)
write.csv(grab_discharge_dat, 'Data/grab_discharge.csv', row.names = FALSE)


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
test <-left_join(grab_dat, grab_discharge, by = c("siteID", "DATE")) #want to mach observances by site and date


#### exploring variation within a site across years (temporal) #### 
HOPB<- filter(grab_dat, siteID == "HOPB")
plot(HOPB$DATE, HOPB$waterTotalOrganicCarbon, type='l')  

#all sites 
ggplot(data = grab_dat, aes(x=DATE, y=waterTotalOrganicCarbon)) + geom_line(aes(colour=siteID))

#select out only 2016-current #seems to messy to actually use 
grab_dat$YEAR<-lubridate::year(grab_dat$DATE)
recent_years<-filter(grab_dat, YEAR >= 2016)
ggplot(data = recent_years, aes(x=DATE, y=waterTotalOrganicCarbon)) + geom_line(aes(colour=siteID))


#### Graphing for exploring variation across sites (spatial) #### 
#graph bar plot with error bars 
#pH 
my_pH <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(pH, na.rm=TRUE),
    sd=sd(pH, na.rm=TRUE)
  )

print(ggplot(my_pH) +
  geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
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

ggplot(my_con) +
  geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
  geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("conductivity")

# "totalSuspendedSolids"
my_tss <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(totalSuspendedSolids, na.rm=TRUE),
    sd=sd(totalSuspendedSolids, na.rm=TRUE)
  )

print(ggplot(my_tss) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("tss") + theme(axis.text.x = element_text(angle = 90))

#"waterNitriteN"  
my_N <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterNitriteN, na.rm=TRUE),
    sd=sd(waterNitriteN, na.rm=TRUE)
  )

print(ggplot(my_N) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("nitrite") + theme(axis.text.x = element_text(angle = 90))

#"waterTotalOrganicCarbon" 
my_orgcarbon <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterTotalOrganicCarbon, na.rm=TRUE),
    sd=sd(waterTotalOrganicCarbon, na.rm=TRUE)
  )

print(ggplot(my_orgcarbon) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("orgcarbon") + theme(axis.text.x = element_text(angle = 90))


#"dissolvedInorganicCarbon"  
my_inorgcarbon <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(dissolvedInorganicCarbon, na.rm=TRUE),
    sd=sd(dissolvedInorganicCarbon, na.rm=TRUE)
  )

print(ggplot(my_inorgcarbon) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("inorganic carbon") + theme(axis.text.x = element_text(angle = 90))

# "dissolvedOrganicCarbon"      
my_diss_org_C <- grab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(dissolvedOrganicCarbon, na.rm=TRUE),
    sd=sd(dissolvedOrganicCarbon, na.rm=TRUE)
  )

print(ggplot(my_diss_org_C) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("dissolved organic carbon") + theme(axis.text.x = element_text(angle = 90))

#"uvAbsorbance250" 
#waterMagnesium


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
