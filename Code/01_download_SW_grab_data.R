#### GRAB SAMPLES #### 
#All 27 sites: Grab samples of surface water chemistry including general chemistry (DOC), anions, cations, and nutrients. streams 26 times per year
nutrients_grab <- loadByProduct(dpID="DP1.20093.001", 
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

#get chl-a data from surface water grab smaples 
chl_grab <- loadByProduct(dpID="DP1.20163.001", 
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

#get isotope data from surface water grab smaples 
iso_grab <- loadByProduct(dpID="DP1.20206.001", 
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

#get nitrate 20033.001
nitrate_grab <- loadByProduct(dpID="DP1.20033.001", 
                              site=c(
                                'HOPB','POSE','KING', 'WALK','LECO','MAYF','PRIN',
                                'BLDE','COMO','MART', 'BIGC','CARI', 'WLOU', "FLNT", 
                                "MCDI", 'LEWI', "BLUE", "TECR", "REDB", "SYCA", 
                                "MCRA", "OKSR", "ARIK", "GUIL", "CUPE", "TOMB", "BLWA"
                              ),
                              startdate="2012-01", 
                              enddate="2019-09", 
                              package="basic", #basic will just give you concentrations, #expanded will give you flags 
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
grab_chl<-read.csv('Data/surface_water_chla_grab.csv', header = TRUE) 
grab_discharge<-read.csv('Data/grab_discharge.csv', header = TRUE)

### Extracts date only to a new column 
grab_dat$DATE<-as.Date(grab_dat$collectDate,format="%Y-%m-%d")
grab_discharge$DATE<-as.Date(grab_discharge$collectDate,format="%Y-%m-%d")

### select out only needed columns from these tables and join all tables into one 
grab_dat<-subset(grab_dat, select = -c(collectDate, laboratoryName, coolerTemp, #this function you list columns to remove
                                       sampleCondition, remarks, shipmentWarmQF, externalLabDataQF, 
                                       receivedBy, shipmentCondition, shipmentLateQF)) 
grab_info<-select(grab_info, decimalLatitude, decimalLongitude, elevation, parentSampleID, dissolvedOxygen, specificConductance, waterTemp)
grab_discharge<-select(grab_discharge, siteID, DATE, totalDischarge)
grab_chl<-select(grab_chl, uid, siteID, collectDate, namedLocation, analyte, analyteConcentration) # function lists columns that I want to save
## need to rotate table with chl data and select out chla
grab_chl2<-tidyr::pivot_wider(grab_chl, id_cols = c(uid, siteID, collectDate, namedLocation), names_from = analyte,
                              values_from = analyteConcentration)
data.table::setnames(grab_chl2, "total chlorophyll a", "Tchla")
data.table::setnames(grab_chl2, "chlorophyll a", "chla")
grab_chl2$DATE<-as.Date(grab_chl2$collectDate,format="%Y-%m-%d")
grab_chl2<-select(grab_chl2, siteID, DATE, namedLocation, Tchla, chla) 
grab_chl3<-grab_chl2[!with(grab_chl2,is.na(Tchla)& is.na(chla)),] #get rid of rows with NAs in both columns
grab_chl3$Tchla[is.na(grab_chl3$Tchla)] <- grab_chl3$chla[is.na(grab_chl3$Tchla)] #fill in the NA of Tchla with the values from Chla

#join tables into one big table #want to mach observances by site and date
sw_data<-left_join(grab_dat, grab_info, by = c("sampleID" = "parentSampleID")) %>%
  left_join(grab_discharge, by = c("siteID", "DATE")) %>%
  left_join(grab_chl3,  by = c("siteID", "DATE", "namedLocation"))

#### exploring variation within a site across years (temporal) #### 
CARI<- filter(sw_data, siteID == "CARI")
plot(CARI$DATE, CARI$totalDischarge)  

#all sites 
ggplot(data = sw_data, aes(x=DATE, y=totalDischarge)) + geom_point(aes(colour=siteID))

#select out only 2016-current #seems too messy to actually use 
sw_data$YEAR<-lubridate::year(sw_data$DATE)
recent_years<-filter(sw_data, YEAR >= 2016)
ggplot(data = recent_years, aes(x=DATE, y=totalDischarge)) + geom_point(aes(colour=siteID))


#### Graphing for exploring variation across sites (spatial) #### 
#graph bar plot with error bars 
#DO
my_DO <- sw_data %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(dissolvedOxygen, na.rm=TRUE),
    sd=sd(dissolvedOxygen, na.rm=TRUE)
  )

print(ggplot(my_DO) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("DO") + theme(axis.text.x = element_text(angle = 90))


#nitrate/nitrite
my_NN <- sw_data %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterNitrateAndNitriteN, na.rm=TRUE),
    sd=sd(waterNitrateAndNitriteN, na.rm=TRUE)
  )

print(ggplot(my_NN) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("NitrateNitrite") + theme(axis.text.x = element_text(angle = 90))


#discharge
my_dis <- sw_data %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(totalDischarge, na.rm=TRUE),
    sd=sd(totalDischarge, na.rm=TRUE)
  )

print(ggplot(my_dis) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("discharge") + theme(axis.text.x = element_text(angle = 90))



#chla
my_chla <- sw_data %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Tchla, na.rm=TRUE),
    sd=sd(Tchla, na.rm=TRUE)
  )

print(ggplot(my_chla) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("chla") + theme(axis.text.x = element_text(angle = 90))


#pH 
my_pH <- sw_data %>%
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

#### Extracts datetime and converts datetime string to POSIXct ####
wqvalues101$startDateTime<-as.POSIXct(wqvalues101$startDateTime,format="%Y-%m-%dT%H:%M:%OS")
wqvalues102$startDateTime<-as.POSIXct(wqvalues102$startDateTime,format="%Y-%m-%dT%H:%M:%OS")

### Plots data ###
plot(wqvalues101$startDateTime,wqvalues101$dissolvedOxygen,type="l",col="blue",main="ARIK DO",xlab="Date",ylab="DO (mg/L)")
lines(wqvalues102$startDateTime,wqvalues102$dissolvedOxygen,col="red")
grid(nx=NULL,ny=NULL, col="lightgray",lty="dotted",lwd=par("lwd"), equilogs=TRUE)
legend("bottomleft",legend=c("Upstream", "Downstream"),
       col=c("blue", "red"), lty=1:1, cex=0.8)

