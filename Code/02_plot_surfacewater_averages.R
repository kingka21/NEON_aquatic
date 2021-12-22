
#### Graphing for exploring variation across sites (spatial) #### 
##written by Katelyn King and modified by Jennifer Edmonds 
#graph bar plot with error bars 
SWgrab_dat<-read.csv('Data/surface_water_grab_QCQC_ph.csv', header = TRUE)


#####pH 
my_pH <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise(
            n=n(),
            mean=mean(pH, na.rm=TRUE),
            sd=sd(pH, na.rm=TRUE)
  )

print(ggplot(my_pH) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("pH") + theme(axis.text.x = element_text(angle = 90)) + ylim (0, 8.4)

#####conductivity
my_con <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(conductivity, na.rm=TRUE),
    sd=sd(conductivity, na.rm=TRUE)
  )

ggplot(my_con) +
  geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
  geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("conductivity")



######NO2.NO3...N      
my_NitrateNitrite <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(NO3.NO2...N, na.rm=TRUE),
    sd=sd(NO3.NO2...N, na.rm=TRUE)
  )


print(ggplot(my_NitrateNitrite) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Nitrate+Nitrite") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.03, 0.941)


#####"TN"  
my_TN <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise(
    n=n(),
    mean=mean(TN, na.rm=TRUE),
    sd=sd(TN, na.rm=TRUE)
  )

print(ggplot(my_TN) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("TN") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.05, 1)


#####dissolvedInorganicCarbon     
my_DIC <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(DIC, na.rm=TRUE),
    sd=sd(DIC, na.rm=TRUE)
  )


print(ggplot(my_DIC) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("DIC") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.7, 105)


#####Calcium      
my_Ca <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Ca, na.rm=TRUE),
    sd=sd(Ca, na.rm=TRUE)
  )


print(ggplot(my_Ca) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Calcium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-2.97, 83)


#####Magnesium     
my_Mg <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Mg, na.rm=TRUE),
    sd=sd(Mg, na.rm=TRUE)
  )


print(ggplot(my_Mg) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Magnesium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-1.7, 41)


#####Potassium      
my_K <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(K, na.rm=TRUE),
    sd=sd(K, na.rm=TRUE)
  )

print(ggplot(my_K) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Potassium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.05, 2.03)


#####Na      
my_Na <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Na, na.rm=TRUE),
    sd=sd(Na, na.rm=TRUE)
  )

print(ggplot(my_Na) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Sodium") + theme(axis.text.x = element_text(angle = 90))


#####Manganese  
my_Mn <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Mn, na.rm=TRUE),
    sd=sd(Mn, na.rm=TRUE)
  )


print(ggplot(my_Mn) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Manganese") + theme(axis.text.x = element_text(angle = 90))  + ylim(-0.005, 0.018)



#####dissolvedOrganicCarbon 
my_disorgcarbon <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(DOC, na.rm=TRUE),
    sd=sd(DOC, na.rm=TRUE)
  )

print(ggplot(my_disorgcarbon) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("DOC") + theme(axis.text.x = element_text(angle = 90))


#####Chlorine
my_Cl <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Cl, na.rm=TRUE),
    sd=sd(Cl, na.rm=TRUE)
  )


print(ggplot(my_Cl) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Chlorine") + theme(axis.text.x = element_text(angle = 90))


#####Iron
my_iron <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Fe, na.rm=TRUE),
    sd=sd(Fe, na.rm=TRUE)
  )

print(ggplot(my_iron) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("iron") + theme(axis.text.x = element_text(angle = 90))




#FOLLOWING PARAMETERS NOT USED IN PCA#
# "waterFluorine"      
my_F <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterFluorine, na.rm=TRUE),
    sd=sd(waterFluorine, na.rm=TRUE)
  )


print(ggplot(my_F) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Fluorine") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.15, 0.68)


# "waterBromine"      
my_Br <- SWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(waterBromine, na.rm=TRUE),
    sd=sd(waterBromine, na.rm=TRUE)
  )

print(ggplot(my_Br) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("Bromine") + theme(axis.text.x = element_text(angle = 90))

