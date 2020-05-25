#### Graphing 13 sites of interest for exploring variation across sites (spatial) #### 
#graph bar plot with error bars 
GWgrab_dat<-read.csv('Data/groundwater_grab_QAQC.csv', header = TRUE)


#####pH 
my_pH <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise(
    n=n(),
    mean=mean(pH, na.rm=TRUE),
    sd=sd(pH, na.rm=TRUE)
  )
my_pH$siteID<-ordered(my_pH$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_pH) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW pH") + theme(axis.text.x = element_text(angle = 90)) + ylim (0, 8.4)

#####conductivity
my_con <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(conductivity, na.rm=TRUE),
    sd=sd(conductivity, na.rm=TRUE)
  )
my_con$siteID<-ordered(my_con$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


ggplot(my_con) +
  geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
  geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("GW conductivity")



######NO2.NO3...N      
my_NitrateNitrite <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(NO3.NO2...N, na.rm=TRUE),
    sd=sd(NO3.NO2...N, na.rm=TRUE)
  )
my_NitrateNitrite$siteID<-ordered(my_NitrateNitrite$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_NitrateNitrite) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Nitrate+Nitrite") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.03, 0.941)


#####"TN"  
my_TN <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise(
    n=n(),
    mean=mean(TN, na.rm=TRUE),
    sd=sd(TN, na.rm=TRUE)
  )
my_TN$siteID<-ordered(my_TN$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))

print(ggplot(my_TN) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW TN") + theme(axis.text.x = element_text(angle = 90))  + ylim(-0.05, 1)


#####dissolvedInorganicCarbon     
my_DIC <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(DIC, na.rm=TRUE),
    sd=sd(DIC, na.rm=TRUE)
  )
my_DIC$siteID<-ordered(my_DIC$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_DIC) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW DIC") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.7, 105)


#####Calcium      
my_Ca <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Ca, na.rm=TRUE),
    sd=sd(Ca, na.rm=TRUE)
  )
my_Ca$siteID<-ordered(my_Ca$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Ca) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Calcium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-2.97, 83)


#####Magnesium     
my_Mg <-GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Mg, na.rm=TRUE),
    sd=sd(Mg, na.rm=TRUE)
  )
my_Mg$siteID<-ordered(my_Mg$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Mg) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Magnesium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-1.7, 41)


#####Potassium      
my_K <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(K, na.rm=TRUE),
    sd=sd(K, na.rm=TRUE)
  )
my_K$siteID<-ordered(my_K$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_K) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Potassium") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.05, 2.03)


#####Na      
my_Na <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Na, na.rm=TRUE),
    sd=sd(Na, na.rm=TRUE)
  )
my_Na$siteID<-ordered(my_Na$siteID, levels=c("MAYF", 'COMO', 'MART', "HOPB", 'TOMB','BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Na) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Sodium") + theme(axis.text.x = element_text(angle = 90))  + ylim(-0.05, 27)


#####Manganese  
my_Mn <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Mn, na.rm=TRUE),
    sd=sd(Mn, na.rm=TRUE)
  )
my_Mn$siteID<-ordered(my_Mn$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Mn) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Manganese") + theme(axis.text.x = element_text(angle = 90))



#####dissolvedOrganicCarbon 
my_disorgcarbon <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(DOC, na.rm=TRUE),
    sd=sd(DOC, na.rm=TRUE)
  )
my_disorgcarbon$siteID<-ordered(my_disorgcarbon$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_disorgcarbon) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW DOC") + theme(axis.text.x = element_text(angle = 90))


#####Chlorine
my_Cl <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Cl, na.rm=TRUE),
    sd=sd(Cl, na.rm=TRUE)
  )
my_Cl$siteID<-ordered(my_Cl$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_Cl) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW Chlorine") + theme(axis.text.x = element_text(angle = 90)) + ylim(-0.05, 32)


#####Iron"\  
my_iron <- GWgrab_dat %>%
  group_by(siteID) %>%
  summarise( 
    n=n(),
    mean=mean(Fe, na.rm=TRUE),
    sd=sd(Fe, na.rm=TRUE)
  )
my_iron$siteID<-ordered(my_iron$siteID, levels=c("MAYF", 'COMO' , 'MART', "HOPB" , 'TOMB' ,'BLWA', 'POSE', 'WALK' , 'CUPE', 'GUIL', 'ARIK', 'KING', 'BLUE'))


print(ggplot(my_iron) +
        geom_bar( aes(x=siteID, y=mean), stat="identity", fill="skyblue") +
        geom_errorbar( aes(x=siteID, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)) + 
  ylab("GW iron") + theme(axis.text.x = element_text(angle = 90))  + ylim(-0.05, 1.5)
