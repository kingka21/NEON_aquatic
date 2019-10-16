#### Script for downloading and organizing NEON datasets ##### 
### code written by Katelyn King 16-OCT-2019 
### code adapted from Bobby Hensley ### 

### load libraries ###
# install neonUtilities
install.packages("neonUtilities")
# load neonUtilities
library(neonUtilities)
?loadByProduct
?getDatatable


### load datasets ###
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

#Water quality Sensor data
#In situ sensor-based specific conductivity, concentration of chlorophyll a, dissolved oxygen content, fDOM concentration, pH, and turbidity, available as one-, five-, and thirty-minute averages in surface water of lakes, wadeable streams, and non-wadeable streams.
waterchem_sensor <- loadByProduct(dpID="DP1.20288.001", 
                           site=c("MOAB","ONAQ"),
                           startdate="2018-05",  #year and month
                           enddate="2018-08")

#Grab samples of surface water chemistry including general chemistry, anions, cations, and nutrients.
nutrients_grab <- loadByProduct(dpID="DP1.20093.001", 
                           site=c("MOAB","ONAQ"),
                           startdate="2018-05", 
                           enddate="2018-08", 
                           package="expanded", #basic will just give you concentrations, #expanded will give you flags 
                           check.size = T)  ### check the size of the file before you download it 


# Creates data frames from tables 
for(i in 1:length(wq)) {assign(names(wq)[i], wq[[i]])}
View(waq_instantaneous)
wqvalues<-data.frame(waq_instantaneous)

#### Seperates measurements by location 
##  Upstream = 101              Downstream = 102
##  US overhanging = 111        DS overhanging = 112
##  Lake/river buoy = 103
##  Lake inlet = 130            Lake outlet = 140
##  NOTE: SUNA and fDOM only at downstream stream station or lake/river buoy
##  NOTE: function removes values, so enter term for locations you DON'T want
wqvalues101<-wqvalues[(wqvalues$horizontalPosition=="101"),]
wqvalues102<-wqvalues[(wqvalues$horizontalPosition=="102"),]

### Extracts datetime and converts datetime string to POSIXct
wqvalues101$startDateTime<-as.POSIXct(wqvalues101$startDateTime,format="%Y-%m-%dT%H:%M:%OS")
wqvalues102$startDateTime<-as.POSIXct(wqvalues102$startDateTime,format="%Y-%m-%dT%H:%M:%OS")

### Plots data ###
plot(wqvalues101$startDateTime,wqvalues101$dissolvedOxygen,type="l",col="blue",main="ARIK DO",xlab="Date",ylab="DO (mg/L)")
lines(wqvalues102$startDateTime,wqvalues102$dissolvedOxygen,col="red")
grid(nx=NULL,ny=NULL, col="lightgray",lty="dotted",lwd=par("lwd"), equilogs=TRUE)
legend("bottomleft",legend=c("Upstream", "Downstream"),
       col=c("blue", "red"), lty=1:1, cex=0.8)
