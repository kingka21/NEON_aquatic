#### continuous discharge analysis 
### Written by Bobby Hensley
library(neonUtilities)
library(plyr)
library(plotly)

#' Set date range. Water years begin and end Oct 1.
startDate<-"2017-10"
endDate<-"2019-09"

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="HOPB", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="HOPB",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-siteSummary

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="LEWI", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="LEWI",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="POSE", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Issue with POSE 2019 regression 2. Ommit these.
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$regressionID!="POSE.2019.reg2"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="POSE",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="FLNT", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="FLNT",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="CUPE", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="CUPE",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="GUIL", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="GUIL",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="KING", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="KING",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="MCDI", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="MCDI",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="LECO", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="LECO",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="WALK", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="WALK",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="BLWA", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="BLWA",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="MAYF", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="MAYF",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="ARIK", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="ARIK",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="BLUE", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="BLUE",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="PRIN", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="PRIN",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="BLDE", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="BLDE",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="COMO", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="COMO",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="WLOU", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="WLOU",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="SYCA", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="SYCA",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="REDB", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="REDB",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="MART", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="MART",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="MCRA", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="MCRA",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="BIGC", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="BIGC",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)

#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="TECR", startdate=startDate, enddate=endDate,package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
#' Remove any quality flagged measurements
csd_continuousDischarge<-csd_continuousDischarge[(csd_continuousDischarge$dischargeFinalQF=="0"),]
#' Plot distribution for quick visual check
csd_continuousDischarge$logQ<-log10(csd_continuousDischarge$maxpostDischarge)
fig<-plot_ly(data=csd_continuousDischarge, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="TECR",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)













#' Download continuous discharge data from NEON API
sampleQ <-neonUtilities::loadByProduct(dpID="DP1.20048.001", site="OKSR", startdate="2016-10", enddate="2021-06",package="basic", check.size = F)
list2env(sampleQ, .GlobalEnv)
sampleQ<-dsc_fieldData
plot_ly(data=sampleQ, x=~collectDate,y=~totalDischarge,name="Q",type="scatter",mode="lines",line=list(color="black",width=2))%>%
  layout(title ="OKSR",titlefont=list(size=30), xaxis = list(title = "Date",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "Q (L/s)",type="log",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
#' Calculates 15 min averages
sampleQ$collectDate<-lubridate::round_date(sampleQ$collectDate,unit="1 month")
sampleQ<-plyr::ddply(sampleQ,c("collectDate"),summarise,meanQ=mean(totalDischarge))
allDates<-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="COMO", startdate="2016-10", enddate="2021-05",package="basic", check.size = F)
list2env(allDates, .GlobalEnv)
allDates<-csd_continuousDischarge[,c("endDate","siteID")]
allDates$endDate<-lubridate::round_date(allDates$endDate,unit="1 month")
allDates<-plyr::ddply(allDates,c("endDate"),summarise,stationHorizontalID=mean(stationHorizontalID))
sampleQ<-merge(allDates,sampleQ, by.x="endDate",by.y="collectDate",all.x=T, all.y=T)
for(i in 1:nrow(sampleQ)){if(is.na(sampleQ[i,3])){sampleQ[i,3]=0}} 
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(sampleQ,c("siteID"),summarise,meanQ=mean(totalDischarge,na.rm=TRUE),sdQ=sd(totalDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)





#' Download continuous discharge data from NEON API
contQ <-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="CARI", startdate="2018-10", enddate="2019-09",package="basic", check.size = F)
list2env(contQ, .GlobalEnv)
contQ<-csd_continuousDischarge[,c("endDate","maxpostDischarge")]
allDates<-neonUtilities::loadByProduct(dpID="DP4.00130.001", site="COMO", startdate="2017-10", enddate="2019-09",package="basic", check.size = F)
list2env(allDates, .GlobalEnv)
allDates<-csd_continuousDischarge[,c("endDate","siteID")]
contQ<-merge(allDates,contQ, by.x="endDate",by.y="endDate",all.x=T, all.y=F)
for(i in 1:nrow(contQ)){if(is.na(contQ[i,3])){contQ[i,3]=0}} 
plot_ly(data=contQ, x=~endDate,y=~maxpostDischarge,name="Q",type="scatter",mode="lines",line=list(color="black",width=2))%>%
  layout(title ="CARI",titlefont=list(size=30), xaxis = list(title = "Date",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "Q (L/s)",type="log",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(csd_continuousDischarge,c("siteID"),summarise,meanQ=mean(maxpostDischarge,na.rm=TRUE),sdQ=sd(maxpostDischarge,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)





#' TOMB pulls data from USGS API
siteNumber <- "02469761"
parameterCd <- "00060"  # Discharge in cfs
startDate <- "2017-10-01" 
endDate <- "2019-10-01" 
usgsData<-dataRetrieval::readNWISuv(siteNumber, parameterCd, startDate, endDate)
#' Convertes discharge from cfs to L/s
usgsData$dischargeLps=usgsData$X_00060_00000*0.0283*1000
usgsData$siteID="TOMB"
#' Plot distribution for quick visual check
usgsData$logQ<-log10(usgsData$dischargeLps)
fig<-plot_ly(data=usgsData, x=~logQ,type="histogram",histnorm="probability",marker=list(color="orange"))
fig<-fig%>% layout(title ="TOMB",titlefont=list(size=30), xaxis = list(title = "LogQ (L/s)",titlefont=list(size=30),tickfont=list(size=30)),yaxis = list (title = "p(x)",titlefont=list(size=30),tickfont=list(size=30)),margin=list(l=50, r=50, t=100, b=100, pad=4))
fig
#' Calculate mean and standard deviation
siteSummary<-plyr::ddply(usgsData,c("siteID"),summarise,meanQ=mean(dischargeLps,na.rm=TRUE),sdQ=sd(dischargeLps,na.rm=TRUE))
#' Append to table of all sites
allSites<-rbind(allSites,siteSummary)




allSites<-na.omit(allSites)




