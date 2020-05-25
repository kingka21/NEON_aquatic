
#Combine groundwater and surface water data by matching the dates, so as to do further correlations#

install.packages(fuzzyjoin)
library(dplyr)
library(fuzzyjoin)

#combine two datasets using exact matching of dates#
SWgrab_chem_dat_subsetFINAL_PIVOT$DATE<-as.Date(SWgrab_chem_dat_subsetFINAL_PIVOT$collectDate, format="%Y-%m")
GWgrab_chem_dat_subsetFINAL_PIVOT$DATE<-as.Date(GWgrab_chem_dat_subsetFINAL_PIVOT$collectDate, format="%Y-%m")
all_data<-left_join(SWgrab_chem_dat_subsetFINAL_PIVOT, GWgrab_chem_dat_subsetFINAL_PIVOT, by = c("siteID", "DATE"))
write.csv(all_data, 'Data/SWgrabPLUSGWgrab13sites_exactmatch.csv', row.names = FALSE)

###### Combine Files with Imperfect Match on Date###############################################################
###### written by Aaron Wong, associate professor at Nevada State college#######################################

surface <- SWgrab_chem_dat_subsetFINAL_PIVOT
groundwater <- GWgrab_chem_dat_subsetFINAL_PIVOT

# Create an empty data frame where the desired information will be placed
combined <- data.frame()

# Convert dates to numeric for time processing
surface$collectDateNumeric = as.numeric(as.Date(surface$collectDate, format="%Y-%m-%d"))
groundwater$collectDateNumeric = as.numeric(as.Date(groundwater$collectDate, format="%Y-%m-%d"))

# Loop through siteIDs
for ( siteID in unique(surface$siteID) )
{
  # Get GW data according to siteID
  surface_site = surface[ surface$siteID == siteID, ]
  groundwater_site = groundwater[ groundwater$siteID == siteID, ]
  
  # Loop through each surface collection from this site
  for ( i in 1:nrow(surface_site ) )
  {
    # Calculate the difference in times from collection to discharge
    groundwater_site$timediff = groundwater_site$collectDateNumeric - surface_site$collectDateNumeric[i]
    
    # Calculate the absolute time difference
    groundwater_site$abstimediff = abs(groundwater_site$collectDateNumeric - surface_site$collectDateNumeric[i])
    
    # Sort the rows in discharge_site by the absolute time difference
    sorted_groundwater_site <- groundwater_site[order( groundwater_site$abstimediff), ]
    
    # Get the row with the smallest absolute time difference
    smallest_diff = sorted_groundwater_site[1,]
    
    # Create a combined row of data by merging the surface_site row with the smallest_diff row
    new_row <- left_join( surface_site[i,], smallest_diff,  by = c("siteID") )
    
    # Add the new row to the data frame
    combined <- rbind( combined, new_row )
  }
}

write.csv(combined, 'Data/CombinedSW_GW.csv', row.names = FALSE)




#####################################################################################################################

install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
my_data <- all_data[,13:24] ### select out what columns you want to compare (I am chasing columns 3 through 14 here) OR you can do it like this and choose the individual column numbers if they aren't in order:  
my_data <- data[, c(3,10,12,13,14)]

#then use this function to plot. the histogram = TRUE part adds a histogram of the variable so you can see the normality. 
chart.Correlation(my_data, histogram=TRUE, pch=19)