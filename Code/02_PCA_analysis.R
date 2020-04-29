#### PCA analysis #### 
install.packages("factoextra") #install the package if you have not used it before
library(factoextra) #load library to perform PCA 
library(dplyr)
#read in data 
grab_dat<-read.csv('Data/surface_water_grab.csv', header = TRUE)

variables<-sw_data[,6:48] #select only columns with variables for PCA, this is columns 6 through 48
var2<-subset(variables, select = -c(DATE, decimalLatitude, decimalLongitude, totalDischarge))
var2$siteID<-paste(sw_data$siteID) #add siteID column back to the table to be able to group later
variables2<-na.omit(var2) # get rid of rows with NA
all.pca<-prcomp(variables2[,-40]) # perform PCA without the siteID which is column 29 

fviz_eig(all.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

# graph results with the site grouping 
fviz_pca_ind(all.pca, label="none", habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.95)

#plot the axis and the contribution of each variable 
fviz_pca_var(all.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


### NEW PLOTS 
variables2<-na.omit(variables) # get rid of rows with NA
data_forPCA<-variables2[,3:14] #select only columns with variables for PCA, this is columns 6 through 48

all.pca<-prcomp(data_forPCA, scale=TRUE) # perform PCA

fviz_eig(all.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

# graph results with the site grouping 
fviz_pca_ind(all.pca, label="var", habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.95)


#plot the axis and the contribution of each variable 
fviz_pca_var(all.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#plot ellises and the lines 
fviz_pca_biplot(all.pca, label="var", habillage=variables2$siteID,  # group by the site
                addEllipses=TRUE, ellipse.level=0.95) 


### MAP 
#decimalLatitude and decimalLongitude 
latlon<-sw_data[!with(sw_data,is.na(decimalLongitude)),] %>% 
        distinct(siteID, .keep_all = TRUE) %>%
  select(siteID, decimalLatitude, decimalLongitude ) %>%
  filter(siteID == "ARIK" | siteID == "COMO" | siteID == "CUPE" | siteID == "GUIL" | siteID == "KING" | siteID == "MAYF" | siteID == "TOMB" | siteID == "BLUE" | siteID == "FLNT" ) 

latlon<-sw_data[!with(sw_data,is.na(decimalLongitude)),] %>% 
  distinct(siteID, .keep_all = TRUE) %>%
  select(siteID, decimalLatitude, decimalLongitude ) 


usa<-map_data("usa")  #pull out the usa map
base_map<-ggplot(data = usa) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 

base_map + geom_point(data = latlon, aes(x = decimalLongitude, y = decimalLatitude)) + 
                geom_label(data = latlon, aes(x = decimalLongitude, y = decimalLatitude, label=siteID),
                           nudge_x = 0,
                           nudge_y = 0 )

#HOPB, POSE, WALK, MAYF, KING, PRIN, COMO, BLDE, REDB, TECR, MART
grab_dat<-read.csv('Data/surface_water_grab.csv', header = TRUE)
var<-select(grab_dat, siteID, 
                  externalConductance, dissolvedOrganicCarbon, waterIron, 
                  waterNitrateAndNitriteN, waterManganese, waterChlorine, 
                  waterFluorine, waterPotassium, waterSodium, 
                  waterCalcium, dissolvedInorganicCarbon, waterMagnesium) %>%
  filter(siteID == "HOPB" | siteID == "POSE" | siteID == "WALK" | siteID == "MAYF" | siteID == "KING" | siteID == "PRIN" | siteID == "COMO" | siteID == "BLDE" | siteID == "REDB" | siteID == "TECR" | siteID == "MART")  

#variables2<-na.omit(var) # get rid of rows with NA
var<-var[!with(var,is.na(dissolvedOrganicCarbon)),]
var<-var[!with(var,is.na(externalConductance)),]

data_forPCA<-var[,2:13] #select only columns with variables for PCA, this is columns 6 through 48

all.pca<-prcomp(data_forPCA) # perform PCA

fviz_eig(all.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

# graph results with the site grouping 
fviz_pca_ind(all.pca, label="none", habillage=var$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.95)


#plot the axis and the contribution of each variable 
fviz_pca_var(all.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)            

fviz_pca_biplot(all.pca, label="var", habillage=var$siteID,  # group by the site
                addEllipses=TRUE, ellipse.level=0.95)


#### PCA with CHL and only 10 sites 
### variables is from merging sensor and grab data 
#select only rows with all of the data points 
variables2<-na.omit(variables) 

data_forPCA<-variables2[,3:14] #select only columns with variables for PCA, this is columns 6 through 48

all.pca<-prcomp(data_forPCA) # perform PCA

fviz_eig(all.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

# graph results with the site grouping 
fviz_pca_ind(all.pca, label="none", habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.95)


#plot the axis and the contribution of each variable 
fviz_pca_var(all.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)            

fviz_pca_biplot(all.pca, label="var", habillage=variables2$siteID,  # group by the site
                addEllipses=TRUE, ellipse.level=0.95)

#### Random Forest 
install.packages('randomForest')
library(randomForest)
library(ggplot2)
RF_chl <- randomForest(chlorophyll_mean ~ dissolvedOrganicCarbon + waterIron + waterNitrateAndNitriteN +
                         waterManganese + waterChlorine + waterFluorine + waterPotassium + waterSodium +
                         waterCalcium + dissolvedInorganicCarbon + waterMagnesium,
                       data=variables2, ntree=10, importance=T, na.action=na.omit)

TP_imp <-randomForest::importance(RF_chl, type=1, scale=FALSE) #mean decrease in accuracy (also called permutation accuracy importance).
imp<-as.data.frame(TP_imp) 
imp$'Pred'   <-rownames(imp)
imp <- structure(imp$`%IncMSE`, names = as.character(imp$Pred))

dotchart(imp[order(imp)], xlab = "mean decrease in accuracy",
         main = set_string[i])
