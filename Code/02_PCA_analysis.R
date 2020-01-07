#### PCA analysis #### 
install.packages("factoextra") #install the package if you have not used it before
library(factoextra) #load library to perform PCA 

#read in data 
grab_dat<-read.csv('Data/surface_water_grab.csv', header = TRUE)

variables<-grab_dat[,13:40] #select only columns with variables for PCA, this is columns 13 through 40
variables$siteID<-paste(grab_dat$siteID) #add siteID column back to the table to be able to group later
variables2<-na.omit(variables) # get rid of rows with NA
all.pca<-prcomp(variables2[,-29], scale=TRUE) # perform PCA without the siteID which is column 29 

fviz_eig(all.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

# graph results with the site grouping 
fviz_pca_ind(all.pca, label="none", habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.95)
fviz_eig(all.pca)  ###eigen values associated with each PC (how much does each explain the variation in the data)
eig.val<-get_eigenvalue(all.pca) ## table of the eigenvalues 

#plot the axis and the contribution of each variable 
fviz_pca_var(all.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
