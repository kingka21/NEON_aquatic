#### PCA analysis, Modified 21May2020 by JWEdmonds #### 
install.packages("factoextra") #install the package if you have not used it before
library(factoextra) #load library to perform PCA 
library("devtools")
library("factoextra")

#################read in data for file with 10 sites (+3 replacement) chosen, sites = 13 ###########################
surface_grab_dat<-read.csv('Data/surface_water_grab_QAQC.csv', header = TRUE)

variables<-surface_grab_dat[,4:34] #select only columns with variables for PCA, this is columns 3 through 34
variables$siteID<-paste(surface_grab_dat$siteID) #add siteID column back to the table to be able to group later

#removing NA rows takes rows of data from 748 to 684#
variables2<-na.omit(variables) # get rid of rows with NA
write.csv(variables2, 'Data/SW_grab_NoNAs.csv', row.names = FALSE) #checking order of variables so as to determine which ones to keep in next step#

Thirteensites_allvariables.pca<-prcomp(variables2[, -32], retx=TRUE,  scale=TRUE) #perform PCA without the siteID which is column 32
# graph results of PCA with the site grouping using 13 sites and ALL variables 
PCA_allvariables<-fviz_pca_ind(Thirteensites_allvariables.pca, label="none", habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)

PCA_allvariables +ylim(-5,3.5)+xlim(-5, 5)

#graph vectors of PCA, 13 sites, all variables
fviz_pca_var(Thirteensites_allvariables.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#Reduced number of variables down to 12
write.csv(variables2, 'Data/variable_list.csv', row.names = FALSE) #checking order of variables so as to determine which ones to keep in next step#
variables_subset<-subset(variables2, select = -c(1,2,3,6,9,11,12,13,15,19,20,22,23,24,25,27,28,29,30,31)) #removing redundant variables#
write.csv(variables_subset, 'Data/variable_subset.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#
Thirteensites_subsetvar.pca<-prcomp(variables_subset[, -12], scale=TRUE) #perform PCA with reduced dataset (13 sites), without the siteID which is column 13 

fviz_eig(Thirteensites_subsetvar.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(Thirteensites_subsetvar.pca) ## table of the eigenvalues 

# graph results with 13 sites, grouping using reduced variable dataset
PCA_reducedvariables<-(fviz_pca_ind(Thirteensites_subsetvar.pca, label="none", habillage=variables_subset$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75))

PCA_reducedvariables +ylim(-3,3.5)+xlim(-4.5, 3.5)

fviz_eig(Thirteensites_subsetvar.pca)  ###eigen values associated with each PC (how much does each explain the variation in the data)
eig.val<-get_eigenvalue(Thirteensites_subsetvar.pca) ## table of the eigenvalues 

#plot the axis and the contribution of each variable using reduced variable dataset (13 sites)
fviz_pca_var(Thirteensites_subsetvar.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# graph results with the site grouping and have vectors displayed, 13 sites, reduced variable set
Biplot_reducedvariables<-fviz_pca_biplot(Thirteensites_subsetvar.pca, label="var", col.var="black", habillage=variables_subset$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)

Biplot_reducedvariables +ylim(-3.5,5)+xlim(-6.5, 3.0)

#Created data.frame reduced to 10 sites, using replacement sites (POSE for WALK; BLWA for TOMB; GUIL for CUPE)
variables_replacement<- filter(variables_subset, !grepl("CUPE", variables_subset$siteID) & 
                                !grepl("WALK", variables_subset$siteID) & 
                                !grepl("TOMB", variables_subset$siteID))
write.csv(variables_replacement, 'Data/variable_list_replacement.csv', row.names = FALSE) #checking order of variables and 574 observations#

Tensites_replacement.pca<-prcomp(variables_replacement[, -12], scale=TRUE) #perform PCA with reduced dataset (10 sites), without the siteID which is column 13 
#Plot PCA with 10 sites (replacement sites included)
PCA_replacement<-(fviz_pca_ind(Tensites_replacement.pca, label="none", habillage=variables_replacement$siteID,  # group by the site
                                    addEllipses=TRUE, ellipse.level=0.75))

PCA_replacement +ylim(-3,3.5)+xlim(-4.5, 3.5)
#Graph results with the site grouping and vectors displayed, 10 sites, reduced variable set
Biplot_replacement<-fviz_pca_biplot(Tensites_replacement.pca, label="var", col.var="black", habillage=variables_replacement$siteID,  # group by the site
                                         addEllipses=TRUE, ellipse.level=0.75)

Biplot_replacement +ylim(-3.0,5)+xlim(-5.6, 3.0)


###################read in surface water data for file with all sites #########################
SWgrab_allsites_dat<-read.csv('Data/surface_water_grab_allsites_PIVOT.csv', header = TRUE)

variables_allsites<-SWgrab_allsites_dat[,4:34] #select only columns with variables for PCA, this is columns 11 through 39
variables_allsites$siteID<-paste(SWgrab_allsites_dat$siteID) #add siteID column back to the table to be able to group later
variables2_allsites<-na.omit(variables_allsites) # get rid of rows with NA
write.csv(variables2_allsites, 'Data/variable_allsites_SW.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#

allsites_allvariables.pca<-prcomp(variables2_allsites[, -32], scale=TRUE) #perform PCA without the siteID which is column 30
# graph results with the site grouping, all variables 
PCAallsites_allvar<-fviz_pca_ind(allsites_allvariables.pca, label="none", habillage=variables2_allsites$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)
PCAallsites_allvar +ylim(-4,5)+xlim(-7.5, 4.0)

variables3_allsites<-subset(variables2_allsites, select = -c(1,2,4,5,6,7,10,13,14,15,16,17,19,20,21,22,23,26,29,31)) #removing redundant variables#
allsites_subsetvar.pca<-prcomp(variables3_allsites[, -12],  scale. = T) #perform PCA with reduced dataset (10 sites), without the siteID which is column 12 
write.csv(variables3_allsites, 'Data/variables_allsites.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#


PCAallsites_subsetvar<-fviz_pca_ind(allsites_subsetvar.pca, label="none", habillage=variables3_allsites$siteID,  # group by the site
             addEllipses=FALSE, ellipse.level=0.75)
PCAallsites_subsetvar +ylim(-2.5,2.5)+xlim(-3, 3.5)

#LEWI, REDB and PRIN have 8 samples, so they were removed even though they inhabit a very unique PCA space
#TECR has very few samples, and overlaps with COMO cluster
#SYC too few samples, overlaps
#BIGC overlaps with HOPB
#WLOU has very few samples and overlaps with HOPB
#OKSR, MCDI, LECO, FLNT, MCRA, BDLE, CARI all overlap quite a bit
variables4_moresites<- filter(variables3_allsites, !grepl("MCDI", variables3_allsites$siteID) & 
                                                    !grepl("OKSR", variables3_allsites$siteID) & 
                                                    !grepl("TECR", variables3_allsites$siteID) & 
                                                    !grepl("LECO", variables3_allsites$siteID) & 
                                                    !grepl("CARI", variables3_allsites$siteID) &
                                                    !grepl("LEWI", variables3_allsites$siteID) &
                                                    !grepl("PRIN", variables3_allsites$siteID) &
                                                    !grepl("REDB", variables3_allsites$siteID) &
                                                    !grepl("WLOU", variables3_allsites$siteID) &
                                                    !grepl("SYC", variables3_allsites$siteID) &
                                                    !grepl("MCRA", variables3_allsites$siteID) &
                                                    !grepl("BLDE", variables3_allsites$siteID) &
                                                    !grepl("BIGC", variables3_allsites$siteID) &
                                                    !grepl("FLNT", variables3_allsites$siteID))
write.csv(variables4_moresites, 'Data/variables4_somesites.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#


testsites_subsetvar.pca<-prcomp(variables4_moresites[, -12],  scale. = T) #perform PCA with reduced dataset (10 sites), without the siteID which is column 13 
PCAtestsites_subsetvar<-fviz_pca_ind(testsites_subsetvar.pca, label="none", habillage=variables4_moresites$siteID,  # group by the site
                                    addEllipses=FALSE, ellipse.level=0.75)
PCAtestsites_subsetvar +ylim(-3,3)+xlim(-2.5, 5)

# graph results with the site grouping and have vectors displayed, reduced variable set
fviz_pca_biplot(allsites_subsetvar.pca, label="var", col.var="black", habillage=variables_subset$siteID,  # group by the site
                addEllipses=TRUE, ellipse.level=0.66)


###################read in Groundwater data for file with 13 sites selected (10 + 3 replacement)#########################

groundwater_grab_dat<-read.csv('Data/groundwater_grab_QAQC.csv', header = TRUE)

variables_GW<-groundwater_grab_dat[,4:37] #select only columns with variables for PCA, this is columns 11 through 39
write.csv(variables_GW, 'Data/variable_GW.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#

variables_GW$siteID<-paste(groundwater_grab_dat$siteID) #add siteID column back to the table to be able to group later
variables2_GW<-na.omit(variables_GW) # get rid of rows with NA
GWsubsetFINAL_allvariables_GW.pca<-prcomp(variables2_GW[, -35], scale=TRUE) #perform PCA without the siteID which is column 30
# graph results with the site grouping using nine sites (with CUPE/GUIL) and all variables 
PCA_allvariables_GW<-fviz_pca_ind(GWsubsetFINAL_allvariables_GW.pca, label="none", habillage=variables2_GW$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)

PCA_allvariables_GW +ylim(-3.5,3)+xlim(-4.5, 4.5)

variables_subsetFINAL_GW<-subset(variables2_GW, select = -c(1,2,3,4,5,7,10,11,14,15,16,18,19,21,22,23,24,25,26,27,31,33,34)) #removing redundant variables#
write.csv(variables_subsetFINAL_GW, 'Data/variable_subset_GW.csv', row.names = FALSE) #checking variables so as to determine you have the correct in next step#
subsetFINAL_subsetvar_GW.pca<-prcomp(variables_subsetFINAL_GW[, -12], scale=TRUE) #perform PCA with reduced dataset (9 sites), without the siteID which is column 12 

fviz_eig(subsetFINAL_subsetvar_GW.pca)  ###Scree plot: how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(subsetFINAL_subsetvar_GW.pca) ## table of the eigenvalues 

# graph GW results with the site grouping using reduced variable dataset (11 sites)
PCA_plot_GW_reducedvar<-fviz_pca_ind(subsetFINAL_subsetvar_GW.pca, label="none", habillage=variables_subsetFINAL_GW$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)
PCA_plot_GW_reducedvar  +ylim(-3,2)+xlim(-3,3.8)

fviz_eig(subsetFINAL_subsetvar_GW.pca)  ###eigen values associated with each PC (how much does each explain the variation in the data)
eig.val<-get_eigenvalue(subsetFINAl_subsetvar_GW.pca) ## table of the eigenvalues 

#plot the axis and the contribution of each variable using reduced variable dataset for groundwater (9 sites with CUPE/GUIL)
fviz_pca_var(subsetFINAL_subsetvar_GW.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# graph groundwater results with the site grouping and have vectors displayed, final 10 sites, reduced variable set
PCA_biplot_GW_reducedvar<-fviz_pca_biplot(subsetFINAL_subsetvar_GW.pca, label="var", col.var="black", habillage=variables_subsetFINAL_GW$siteID,  # group by the site
                addEllipses=TRUE, ellipse.level=0.75)

PCA_biplot_GW_reducedvar +ylim(-3,2)+xlim(-3,3.8)

#Created data.frame reduced to 10 sites, using replacement sites (POSE for WALK; BLWA for TOMB; GUIL for CUPE)
GWvariables_replacement<- filter(variables_subsetFINAL_GW, !grepl("CUPE", variables_subsetFINAL_GW$siteID) & 
                                 !grepl("WALK", variables_subsetFINAL_GW$siteID) & 
                                 !grepl("TOMB", variables_subsetFINAL_GW$siteID))
write.csv(GWvariables_replacement, 'Data/GWvariable_replacement.csv', row.names = FALSE) #checking order of variables and 574 observations#

GWTensites_replacement.pca<-prcomp(GWvariables_replacement[, -12], scale=TRUE) #perform PCA with reduced dataset (10 sites), without the siteID which is column 12 
#Plot PCA with 10 sites (replacement sites included)
PCA_GWreplacement<-(fviz_pca_ind(GWTensites_replacement.pca, label="none", habillage=GWvariables_replacement$siteID,  # group by the site
                               addEllipses=TRUE, ellipse.level=0.75))

PCA_GWreplacement +ylim(-2.5,1.65)+xlim(-3.5, 3.6)

#Graph results with the site grouping and vectors displayed, 10 sites, reduced variable set
Biplot_GWreplacement<-fviz_pca_biplot(GWTensites_replacement.pca, label="var", col.var="black", habillage=GWvariables_replacement$siteID,  # group by the site
                                    addEllipses=TRUE, ellipse.level=0.75)

Biplot_GWreplacement  +ylim(-4.0,1.5)+xlim(-3.0, 4)


