# PCA analysis, Modified 6Aug2020 by JWEdmonds ## 
#modified by KK on 13Sep2020
#modified again by JWE on 20Sep20 when discovered that water chemistry parameters not in same order
install.packages("factoextra") #install the package if you have not used it before
library(factoextra) #load library to perform PCA 
library(devtools)
library(corrplot)
library(ggpubr)
library(scales)
library(viridis)  
library(ggsci)
library(ggplot2)
library(gridExtra)
library(FactoMineR)

#### All 13 sites and all variables for investigation ####
#read in data for file with 10 sites (+3 replacement) chosen, sites = 13 #
#Had to manually remove three sites because could not use grep function correctly when reran code
surface_grab_dat<-read.csv('Data/surface_water_grab_QAQC.csv', header = TRUE) 

variables<-surface_grab_dat[,4:37] #select only columns with variables for PCA
variables$siteID<-paste(surface_grab_dat$siteID) #add siteID column back to the table to be able to group later

#removing NA rows takes rows of data from 748 to 684#
variables2<-na.omit(variables) # get rid of rows with NA

#perform PCA without the siteID which is column 35
Thirteensites_allvariables.pca<-prcomp(variables2[, -35], retx=TRUE,  scale=TRUE) 
# graph results of PCA with the site grouping using 13 sites and ALL variables 
PCA_allvariables<-fviz_pca_ind(Thirteensites_allvariables.pca, label="none", 
                               habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)

PCA_allvariables +ylim(-5,3.5)+xlim(-5, 5)

#graph vectors of PCA, 13 sites, all variables
fviz_pca_var(Thirteensites_allvariables.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#shows the top variables contributing to the 1st and 2nd axis
#red dotted line shows random
fviz_contrib(Thirteensites_allvariables.pca, 
             choice = c("var"), 
             axes = 1:2, #first PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10 

fviz_contrib(Thirteensites_allvariables.pca, 
             choice = c("var"), 
             axes = 2, #2nd PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10

####Reduced number of variables down to 11 ####
#select out Mg, Mn, DOC, TN, Cl, Ca, DIC, NO3NO2, K, Fe, Na, siteID
variables_subset<-select(variables2, Mg, Mn, DOC, TN, Cl, Ca, DIC, 'NO3.NO2...N', K, Fe, Na, siteID) #removing redundant variables#
Thirteensites_subsetvar.pca<-prcomp(variables_subset[, -12], scale=TRUE) #perform PCA with reduced dataset (13 sites), without the siteID which is column 13 

#Scree plot
fviz_eig(Thirteensites_subsetvar.pca)  ## how much does each axis explain the variation in the data
eig.val<-get_eigenvalue(Thirteensites_subsetvar.pca) ## table of the eigenvalues 

# graph results with 13 sites, grouping using reduced variable dataset
PCA_reducedvariables<-fviz_pca_ind(Thirteensites_subsetvar.pca, label="none", 
                                    habillage=variables_subset$siteID,  # group by the site
                                    addEllipses=TRUE, ellipse.level=0.75)

PCA_reducedvariables +ylim(-3,4)+xlim(-3.5, 4.5)

#plot the axis and the contribution of each variable using reduced variable dataset (13 sites)
fviz_pca_var(Thirteensites_subsetvar.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#look at correlation and representation of variables Cos2
var<-get_pca_var(Thirteensites_subsetvar.pca)

corrplot(var$cos2, is.corr=FALSE) 
head(var$cos2, 4)

fviz_cos2(Thirteensites_subsetvar.pca, choice = "var", axes = 1:2)
fviz_pca_var(Thirteensites_subsetvar.pca, alpha.var = "cos2")

set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(Thirteensites_subsetvar.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")


fviz_pca_ind(Thirteensites_subsetvar.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE # Avoid text overlapping (slow if many points)
)


# graph results with the site grouping and have vectors displayed, 13 sites, reduced variable set
Biplot_reducedvariables<-fviz_pca_biplot(Thirteensites_subsetvar.pca, font.legend=c(12,"bold","black"), font.x=c(14,"plain", "black"), font.y=c(14,"plain", "black"), 
             arrowsize = 0.5, labelsize = 5, label="var", col.var="black", axes.linetype = "solid", habillage=variables_subset$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75, 
             font.main=24, subtitle = "Principal Component Analysis", caption = "Source: factoextra",
             xlab = "PC1 (30.6%)", ylab = "PC2 (18.5%)", 
             legend.title = "Sites", legend.position = "top", 
             title = "NEON Surface Water Chemistry", repel=TRUE)

Biplot_reducedvariables + ylim(-3,5.5) + xlim(-3, 6.5) + 
          scale_color_ucscgb() + 
           scale_fill_ucscgb() +
          border(color = "black", size = 0.8, linetype = NULL) + 
          theme( # Hide panel borders and remove grid lines
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
  axis.line = element_line(colour = "black") # Change axis line
)


corrplot(var$contrib, is.corr=FALSE) #plot to show contribution for axis
#contribution to 1st axis
fviz_contrib(Thirteensites_subsetvar.pca, choice = "var", axes = 1, top = 10)
#contribution to 2nd axis
fviz_contrib(Thirteensites_subsetvar.pca, choice = "var", axes = 2, top = 10)
#contribution to both 
fviz_contrib(Thirteensites_subsetvar.pca, choice = "var", axes = 1:2, top = 10)

#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%. 
#For a given component, a variable with a contribution larger than this cutoff could be 
#considered as important in contributing to the component.


Thirteensites_subsetvar.pca$quanti.sup ### ??? dont know what this is 

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- Thirteensites_subsetvar.pca$rotation
sdev <- Thirteensites_subsetvar.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
print(var.coord[, 1:2])

#### reduce data to 10 sites ####
#Created data.frame reduced to 10 sites, using replacement sites (POSE for WALK; BLWA for TOMB; GUIL for CUPE)
tensites<- filter(variables_subset, 
                               !grepl("CUPE", siteID) & 
                                 !grepl("WALK", siteID) & 
                              !grepl("TOMB", siteID))

Tensites.pca<-prcomp(tensites[, -12], scale=TRUE) #perform PCA with reduced dataset (10 sites), without the siteID which is column 13 
#Plot PCA with 10 sites (replacement sites included)
PCA_replacement<-fviz_pca_ind(Tensites.pca, label="none", 
                               habillage=tensites$siteID,  # group by the site
                                    addEllipses=TRUE, ellipse.level=0.75)

PCA_replacement +ylim(-3,3.5)+xlim(-4.5, 3.5)
#Graph results with the site grouping and vectors displayed, 10 sites, reduced variable set
Biplot_replacement<-fviz_pca_biplot(Tensites.pca, label="var", col.var="black", 
                                    habillage=tensites$siteID,  # group by the site
                                         addEllipses=TRUE, ellipse.level=0.75)

#Biplot_replacement +ylim(-3.0,5)+xlim(-5.6, 3.0)

#### Hierarchical clustering #### 
pca_result<-PCA(tensites[, -12], scale.unit = TRUE, graph=FALSE) # then data are scaled to unit variance
tensites.hc<-HCPC(pca_result, graph = FALSE)

fviz_dend(tensites.hc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

fviz_cluster(tensites.hc,
             repel = FALSE,            #  label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

## clust number is assigned to each observation 
clusters<-tensites.hc$data.clust
clusters$siteID<-tensites$siteID #add back in the siteID to see what sites are in each cluster

##look at the variables that explain clustering the most 
tensites.hc$desc.var$quanti
##look at some sites within each cluster, for each cluster, the top 5 closest individuals to the cluster center is shown
tensites.hc$desc.ind$para

###################run PCA for surface water data with all sites #########################

SWgrab_allsites_dat<-read.csv('Data/sw_grab_allsites.csv', header = TRUE)

variables_allsites<-SWgrab_allsites_dat[,4:37] #select only columns with variables for PCA, this is columns 11 through 37 (row 38 has all NA so cancels out all rows if included)
variables_allsites$siteID<-paste(SWgrab_allsites_dat$siteID) #add siteID column back to the table to be able to group later
variables_allsites<-na.omit(variables_allsites) # get rid of rows with NA

summary(as.factor(variables_allsites$siteID))

#remove sites that have 3 pointsor less, not enough for an ellipse
variables_allsites$siteID<-as.factor(variables_allsites$siteID)
variables2<- dplyr::filter(variables_allsites, 
                           siteID == "ARIK" |siteID == "BIGC" |
                             siteID == "BLDE" |siteID == "BLUE"|
                             siteID == "BLWA"|siteID == "CARI"|
                             siteID == "COMO"| siteID == "CUPE"| 
                             siteID == "FLNT" |siteID == "HOPB"| 
                             siteID == "LECO"| siteID == "LEWI" | 
                             siteID == "MAYF"| siteID == "MCDI"|
                             siteID == "MCRA" | siteID == "POSE" |
                             siteID == "PRIN"|siteID == "WLOU") %>% 
  droplevels()
                    
                    
summary(variables2$siteID)             

#run PCA
allsites_allvariables.pca<-prcomp(variables2[, -35], scale=TRUE) #perform PCA without the siteID which is column 35

# graph results with the site grouping, all variables 
PCAallsites_allvar<-fviz_pca_ind(allsites_allvariables.pca, label="none", 
                                 habillage=variables2$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)
#PCAallsites_allvar +ylim(-4,5)+xlim(-7.5, 4.0) #adjust axis if needed

#graph vectors of PCA
fviz_pca_var(allsites_allvariables.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#shows the top variables contributing to the 1st and 2nd axis
#red dotted line shows random
fviz_contrib(allsites_allvariables.pca, 
             choice = c("var"), 
             axes = 1:2, #both PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10 
fviz_contrib(allsites_allvariables.pca, 
             choice = c("var"), 
             axes = 1, #1st PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10
fviz_contrib(allsites_allvariables.pca, 
             choice = c("var"), 
             axes = 2, #2nd PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10

#HCPC cluster analysis 
pca_result<-PCA(variables2[, -35], scale.unit = TRUE, graph=FALSE) # then data are scaled to unit variance
allsites.hc<-HCPC(pca_result, graph = TRUE) #can click on where to cut graph

fviz_dend(allsites.hc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

fviz_cluster(allsites.hc,
             repel = FALSE,            #  label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

## clust number is assigned to each observation 
clusters<-allsites.hc$data.clust
clusters$siteID<-variables2$siteID #add back in the siteID to see what sites are in each cluster

##look at the variables that explain clustering the most 
table<-allsites.hc$desc.var$quanti
##look at some sites within each cluster, for each cluster, the top 5 closest individuals to the cluster center is shown
allsites.hc$desc.ind$para

## observations from previous selection of 10 sites
#LEWI, REDB and PRIN have 8 samples, so they were removed even though they inhabit a very unique PCA space
#TECR has very few samples, and overlaps with COMO cluster
#SYC too few samples, overlaps
#BIGC overlaps with HOPB
#WLOU has very few samples and overlaps with HOPB
#OKSR, MCDI, LECO, FLNT, MCRA, BDLE, CARI all overlap quite a bit


###################read in Groundwater data for file with 13 sites selected (10 + 3 replacement)#########################

groundwater_grab_dat<-read.csv('Data/groundwater_grab_QAQC_3gone.csv', header = TRUE)

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
PCA_biplot_GW_reducedvar<-fviz_pca_biplot(subsetFINAL_subsetvar_GW.pca, font.legend=c(12,"bold","black"), font.x=c(14,"plain", "black"), font.y=c(14,"plain", "black"), 
                                          arrowsize = 0.5, labelsize = 5, label="var", col.var="black", axes.linetype = "solid", habillage=variables_subsetFINAL_GW$siteID,  # group by the site
                                          addEllipses=TRUE, ellipse.level=0.75, font.main=24, subtitle = "Principal Component Analysis", caption = "Source: factoextra",
                                          xlab = "PC1 (29%)", ylab = "PC2 (19.8%)", legend.title = "Sites", legend.position = "top", title = "NEON Groundwater Chemistry", repel=TRUE)


PCA_biplot_GW_reducedvar +xlim(-3, 4.5) +scale_y_reverse(limits=c(1.5,-4)) +scale_color_ucscgb() +scale_fill_ucscgb() +border(color = "black", size = 0.8, linetype = NULL) + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"),
  
)

print(variables_subsetFINAL_GW)
eig.val <- get_eigenvalue(subsetFINAL_subsetvar_GW.pca)
eig.val
var <- get_pca_var(subsetFINAL_subsetvar_GW.pca)
var
corrplot(var$contrib, is.corr=FALSE) 
fviz_contrib(subsetFINAL_subsetvar_GW.pca, choice = "var", axes = 1, top = 10)
#The red dashed line on the graph above indicates the expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%. 
#For a given component, a variable with a contribution larger than this cutoff could be 
#considered as important in contributing to the component.
fviz_contrib(subsetFINAL_subsetvar_GW.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(subsetFINAL_subsetvar_GW.pca, choice = "var", axes = 1:2, top = 10)

subsetFINAL_subsetvar_GW.pca$quanti.sup

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- subsetFINAL_subsetvar_GW.pca$rotation
sdev <- subsetFINAL_subsetvar_GW.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
print(var.coord[, 1:4])

#Was not working for me to remove this way on 6 July 2020 #Created data.frame reduced to 10 sites, using replacement sites (POSE for WALK; BLWA for TOMB; GUIL for CUPE)
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


