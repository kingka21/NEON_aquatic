# PCA analysis 
#Created by KK  March 30, 2020
#Modified 6Aug2020 by JWEdmonds, 20Sep20  ## 
#modified by KK on 13Sep2020; 28Sep 2020, Oct 5, 2020 


install.packages("factoextra") #install the package if you have not used it before
library(factoextra) #library for great visualization of PCA results
library(FactoMineR) #library to perfomr PCA and HCPC clusters
library(ggplot2)
library(dplyr) #library for data maniupulation 
library(corrplot)
library(ggpubr)
library(scales)
library(viridis)  
library(ggsci)
library(gridExtra)
install.packages("psych")
library(psych)
###################run PCA for surface water data with all sites #########################

#bring in most recent data, QAQC and all sites updated Sep 28 
sw_all<-read.csv('Data/surface_water_grab_QCQC_ph.csv', header = TRUE)
SWgrab_allsites_dat<-select(sw_all, -c(sampleID, collectDate, pH, TDS, Br)) %>% #remove columns I don't need
                     data.table:: setnames(old=c("NH4...N", "NO3.NO2...N", "initialSamplepH"), new=c("NH4N", "NO3NO2N", "pH"))

summary(SWgrab_allsites_dat$siteID)

#keep all sites - they all have >20 obs! 
#sites_w_data<- dplyr::filter(SWgrab_allsites_dat, 
 #                          siteID == "ARIK" |siteID == "BIGC" |
  #                           siteID == "BLDE" |siteID == "BLUE"|
   #                          siteID == "BLWA"|siteID == "CARI"|
    #                         siteID == "COMO"| 
     #                        siteID == "FLNT" | siteID == "LECO"|
      #                      siteID == "LEWI" | siteID == "MART" | 
       #                      siteID == "MAYF"| siteID == "MCDI"|
        #                     siteID == "REDB"|
         #                    siteID == "TOMB") %>% 
      #                        droplevels()
                    

#*run PCA ####
SWgrab_allsites_dat<-na.omit(SWgrab_allsites_dat) # rows where lab pH was not available (from 1504 to 1385)
allsites.pca<-prcomp(SWgrab_allsites_dat[, -1], scale=TRUE) #perform PCA without the siteID 

#rotation 
rawLoadings     <- allsites.pca$rotation[,1:3] %*% diag(allsites.pca$sdev, 3, 3)
rotatedLoadings <- varimax(rawLoadings)$loadings
scores <- scale(allsites.pca$x[,1:3]) %*% varimax(rawLoadings)$rotmat
print(scores[1:5,])  

#eigenvalues
fviz_eig(allsites.pca)  ###eigen values associated with each PC (how much does each explain the variation in the data)
eig.val<-get_eigenvalue(allsites.pca) ## table of the eigenvalues 

# graph results with the site grouping, all variables 
fviz_pca_ind(allsites.pca, label="none", 
                                 habillage=SWgrab_allsites_dat$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)
#PCAallsites_allvar +ylim(-4,5)+xlim(-7.5, 4.0) #adjust axis if needed

#graph vectors of PCA
vectors_SW<-fviz_pca_var(allsites.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             title = "Surface water")     # Avoid text overlapping

#shows the top variables contributing to the 1st and 2nd axis
#red dotted line shows expected average contribution. 
#If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%. 
#For a given component, a variable with a contribution larger than this cutoff could be considered as important in contributing to the component.
fviz_contrib(allsites.pca, 
             choice = c("var"), 
             axes = 1:2, #both PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10 
fviz_contrib(allsites.pca, 
             choice = c("var"), 
             axes = 1, #1st PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10
fviz_contrib(allsites.pca, 
             choice = c("var"), 
             axes = 2, #2nd PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10
fviz_contrib(allsites.pca, 
             choice = c("var"), 
             axes = 3, #3rd PCA axis 
             sort.val = c("desc"), #descending order 
             top = 10) #show top 10

#variable loadings table
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}

loadings <- allsites.pca$rotation
sdev <- allsites.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
print(round(var.coord[, 1:3], 3))

#*HCPC cluster analysis ####
pca_result<-PCA(SWgrab_allsites_dat[, -1], scale.unit = TRUE, graph=FALSE) # then data are scaled to unit variance
allsites.hc<-HCPC(pca_result, graph = TRUE) #can click on where to cut graph

fviz_dend(allsites.hc, 
          show_labels = FALSE,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco"           # Rectangle color
          #labels_track_height = 0.8      # Augment the room for labels
)

#make a color pallete 
fviz_cluster(allsites.hc,
             repel = FALSE,            #  label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal()
)

## clust number is assigned to each observation 
clusters<-allsites.hc$data.clust
clusters$siteID<-SWgrab_allsites_dat$siteID #add back in the siteID to see what sites are in each cluster

##look at the variables that explain clustering the most 
allsites.hc$desc.var$quanti
##look at some sites within each cluster, for each cluster, the top 5 closest individuals to the cluster center is shown
allsites.hc$desc.ind$para

#find out what sites are in each cluster
frequency_clust<-rename(count(clusters, clust, siteID), Freq = n)

#new figure 
fviz_pca_ind(allsites.pca, label="none", 
             habillage=clusters$clust,  # group by cluster
             addEllipses=FALSE ,
             palette = c("darkblue", "#F6Be00",  "#552586", "grey", "darkgreen", "darkred"), 
             title = "surface water")


cowplot::save_plot("Figures/vectors.png", Figure_vectors, base_width = 7,
                   base_height = 4, dpi=600)


#potential plot for pub
#PCA_biplot_GW_reducedvar +xlim(-3, 4.5) +scale_y_reverse(limits=c(1.5,-4)) +scale_color_ucscgb() +scale_fill_ucscgb() +border(color = "black", size = 0.8, linetype = NULL) + theme(
  # Hide panel borders and remove grid lines
 # panel.border = element_blank(),
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  #axis.ticks = element_blank(),
  #axis.text.x = element_blank(),
  #axis.text.y = element_blank(),
  # Change axis line
#  axis.line = element_line(colour = "black")
#)


#### trying to figure out varimax rotation #### 
#using a package called "psych"
#rotation 
pc <- psych::principal(SWgrab_allsites_dat[, -1], 4,rotate="varimax") #have to chose the # of components to extract 
print(pc$scores[1:5,])  # Scores returned by principal()

#compare to no rotation 
pc2 <- psych::principal(SWgrab_allsites_dat[, -1], 4,rotate="none") #have to chose the # of components to extract 
print(pc2$scores[1:5,])  # Scores returned by principal()

print(pc$loadings)
print(pc2$loadings)

# using  the prcomp output 
allsites.pca <- prcomp(SWgrab_allsites_dat[, -1], center=T, scale=T)
rawLoadings     <- allsites.pca$rotation[,1:4] %*% diag(allsites.pca$sdev, 4, 4)
scores <- scale(allsites.pca$x[,1:4]) %*% varimax(rawLoadings)$rotmat
print(scores[1:5,])     

#other option using prcomp output 
rawLoadings     <- allsites.pca$rotation[,1:4] %*% diag(allsites.pca$sdev, 4, 4)
rotatedLoadings <- varimax(rawLoadings)$loadings
invLoadings     <- t(pracma::pinv(rotatedLoadings)) #have to install the pracma package to run this line 
scores          <- scale(SWgrab_allsites_dat[, -1]) %*% invLoadings
print(scores[1:5,])   
