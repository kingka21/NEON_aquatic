# PCA analysis code
#Created by KK  March 30, 2020
#Modified 6Aug2020 by JWEdmonds, 20Sep20  ## 
#modified by KK on 13Sep2020; 28Sep 2020, Oct 5, 2020, Oct 14 
#modified by JWE on 7Dec2020


install.packages("factoextra") #install the package if you have not used it before
install.packages('fmsb')
install.packages("moments")
install.packages("rcompanion")
install.packages("REdaS")
install.packages("parameters")
install.packages("RNHANES")
library(factoextra) #library for great visualization of PCA results
library(FactoMineR) #library to perform PCA and HCPC clusters
library(ggplot2) #library for plotting 
library(dplyr) #library for data maniupulation 
library(agricolae) #library for skewness test
library(rcompanion) #transformTukey function 
library(parameters) # check KMO function 
library(REdaS) # for Bartlett test of sphericity 
library(car) #leveneTest
library(corrplot)
library(ggpubr)
library(scales)
library(viridis)  
library(ggsci)
library(gridExtra)
library(fmsb)
library(ggpubr)
library(moments)
library(tidyverse)
library(RNHANES)

###################run PCA for surface water data with all sites #########################

#bring in most recent data, QAQC and all sites updated Sep 28 
sw_all<-read.csv('Data/surface_water_grab_QCQC_ph.csv', header = TRUE)
SWgrab_allsites_dat<-select(sw_all, -c(sampleID, collectDate, pH, TDS, Br)) %>% #remove columns I don't need
                     data.table:: setnames(old=c("NH4...N", "NO3.NO2...N", "initialSamplepH"), new=c("NH4N", "NO3NO2N", "pH"))

summary(SWgrab_allsites_dat$siteID)

#Remove rows with na
SWgrab_allsites_dat<-na.omit(SWgrab_allsites_dat) # rows where lab pH was not available (from 1504 to 1385)

#Remove one outlier value with a "bad" MN concentration
SWgrab_allsites_dat<-(SWgrab_allsites_dat [-245, ])


#Test of normality for each paramenter using Shapiro-Wilks test
#values >0.5 are normal 
shapiro.test(SWgrab_allsites_dat$Ca)
shapiro.test(SWgrab_allsites_dat$Cl)
shapiro.test(SWgrab_allsites_dat$DIC)
shapiro.test(SWgrab_allsites_dat$DOC)
shapiro.test(SWgrab_allsites_dat$F)
shapiro.test(SWgrab_allsites_dat$Fe)
shapiro.test(SWgrab_allsites_dat$K)
shapiro.test(SWgrab_allsites_dat$Mg)
shapiro.test(SWgrab_allsites_dat$Mn)
shapiro.test(SWgrab_allsites_dat$Na)
shapiro.test(SWgrab_allsites_dat$NH4N)
shapiro.test(SWgrab_allsites_dat$NO3NO2N)
shapiro.test(SWgrab_allsites_dat$Si)
shapiro.test(SWgrab_allsites_dat$SO4)
shapiro.test(SWgrab_allsites_dat$TDP)
shapiro.test(SWgrab_allsites_dat$pH)


#check data for skewness
skewness(SWgrab_allsites_dat$pH)

#### transform data and create a column for the transformed data ####

#Q-Q Plot is bent at top after transformation
SWgrab_allsites_dat$transCa<-transformTukey(
  SWgrab_allsites_dat$Ca,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is great after transformation
SWgrab_allsites_dat$transCl<-transformTukey(
  SWgrab_allsites_dat$Cl,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is bent at top after transformation
SWgrab_allsites_dat$transDIC<-transformTukey(
  SWgrab_allsites_dat$DIC,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)


#Q-Q Plot is great after transformation
SWgrab_allsites_dat$transDOC<-transformTukey(
  SWgrab_allsites_dat$DOC,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)


#Q-Q Plot is bent at bottom, badly, after transformation
SWgrab_allsites_dat$transF<-transformTukey(
  SWgrab_allsites_dat$F,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is bent at bottom, badly, after transformation
SWgrab_allsites_dat$transFe<-transformTukey(
  SWgrab_allsites_dat$Fe,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is bent at bottom a little after transformation
SWgrab_allsites_dat$transK<-transformTukey(
  SWgrab_allsites_dat$K,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is bent at top, badly, after transformation
SWgrab_allsites_dat$transMg<-transformTukey(
  SWgrab_allsites_dat$Mg,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q Plot is bent badly at bottom after transformation
SWgrab_allsites_dat$transMn<-transformTukey(
  SWgrab_allsites_dat$Mn,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks good after transformation
SWgrab_allsites_dat$transNa<-transformTukey(
  SWgrab_allsites_dat$Na,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks good after transformation
SWgrab_allsites_dat$transNH4N<-transformTukey(
  SWgrab_allsites_dat$NH4N,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks OK, bent a little at top after transformation
SWgrab_allsites_dat$transNO3NO2N<-transformTukey(
  SWgrab_allsites_dat$NO3NO2N,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks good after transformation
SWgrab_allsites_dat$transSi<-transformTukey(
  SWgrab_allsites_dat$Si,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks good after transformation
SWgrab_allsites_dat$transSO4<-transformTukey(
  SWgrab_allsites_dat$SO4,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)

#Q-Q plot looks bent at bottom after transformation
SWgrab_allsites_dat$transTDP<-transformTukey(
  SWgrab_allsites_dat$TDP,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 2,
  returnLambda = FALSE
)


#visual assessment of frequency of data binned
plotNormalHistogram(SWgrab_allsites_dat$transFe)

#check data for skewness in transformed data
skewness(SWgrab_allsites_dat, na.rm = TRUE)

#create dataframe with just values for transformed data
SWgrab_allsites_dat_trans<-select(SWgrab_allsites_dat, -c(Ca, Cl, DIC, DOC, F, Fe, K, Mg, Mn, Na, NH4N, NO3NO2N, Si, SO4, TDP))

#*run the Kaiser-Meyer-Olkin (KMO) Test - how suited the data are for factor analysis ####
#values <0.6 are not adequate. Values >0.8 are ideal
#checking both transformed and untransformed data
check_kmo(SWgrab_allsites_dat [, -1])
check_kmo(SWgrab_allsites_dat_trans [, -1])

#*Bartlett's Test Of Sphericity tests the hypothesis that your correlation matrix is an identity matrix ####
#values <0.05 indicate factor analysis is useful
bart_spher(SWgrab_allsites_dat_trans [,-1], use = c("pairwise.complete.obs"))

#*look at homogeneity variance  ####
#transform to long format
long<-pivot_longer(SWgrab_allsites_dat_trans, 
                   cols= pH:transTDP,
                   names_to = 'analyte', 
                   values_to = 'value')

#Levene's Test and Fligner-Killeen test of homogeneity of variance
#Levene's test is statistically significant, then the null hypothesis, that the groups have equal variances, is rejected. i.e. you have unequal variances 
car::leveneTest(value ~ analyte, data= long)
fligner.test(value ~ analyte, data= long)

#convert values to z-scores and test homogeneity of variance
scaled<- scale(SWgrab_allsites_dat_trans[,-1], scale=TRUE) #sample mean zero, and scaled to have sample standard deviation one.
scaled<-as.data.frame(scaled)
longscale<-pivot_longer(scaled, 
                        cols= pH:transTDP,
                        names_to = 'analyte', 
                        values_to = 'value')
longscale$analyte<-as.factor(longscale$analyte)
car::leveneTest(value ~ analyte, data= longscale)
fligner.test(value ~ analyte, data= longscale)

#visualization of scatter in datasets
ggplot(longscale, aes(x = analyte, y = value)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Transformed and Scaled values for each analyte")

#Visual representation of variable's distribution
ggdensity(SWgrab_allsites_dat_trans, x = "pH", fill = "lightgray") +
  scale_x_continuous() +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

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
                    

##### run PCA ####
allsites.pca_trans<-prcomp(SWgrab_allsites_dat_trans [, -1], scale=TRUE) #perform PCA without the siteID and without line 268, a CUPE site with outlier Mn value

#alternative PCA command with retx=TRUE
allsites.pca_trans_ROTATE<-prcomp(SWgrab_allsites_dat_trans[, -1], retx=TRUE, scale=TRUE)

#rotation 
#rawLoadings     <- allsites.pca$rotation[,1:4] %*% diag(allsites.pca$sdev, 4, 4)
#rotatedLoadings <- varimax(rawLoadings)$loadings
#scores <- scale(allsites.pca$x[,1:4]) %*% varimax(rawLoadings)$rotmat
#print(scores[1:5,])  

#eigenvalues
fviz_eig(allsites.pca_trans)  ###eigen values associated with each PC (how much does each explain the variation in the data)
eig.val_trans<-get_eigenvalue(allsites.pca_trans) ## table of the eigenvalues 

# graph results with the site grouping, all variables 
fviz_pca_ind(allsites.pca_trans, label="none", 
                                 habillage=SWgrab_allsites_dat_trans$siteID,  # group by the site
             addEllipses=TRUE, ellipse.level=0.75)
#PCAallsites_allvar +ylim(-4,5)+xlim(-7.5, 4.0) #adjust axis if needed

#graph vectors of PCA
vectors_SW<-fviz_pca_var(allsites.pca_trans,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             title = "Surface Water PCA")     # Avoid text overlapping

vectors_SW

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

loadings <- allsites.pca_trans$rotation
sdev <- allsites.pca_trans$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
print(round(var.coord[, 1:3], 3))

#*HCPC cluster analysis ####
pca_result_trans<-PCA(SWgrab_allsites_dat_trans[, -1], scale.unit = TRUE, graph=FALSE) # then data are scaled to unit variance
allsites.hc_trans<-HCPC(pca_result_trans, graph = TRUE) #can click on where to cut graph

# makes dendrogram figure - takes a bit to compute
dendrogram_figure<-fviz_dend(allsites.hc_trans, 
          show_labels = FALSE,                     # Label size
          palette = "nejm",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "nejm" ,          # Rectangle color
          title = '')

#shows the clusters 
fviz_cluster(allsites.hc_trans,
             repel = FALSE,            #  label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "nejm",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal()
)

## clust number is assigned to each observation 
clusters_trans<-allsites.hc_trans$data.clust
clusters_trans$siteID<-SWgrab_allsites_dat_trans$siteID #add back in the siteID to see what sites are in each cluster

##look at the variables that explain clustering the most 
allsites.hc_trans$desc.var$quanti

##look at the axes that clusters fall on
allsites.hc_trans$desc.axes$quanti

##look at some sites within each cluster, for each cluster, the top 5 closest individuals to the cluster center is shown
allsites.hc_trans$desc.ind$para

#find out what sites are in each cluster
frequency_clust_trans<-rename(count(clusters_trans, clust, siteID), Freq = n)

#new figure for drawing ellipses by groups 
cluster_figure<-fviz_pca_ind(pca_result_trans, label="none", 
             habillage=clusters_trans$clust,  # group by cluster
             addEllipses=TRUE ,
             palette = "nejm", 
             title = '') +  
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"), 
                   legend.position = "top",  legend.box = "horizontal") + 
            xlab("Component 1 (35.2%)") + 
            ylab(" Component 2 (12.4%)")  
  
             
#cluster_figure +  guides(palette = guide_legend(nrow = 1))

Fig3<-cowplot::plot_grid(cluster_figure, dendrogram_figure,
                         labels = c("a", "b")
)
cowplot::save_plot("Figures/Fig3.png", Fig3, base_width = 8,
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
#install.packages("psych")
library(psych)
#rotation 
pc <- psych::principal(SWgrab_allsites_dat[, -1], 4,rotate="varimax") #have to chose the # of components to extract 
print(pc$scores[1:5,])  # Scores returned by principal()

#compare to no rotation 
pc2 <- psych::principal(SWgrab_allsites_dat[, -1], 5,rotate="none") #have to chose the # of components to extract 
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
