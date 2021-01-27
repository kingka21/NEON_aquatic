###### MAP CODE ###### 
#load libraries
library(rgdal) # read in shapefiles 
library(ggplot2)
library(dplyr)
library(raster) # to use the crs function to check projection  
library(sp) #spatial data handling 
library(cowplot) #arranging maps
library(ggsn) #add scale bar and north arrow 

#load shapefile 
poly<-readOGR(dsn = "/Users/katelynking/Desktop/NEON_Domains", layer = "NEON_Domains")
crs(poly) #shows me the projection of the shapefiles so that I can project the same to the points 

#read in lat/lon of points 
sites<-read.csv("/Data/neon_sites_latlon.csv")
sites$group<-as.factor(sites$group)

## points over NEON domains map 
# this is super zoomed out because of Alaska and PR 
jpeg('Figures/site_map.jpeg',width = 6,height = 4,units = 'in',res=600)
site_map<-ggplot(sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=2) +
  geom_point(shape = 1,size = 2,colour = "black") + 
  geom_path(data= poly, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF"),
                     name='group')+
  theme_bw() + 
  theme(#axis.text = element_blank(),
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        #axis.title = element_blank(),
        legend.text=element_text(colour='black', size=7),
        legend.title = element_text(color = "black", size = 7),
        legend.position = c(0.15, 0.35)) +
  guides(color = guide_legend(override.aes = list(size=1.5))) #increase legend point size
site_map
dev.off()


#### select out just the 48 states, then Alaska then PR ####

###change to an sf object   
sf_poly<-sf::st_as_sf(poly)

#select out contiguous states
poly_contig<-dplyr::filter(sf_poly, DomainID == "1" | DomainID == "2" | DomainID =="3"| DomainID =="4"|
                              DomainID =="5" |DomainID =="6" |DomainID =="7"| 
                             DomainID == "8" | DomainID == "9" | DomainID == "10" | DomainID == "11" | 
                             DomainID =="12" | DomainID == "13" | DomainID == "14" | 
                             DomainID == "15"  | DomainID =="16"| DomainID =="17" 
                            )  

poly_contig$DomainID<-as.factor(poly_contig$DomainID)
poly_contig$DomainID <- droplevels(poly_contig$DomainID)
#select out contiguous sites 
contig_sites<-filter(sites, site != "OKSR") %>%
  filter( site != "CARI") %>%
  filter( site != "GUIL") %>%
  filter( site != "CUPE") 

#select out PR 
poly_PR<-dplyr::filter(sf_poly, DomainID == "4") 
poly_PR$DomainID<-as.factor(poly_PR$DomainID)
poly_PR$DomainID <- droplevels(poly_PR$DomainID)
poly_PR<-dplyr::filter(poly_PR, OBJECTID == "0") 
#select out PR sites 
PR_sites<-filter(sites, site == "GUIL" | site == "CUPE")


#select out ALASKA 
poly_AL<-dplyr::filter(sf_poly, DomainID == "18" | DomainID == "19" ) 
poly_AL$DomainID<-as.factor(poly_AL$DomainID)
poly_AL$DomainID <- droplevels(poly_AL$DomainID)
#select out AL sites 
AL_sites<-filter(sites, site == "OKSR" | site == "CARI")

#convert it back to a spatial polygons data frame
poly_contig_sp<-as(poly_contig, "Spatial")
poly_PR_sp<-as(poly_PR, "Spatial")
poly_AL_sp<-as(poly_AL, "Spatial")

# CONTIG MAP
contig_plot<-ggplot(contig_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group, shape=group), size=4) +
  #geom_point(shape=1, size = 2,colour = "black") +
  scale_shape_manual(values=c(1,2,3,4,5,6,7)) +
  geom_path(data= poly_contig_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF"),
                     name='group')+
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size=3))) + #increase legend point size
  coord_sf(xlim = c(-129, -65), #change the extent to map and zoom in on contiguous US
           ylim = c(20, 50)) + 
  north( x.min = -129, x.max=-65, y.min = 20, y.max=50, symbol=3, scale=.1, location = "topright", anchor = c(x = -63, y = 50) ) + # anchor can change location
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "topleft", st.size = 2,
                 x.min = -129, x.max=-65, y.min = 20, y.max=50)

# ALASKA MAP
al_plot<-ggplot(AL_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=2) +
  geom_path(data= poly_AL_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#20854EFF"))+
  theme_bw() + 
   theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright", st.size = 2, st.dist = 0.03, border.size = 0.5,
                 x.min = -170, x.max=-140, y.min = 55, y.max=75)


  
# Puerto Rico MAP
pr_plot<-ggplot(PR_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=2) +
  geom_path(data= poly_PR_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#7876B1FF")) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  coord_sf(xlim = c(-67.4, -65.5), #change the extent to map and zoom in on contiguous US
           ylim = c(17.8, 18.5)) +
  ggsn::scalebar(dist = 20, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright", anchor = c(x = -65.5, y = 17.9), st.size = 2, st.dist = 0.04, border.size = 0.5,
                 x.min = -67.4, x.max=-65.5, y.min = 17.8, y.max=18.5)

# inset map 
inset_map<-ggdraw() +
  draw_plot(contig_plot) +
  draw_plot(al_plot, x = 0.07, y = 0.27, width = 0.25, height = 0.25) + 
  draw_plot(pr_plot, x = 0.68, y = 0.25, width = 0.2, height = 0.25)

inset_map
ggsave(filename = "inset_map.png", 
       plot = inset_map,
       path = "Figures/", 
       width = 7, 
       height = 6,
       dpi = 150)

