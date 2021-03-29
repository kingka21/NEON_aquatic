###### MAP CODE ###### 
#written by K King March, 2021
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
crs(poly) #shows me the projection of the shapefiles 
#project to USA Contiguous albers equal area
us.aea<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
poly.aea <- spTransform(poly, us.aea)
crs(poly.aea)

#read in lat/lon of points 
sites<-read.csv("Data/neon_sites_latlon.csv")
sites$group<-as.factor(sites$group)
#projecting
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
sites.proj <- SpatialPoints(coords=sites[,c(3,2)], #lon then lat 
                       proj4string=wgs1984.proj)
# Reproject to Albers Equal Area Projection so it matches the US Map
sites.aea <- spTransform(sites.proj, us.aea)
sites.plot<-as(sites.aea, "data.frame") #change to a data frame for ggplot 
sites.plot$group<-as.factor(sites$group)
sites.plot$site<-sites$site

## points over NEON domains map 
# this is super zoomed out because of Alaska and PR 
jpeg('Figures/site_map.jpeg',width = 6,height = 4,units = 'in',res=600)
site_map<-ggplot(sites.plot, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=2) +
  geom_point(shape = 1,size = 2,colour = "black") + 
  geom_path(data= poly.aea, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
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

### try to limit it to Alaska, US, PR
site_map<-ggplot(sites.plot, aes(x=lon,y=lat))+
  geom_point(aes(colour=group, shape=group), size=3) +
  #geom_point(shape = 1,size = 2,colour = "black") + 
  scale_shape_manual(values=c(16,17,3,4,18,25,15)) +
  geom_path(data= poly.aea, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF"),
                     name='group')+
  theme_bw() + 
  theme(axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.text=element_text(colour='black', size=7),
    legend.title = element_text(color = "black", size = 7),
    legend.position = c(0.15, 0.35)) +
  guides(color = guide_legend(override.aes = list(size=1.5))) + #increase legend point size
  coord_sf(xlim = c(-4000000, 2500000), #change the extent to map and zoom in on contiguous US
           ylim = c(-1500000, 4500000)) + 
  north(x.min = -4000000, x.max=2500000, y.min = -1500000, y.max=4500000,  symbol=3, scale=.1, location = "topright" ) + # anchor can change location # anchor = c(x = -63, y = 50) 
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = FALSE, model = "GRS80", location = "bottomleft", st.size = 2,  
                 x.min = -4000000, x.max=2500000, y.min = -1500000, y.max=4500000) #
                 
site_map


# Puerto Rico MAP
pr_plot<-ggplot(PR_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group, shape=group), size=4) +
  scale_shape_manual(values=c(18)) +
  geom_path(data= poly_PR_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#7876B1FF")) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  #coord_sf(xlim = c(3100000, 3270000), #change the extent to map 
    #       ylim = c(-1675000, -1575000)) +
  ggsn::scalebar(data=PR_sites, dist = 100, dist_unit= "mi", transform = TRUE, model = "GRS80", location = "bottomright", st.size = 3, st.dist = 0.06, height = 0.1) 
                 #x.min = 3100000, x.max=3270000, y.min = -1675000, y.max=-1575000)

pr_plot

# inset map with symbols 
jpeg('Figures/site_map.jpeg',width = 6,height = 4,units = 'in',res=600)
inset_map<-ggdraw() +
  draw_plot(site_map) +
  draw_plot(pr_plot, x = 0.5, y = 0.6,width = 0.2, height = 0.25)
inset_map
dev.off()
################################################################
#### select out just the 48 states, then Alaska then PR ####

###change to an sf object   
sf_poly<-sf::st_as_sf(poly.aea)

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
contig_sites<-filter(sites.plot, site != "OKSR") %>%
  filter( site != "CARI") %>%
  filter( site != "GUIL") %>%
  filter( site != "CUPE") 

#select out PR 
poly_PR<-dplyr::filter(sf_poly, DomainID == "4") 
poly_PR$DomainID<-as.factor(poly_PR$DomainID)
poly_PR$DomainID <- droplevels(poly_PR$DomainID)
poly_PR<-dplyr::filter(poly_PR, OBJECTID == "0") 
#select out PR sites 
PR_sites<-filter(sites.plot, site == "GUIL" | site == "CUPE")


#select out ALASKA 
poly_AL<-dplyr::filter(sf_poly, DomainID == "18" | DomainID == "19" ) 
poly_AL$DomainID<-as.factor(poly_AL$DomainID)
poly_AL$DomainID <- droplevels(poly_AL$DomainID)
#select out AL sites 
AL_sites<-filter(sites.plot, site == "OKSR" | site == "CARI")

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
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "topleft", anchor = c(x = -129, y = 50.5), st.size = 2,
                 x.min = -129, x.max=-65, y.min = 20, y.max=50)

# ALASKA MAP
al_plot<-ggplot(AL_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group, shape=group), size=4) +
  scale_shape_manual(values=c(1,4)) +
  geom_path(data= poly_AL_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#20854EFF"))+
  theme_bw() + 
   theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  coord_sf(xlim = c(-170, -139), #change the extent to map 
           ylim = c(55, 75)) +
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright",  anchor = c(x = -141, y = 56), st.size = 1.5, st.dist = 0.04, border.size = 0.5,
                 x.min = -170, x.max=-139, y.min = 55, y.max=75)


  
# Puerto Rico MAP
pr_plot<-ggplot(PR_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group, shape=group), size=4) +
  scale_shape_manual(values=c(18)) +
  geom_path(data= poly_PR_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#7876B1FF")) +
  theme_bw() +
 theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  coord_sf(xlim = c(-67.4, -65.5), #change the extent to map 
           ylim = c(17.8, 18.6)) +
  ggsn::scalebar(dist = 50, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright", anchor = c(x = -65.6, y = 17.9), st.size = 1.5, st.dist = 0.06, border.size = 0.25,
                 x.min = -67.4, x.max=-65.5, y.min = 17.8, y.max=18.5)

# inset map with symbols 
inset_map<-ggdraw() +
  draw_plot(contig_plot) +
  draw_plot(al_plot, x = 0.07, y = 0.22, width = 0.25, height = 0.25) + 
  draw_plot(pr_plot, x = 0.68, y = 0.25, width = 0.2, height = 0.25)

inset_map
ggsave(filename = "inset_map.png", 
       plot = inset_map,
       path = "Figures/", 
       width = 7, 
       height = 4,
       dpi = 150)

### map with colors only ####
# CONTIG MAP
contig_plot<-ggplot(contig_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=4) +
  geom_point(shape=1, size=4,colour = "black") +
  geom_path(data= poly_contig_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF"),
                     name='group')+
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size=3))) + #increase legend point size
  coord_sf(xlim = c(-129, -65), #change the extent to map and zoom in on contiguous US
           ylim = c(20, 50)) + 
  north( x.min = -129, x.max=-65, y.min = 20, y.max=50, symbol=3, scale=.1, location = "topright", anchor = c(x = -63, y = 50) ) + # anchor can change location
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "topleft", anchor = c(x = -129, y = 50.5), st.size = 2,
                 x.min = -129, x.max=-65, y.min = 20, y.max=50)

# ALASKA MAP
al_plot<-ggplot(AL_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=4) +
  geom_point(shape=1, size=4,colour = "black") +
  geom_path(data= poly_AL_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#BC3C29FF", "#20854EFF"))+
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  coord_sf(xlim = c(-170, -139), #change the extent to map 
           ylim = c(55, 75)) +
  ggsn::scalebar(dist = 500, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright",  anchor = c(x = -141, y = 56), st.size = 1.5, st.dist = 0.04, border.size = 0.5,
                 x.min = -170, x.max=-139, y.min = 55, y.max=75)



# Puerto Rico MAP
pr_plot<-ggplot(PR_sites, aes(x=lon,y=lat))+
  geom_point(aes(colour=group), size=4) +
  geom_point(shape=1, size=4,colour = "black") +
  geom_path(data= poly_PR_sp, aes(long,lat,group=group),colour='black', size=0.2) + coord_equal() +
  scale_color_manual(values=c("#7876B1FF")) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) + 
  coord_sf(xlim = c(-67.4, -65.5), #change the extent to map 
           ylim = c(17.8, 18.6)) +
  ggsn::scalebar(dist = 50, dist_unit= "km", transform = TRUE, model = "WGS84", location = "bottomright", anchor = c(x = -65.6, y = 17.9), st.size = 1.5, st.dist = 0.06, border.size = 0.25,
                 x.min = -67.4, x.max=-65.5, y.min = 17.8, y.max=18.5)

# inset map with symbols 
inset_map_colors<-ggdraw() +
  draw_plot(contig_plot) +
  draw_plot(al_plot, x = 0.07, y = 0.22, width = 0.25, height = 0.25) + 
  draw_plot(pr_plot, x = 0.68, y = 0.25, width = 0.2, height = 0.25)

inset_map_colors
ggsave(filename = "inset_map_colors.png", 
       plot = inset_map_colors,
       path = "Figures/", 
       width = 7, 
       height = 4,
       dpi = 150)
