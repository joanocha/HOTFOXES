setwd("~/Documents/work/PhD/Tasks/Task2_Phenotyping/Analysis/phyenotypic_Analyses/R/myPhysiology/META/")

library(ggplot2)  # ggplot() fortify()
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
library(raster)  # intersect()
library(ggsn)  # north2() scalebar()
library(rworldmap)  # getMap()
library(sf)
library(raster)
library(dplyr)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
require(maptools)
library(maps)
library("ggspatial")
shpdata_Rueppell<- readShapeLines(fn = "aridity_index/vulpes/vrueppellii/data_0.shp")
shpdata_Pale<- readShapeLines(fn = "aridity_index/vulpes/vpallida/data_0.shp")
shpdata_Fennec<- readShapeLines(fn = "aridity_index/vulpes/vzerda/data_0.shp")
shpdata_Red<- readShapeLines(fn = "aridity_index/vulpes/vvulpes_2021_global/Vulpes_global_revised.shp")
shpdata_Red_NA<- readShapeLines(fn = "aridity_index/vulpes/vvulpes_2021_NAfrica/Vulpes_Africa.shp")
shpdata_Red_EU<- readShapeLines(fn = "aridity_index/vulpes/vvulpes_2021_Eurasia/Vulpes_Eurasia_revised.shp")

world = world
world = world %>% filter(continent != 'Antarctica') 
Europe= world %>% filter(is.na(geom) | subregion == 'Southern Europe')
Africa= world %>% filter(is.na(geom) | continent == 'Africa')
Asia= world %>% filter(is.na(geom) | continent == 'Asia')
ai <- raster(readGDAL("aridity_index/AI_annual_old/ai_yr/w001001.adf"))
plot(ai)
starbucks <- read.csv("allvulpesDataCoords.csv")
starbucksR = starbucks  %>% filter(Species == "V. rueppellii")
starbucksRF_NA = starbucks  %>% filter(Region == "North Africa")
starbucksRF_EU = starbucks  %>% filter(Region == "Europe")
starbucksRF_ME = starbucks  %>% filter(Region == "Middle East")
starbucksF = starbucks  %>% filter(Species == "V. zerda")
starbucksP = starbucks  %>% filter(Species == "V. pallida")
colours <- colours <- c('#edc9af','#C1E1C1')

#--- Plot the data ---#
ppt <- stack(x = 'ppt.nc')
pet <- stack(x = 'pet.nc')
ppt_mean <- calc(ppt, # RasterStack object
                 fun = mean, # Function to apply across the layers
                 na.rm = TRUE)
pet_mean <- calc(pet,
                 fun = mean, 
                 na.rm = TRUE)
ext <- extent(c(xmin = -180, xmax = 180, 
                ymin = -60, ymax = 90))
ppt_mean <- crop(x = ppt_mean, 
                 y = ext)
pet_mean <- crop(x = pet_mean, 
                 y = ext)
aridity_index <- overlay(x = ppt_mean, # Raster object 1
                         y = pet_mean, # Raster object 2
                         fun = function(x, y){return(x / y)}) # Function to apply
aridity_index_matrix <- rasterToPoints(aridity_index)
aridity_index_df <- as.data.frame(aridity_index_matrix)
saturation<-0.2
aridity_index_df$layer[aridity_index_df$layer == "Inf"] <- saturation
aridity_index_df$layer[aridity_index_df$layer > saturation ] <- saturation

Figure1B<-ggplot(world)+geom_sf(fill = "#E6E6E6", colour="#E6E6E6")  + 
  geom_polygon(data = shpdata_Rueppell,aes(x = long, y = lat, group = group), fill = "#e31a1c", color="#e31a1c", alpha = 0.7, size = 0.2, show.legend = TRUE) +
  geom_polygon(data = shpdata_Red_NA,aes(x = long, y = lat, group = group), fill = "#ff7f00", color="#ff7f00", alpha = 0.4, size = 0.2, show.legend = FALSE) +
  #geom_polygon(data = shpdata_Red_EU,aes(x = long, y = lat, group = group), fill = "#33A02C", color="#33A02C", alpha = 0.3, size = 0.3, show.legend = FALSE) +
  geom_polygon(data = shpdata_Red_EU,aes(x = long, y = lat, group = group), fill = "#ff7f00", color="#ff7f00", alpha = 0.4, size = 0.2, show.legend = FALSE) +
  #geom_polygon(data = shpdata_Red_ME,aes(x = long, y = lat, group = group), fill = "#ff7f00", color="#ff7f00", alpha = 0.3, size = 0.3, show.legend = FALSE) +
  geom_polygon(data = shpdata_Pale,aes(x = long, y = lat, group = group), fill = "#fb9a99", color="#fb9a99", alpha = 0.7, size = 0.2, show.legend = FALSE) +
  geom_polygon(data = shpdata_Fennec,aes(x = long, y = lat, group = group), fill = "#6C1817", color="#6C1817", alpha = 0.6, size = 0.2, show.legend = FALSE) +
  theme_void() + theme(legend.position = "bottom") 
ggsave("Figure_1B.pdf", Figure1B, width = 5, height = 2, dpi=300)

Figure1B_zoom<-ggplot(world) +
  geom_raster(data = aridity_index_df, aes(y = y,
                                           x = x,
                                           fill = layer)) +
  annotation_scale(location = "bl", width_hint = 0.1) + 
  coord_sf(xlim = c(-20, 70), ylim = c(3, 50), expand = FALSE) + xlab("") + ylab("") +
  scale_fill_gradient(low = "#edc9af", high = "#E5E5E5") + #light grey is #E5E5E5 and arid color is #edc9af +
  geom_point(data = starbucksR, aes(x = Longitude, y = Latitude), colour = "black", fill="#e31a1c", size = 1, shape=21) +
  geom_point(data = starbucksRF_NA, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1, shape = 21) +
  geom_point(data = starbucksRF_ME, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1, shape = 21) +
  geom_point(data = starbucksRF_EU, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1, shape = 21) + 
  geom_point(data = starbucksP, aes(x = Longitude, y = Latitude), colour = "black", fill="#fb9a99", size = 1, shape = 21) + 
  geom_point(data = starbucksF, aes(x = Longitude, y = Latitude), colour = "black", fill="#6C1817", size = 1, shape = 21) + 
  theme_void() + theme(legend.position = "bottom") 
ggsave("Figure_1B_zoomIn.pdf", Figure1B_zoom, width = 5, height = 3, dpi=300)

                   
Figure2A_zoom<-ggplot(world) +
  geom_raster(data = aridity_index_df, aes(y = y,
                                           x = x,
                                           fill = layer)) +
  annotation_scale(location = "bl", width_hint = 0.1) + 
  coord_sf(xlim = c(-20, 70), ylim = c(3, 50), expand = FALSE) + xlab("") + ylab("") +
  scale_fill_gradient(low = "#edc9af", high = "#E5E5E5") + #light grey is #E5E5E5 and arid color is #edc9af +
  geom_point(data = starbucksR, aes(x = Longitude, y = Latitude), colour = "black", fill="#e31a1c", size = 1.5, shape=21) +
  geom_point(data = starbucksRF_NA, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1.5, shape = 21) +
  geom_point(data = starbucksRF_ME, aes(x = Longitude, y = Latitude), colour = "black", fill="#33a02c", size = 1.5, shape = 21) +
  geom_point(data = starbucksRF_EU, aes(x = Longitude, y = Latitude), colour = "black", fill="#1f78b4", size = 1.5, shape = 21) + 
  theme_void() + theme(legend.position = "bottom") 
ggsave("Figure_2A_zoomIn.pdf", Figure2A_zoom, width = 4, height = 2, dpi=300)
#ggsave("rueppellredfennecRed.pdf", plot, width = 5, height = 4, units = "in", dpi=300)


Figure4A_right<-ggplot(world) +
  geom_raster(data = aridity_index_df, aes(y = y,
                                           x = x,
                                           fill = layer)) +
  annotation_scale(location = "bl", width_hint = 0.1) + 
  coord_sf(xlim = c(-20, 70), ylim = c(3, 50), expand = FALSE) + xlab("") + ylab("") +
  scale_fill_gradient(low = "#edc9af", high = "#E5E5E5") + #light grey is #E5E5E5 and arid color is #edc9af +
  geom_point(data = starbucksRF_NA, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1.5, shape = 21) +
  theme_void() + theme(legend.position = "bottom") 
ggsave("Figure4A_right.pdf", Figure4A_right, width = 4, height = 2, dpi=300)

Figure4A_left<-ggplot(world) +
  geom_raster(data = aridity_index_df, aes(y = y,
                                           x = x,
                                           fill = layer)) +
  annotation_scale(location = "bl", width_hint = 0.1) + 
  coord_sf(xlim = c(-20, 70), ylim = c(3, 50), expand = FALSE) + xlab("") + ylab("") +
  scale_fill_gradient(low = "#E5E5E5", high = "#E5E5E5") + #light grey is #E5E5E5 and arid color is #edc9af +
  #geom_point(data = starbucksRF_NA, aes(x = Longitude, y = Latitude), colour = "black", fill="#ff7f00", size = 1.5, shape = 21) +
  theme_void() + theme(legend.position = "bottom") 
ggsave("Figure4A_left.pdf", Figure4A_left, width = 4, height = 2, dpi=300)

