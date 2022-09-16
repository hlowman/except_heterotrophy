library(sf)
library(ggplot2)
library(tidyverse)
library(nhdplusTools)
library(nhdR)
library(sp)
library(rgdal)
library(mapview)

us_dams <- readOGR("C:/Users/cbarbosa/OneDrive - University of Wyoming/Datasets/US_dams/US_dams.shp",stringsAsFactor = F)
watersheds <- read.csv("data_356rivers/watershed_summary_data.csv")

coords <- coordinates(us_dams)%>%
  as.data.frame(coords)%>%
  rename(lon = V1, lat = V2)

crs <- crs

us_dams_df <- as.data.frame(us_dams, xy = TRUE)

us_dams_latlon <- cbind(us_dams_df, coords)

us_dams_latlon3 <-us_dams_latlon %>%
  mutate(crs = WGS84)


#what is the crs and extent?
st_crs(us_dams)$proj4string
st_bbox(us_dams)

WGS84 = 4326 #EPSG code for coordinate reference system

#plotting us dams
ggplot()+
  geom_sf(data=us_dams, size= 3, color= "black", fill= "cyan1")+
  ggtitle("US dams")+
  coord_sf()

#using the same reference
longlat <- st_crs(us_dams)

#convert df to spatial object
watershed_summary_data_geo <- st_as_sf(watersheds,
                                   coords= c("lon", "lat"),
                                   crs=longlat)

class(watershed_summary_data_geo)
class(us_dams)

#plotting all
ggplot()+
  geom_sf(data = us_dams_latlon)+
  geom_sf(data = watershed_summary_data_geo, size= 3, color= "yellow")+
  ggtitle("Map of dams and watersheds locations")

mapview(us_dams, color= "black") + mapview(watershed_summary_data_geo, color = "yellow")


# site data must include latitude and longitude columns (decimal degrees)
# dataframe generated from filter_powell_estimates scipt
sites = read_tsv('data_356rivers/site_data.tsv')
NAD83 = 4269 #EPSG code for coordinate reference system
comid_from_point = function(lat, lon, CRS) {
  pt = sf::st_point(c(lon, lat))
  ptc = sf::st_sfc(pt, crs=CRS)
  COMID = nhdplusTools::discover_nhdplus_id(ptc)
  if(! length(COMID)) COMID = NA
  return(COMID)
}
vpu_from_point = function(lat, lon, CRS) {
  pt = sf::st_point(c(lon, lat))
  ptc = sf::st_sfc(pt, crs=CRS)
  VPU = nhdR::find_vpu(ptc)
  return(VPU)
}

#COMID <- unlist(mapply(comid_from_point, watersheds$lat, watersheds$lon, NAD83))

COMID <- unlist(mapply(comid_from_point, us_dams_latlon3$lat, us_dams_latlon3$lon, us_dams_latlon3$crs))

