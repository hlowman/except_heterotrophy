library(sf)
library(ggplot2)
library(tidyverse)
library(nhdplusTools)
library(nhdR)
library(sp)
library(rgdal)
library(mapview)

dam_coords <- read_csv('data_working/US_dams.csv')
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
WGS84 = 4326 #EPSG code for coordinate reference system
comid_from_point = function(lat, long, CRS) {
  pt = sf::st_point(c(long, lat))
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
#this calculates how far along a reach any given point falls. That way when we pull in
#watershed summary data for a reach, we can adjust it according to how much
#of the total upstream area actually contributes to the point in question.
# A value of 0 means upstream end; 1 means downstream end.
calc_reach_prop = function(VPU, COMID, lat, long, CRS, quiet=FALSE){

  if(! quiet){
    message(paste0('The nhdR package downloads NHDPlusV2 components to ',
                   nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
                   ' Fortunately, each component need only be downloaded once.'))
  }

  fl = nhdR::nhd_plus_load(vpu=VPU, component='NHDSnapshot',
                           dsn='NHDFlowline', approve_all_dl=TRUE)
  fl_etc = nhdR::nhd_plus_load(vpu=VPU, component='NHDPlusAttributes',
                               dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

  colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
  colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
  fl = fl[fl$COMID == COMID,]
  fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
                 by=c('COMID'='ComID'))

  pt = sf::st_point(c(long, lat))
  ptc = sf::st_sfc(pt, crs=CRS)
  ptct = sf::st_transform(ptc, crs=4269) #CRS for NAD 83
  x = try(suppressWarnings(nhdplusTools::get_flowline_index(fl, points=ptct)))
  if(inherits(x, 'try-error')) return(NA_real_)
  out = 1 - x$REACH_meas / 100 #0=upstream end; 1=downstream end

  return(out)
}

dam_coords$CRS <- 'WGS84'
dam_coords$COMID = NA_real_
dam_coords$VPU = NA_character_

for(i in 1:nrow(dam_coords)){
  dam_coords$COMID[i] <- try(comid_from_point(dam_coords$lat[i],
                                              dam_coords$lon[i], WGS84))
  dam_coords$VPU[i] <- try(vpu_from_point(dam_coords$lat[i],
                                          dam_coords$lon[i], WGS84))
  x <- try(calc_reach_prop(dam_coords$VPU[i], dam_coords$COMID[i],
                           dam_coords$lat[i], dam_coords$lon[i], WGS84,
                           quiet = TRUE))
  if(inherits(x, 'try-error')) x <- NA
  if(length(x)==0) x <- NA
  dam_coords$reach_proportion[i] <- x

  if(i%%100 == 0)
    print(i/nrow(dam_coords))
}

# if the reach proportion cannot be calculated, assume the entire reach is included
dam_coords$reach_proportion[is.na(dam_coords$reach_proportion)] <- 1
write_csv(dam_coords, 'data_working/US_dams.csv')
# COMID <- unlist(mapply(comid_from_point, us_dams_latlon3$lat, us_dams_latlon3$lon, us_dams_latlon3$crs))

