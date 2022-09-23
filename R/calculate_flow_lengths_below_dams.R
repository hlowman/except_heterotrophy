# Calculate stream lengths below dams for each site using NHD

# 9/2022
# A Carter

library(tidyverse)
library(nhdplusTools)
library(mapview)
library(sf)
dams <- read_csv('data_working/US_dams.csv') %>% 
    st_as_sf(coords = c('lon', 'lat'))
sites <- read_csv('data_356rivers/watershed_summary_data.csv')

# functions for working with dams dataset:

a <- nhdplusTools::map_nhdplus(outlets = list(COMID),#dam_coords$COMID[i]),
                               flowline_only = TRUE)
flowline <- nhdplusTools::navigate_nldi(list(featureSource = 'comid',
                                             featureID = COMID),
                                        mode = 'upstreamTributaries',
                                        data_source = '')
st_as_sf(dam_coords[i,], coords = c('lon', 'lat'), crs = WGS84) %>%
  mapview::mapview()+ 
  mapview::mapview(flowline)

subset_dams_in_network <- function(COMID, dams){
  
    data = nhdR::nhd_plus_load(vpu =VPU, 
                               component = 'NHDPlusAttributes',
                               dsn = 'PlusFlowlineVAA',
                               approve_all_dl=TRUE)
  
    comids <- nhdplusTools::get_UT(data, sites$COMID[1])
    d <- filter(dams, COMID %in% comids)
    dam_comids <- c()
    for(i in 1:nrow(d)){
        dd <- nhdplusTools::get_UT(data, d$COMID[i])
        dd <- Filter(function(x) x != d$COMID[i], dd)
        dam_comids <- append(dam_comids, dd)
    }
    
    rp <- bind_rows(d, data.frame(COMID = sites$COMID[1],
                                 reach_proportion = 1-sites$reach_proportion[1])) %>%
      mutate(reach_proportion = 1 - reach_proportion) %>%
      select(ComID = COMID, reach_proportion)
    comids <- purrr::keep(comids, ~ !. %in% dam_comids)
    
    flines <- data %>% filter(ComID %in% comids) %>%
      left_join(rp, by = 'ComID') %>%
      mutate(reach_proportion = case_when(is.na(reach_proportion) ~ 1,
                                          TRUE ~ reach_proportion),
             sum_length = reach_proportion * LengthKM)
    
    flow_length = sum(flines$sum_length)
    # change this code to also get the Arbolate Sum distance of each subset
    }
