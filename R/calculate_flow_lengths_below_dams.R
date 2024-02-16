# Calculate stream lengths below dams for each site using NHD

# 9/2022
# A Carter

library(tidyverse)
library(nhdplusTools)
library(mapview)
library(sf)
WGS84 = 4326
dams <- read_csv('data_ignored/nabd_dams.csv') %>% 
    st_as_sf(coords = c('lon', 'lat'), crs = WGS84)
sites <- read_csv('data_356rivers/watershed_summary_data.csv')

# functions for working with dams dataset:
COMID <- sites$COMID[1]
VPU <- sites$VPU[1]

a <- nhdplusTools::map_nhdplus(outlets = list(COMID),
                               flowline_only = TRUE)
flowline <- nhdplusTools::navigate_nldi(list(featureSource = 'comid',
                                             featureID = COMID),
                                        mode = 'upstreamTributaries',
                                        data_source = '')
st_as_sf(sites[1,], coords = c('lon', 'lat'), crs = WGS84) %>%
  mapview::mapview()+ 
  mapview::mapview(flowline)


subset_dams_in_network <- function(site, dams){
  
    data = nhdR::nhd_plus_load(vpu =site$VPU, 
                               component = 'NHDPlusAttributes',
                               dsn = 'PlusFlowlineVAA',
                               approve_all_dl=TRUE)
    if('LengthKm' %in% colnames(data)) data = rename(data, LengthKM = LengthKm)
  
    comids <- nhdplusTools::get_UT(data, site$COMID)
    if(length(comids) < 1) return(list(fl = NA, al = NA))
    d <- filter(dams, COMID %in% comids) # subset to only include dams upstream of site
    dam_comids <- c()
    while(nrow(d) >= 1){
        cc = d$COMID[1]
        dd <- nhdplusTools::get_UT(data, cc)
        d <- filter(d, !(COMID %in% dd))
        dd <- Filter(function(x) x != cc, dd)
        dam_comids <- append(dam_comids, dd)
    }
    
    
    d <- filter(dams, COMID %in% comids) # subset to only include dams upstream of site
    
    if(nrow(d) != 0){
      rp <- bind_rows(d, data.frame(COMID = site$COMID,
                                   reach_proportion = site$reach_proportion)) %>%
        mutate(reach_proportion = 1 - reach_proportion) %>%
        select(ComID = COMID, reach_proportion)
    }else{
      rp = data.frame(ComID = site$COMID,
                      reach_proportion = 1- site$reach_proportion)
      }
    comids <- purrr::keep(comids, ~ !. %in% dam_comids)
    
    flines <- data %>% filter(ComID %in% comids) %>%
      left_join(rp, by = 'ComID') %>%
      mutate(reach_proportion = case_when(is.na(reach_proportion) ~ 1,
                                          TRUE ~ reach_proportion),
             sum_length = reach_proportion * LengthKM)
    
    flow_length = sum(flines$sum_length)
    arbolate_length = flines$ArbolateSu[flines$ComID == site$COMID] 
    
    return(list(fl = flow_length, 
                al = arbolate_length))
}

for(i in 1:nrow(sites)){
    lengths <- subset_dams_in_network(sites[i,], dams)
    sites$connected_flow_length[i] <- lengths$fl
    sites$total_flow_length[i] <- lengths$al
    
    if(i %% 5 == 0) print(i/nrow(sites))
}
sites <- sites %>%
    mutate(connected_flow_length = case_when(connected_flow_length < 0 ~ NA_real_,
                                             is.infinite(connected_flow_length) ~ NA_real_,
                                             TRUE ~ connected_flow_length),
           total_flow_length = case_when(total_flow_length < 0 ~ NA_real_,
                                         is.infinite(total_flow_length) ~ NA_real_,
                                         TRUE ~ total_flow_length),
           drainage_density = total_flow_length/NHD_TOTDASQKM_corr,
           drainage_density_connected = connected_flow_length/NHD_TOTDASQKM_corr)

write_csv(sites, 'data_356rivers/watershed_summary_data.csv')


# compare dam metrics:

plot(sites$DamDensWs, sites$connected_flow_length/sites$total_flow_length)
plot(sites$DamNIDStorWs, sites$DamDensWs, xlim = c(0, 1000000))

