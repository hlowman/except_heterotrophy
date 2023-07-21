# Clean up variables for literature annual site summaries and extract relevant 
# ones for figures and analyses.
# A Carter

library(tidyverse)
library(nhdplusTools)
library(mapview)
library(sf)

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/literature_streams_watershed_summary_data_1.csv')

# Condense categories and subset to useful columns ####
# Combine streamcat similar land use categories:
sc_vars <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')

nlcd_vars <- sc_vars %>% filter(grepl('^NLCD201[1,6]', Data.Location),
                              grepl('Ws$', Metric.Name)) %>%
  select(-Data.Location)
# Groups: Urban (low, medium, high-intensity and open space land use)
#         Agriculture (crop, hay land cover)
#         Forest (deciduous, evergreen, mixed forest)
#         Barren (barren, ice/snow)
#         Water (open water)
#         wetland (herbaceous, woody wetland)
#         grassland (shrub/scrub, grassland/herbaceous)

d <- dat %>% 
  mutate(NLCD_PctAgriculture = case_when(year < 2016 ~ PctCrop2011Ws + PctHay2011Ws,
                                year >= 2016 ~ PctCrop2016Ws + PctHay2016Ws),
         NLCD_PctForest = case_when(year < 2016 ~ PctDecid2011Ws + 
                                      PctConif2011Ws + PctMxFst2011Ws,
                                    year >=2016 ~ PctDecid2016Ws +
                                      PctConif2016Ws + PctMxFst2016Ws),
         NLCD_PctBarren = case_when(year < 2016 ~ PctBl2011Ws + PctIce2011Ws,
                                    year >= 2016 ~PctBl2016Ws + PctIce2016Ws),
         NLCD_PctWater = case_when(year < 2016 ~ PctOw2011Ws,
                                   year >= 2016 ~ PctOw2016Ws),
         NLCD_PctWetland = case_when(year < 2016 ~ PctHbWet2011Ws + PctWdWet2011Ws,
                                     year >= 2016 ~ PctHbWet2016Ws + PctWdWet2016Ws),
         NLCD_PctGrassland = case_when(year < 2016 ~ PctShrb2011Ws + PctGrs2011Ws,
                                       year >= 2016 ~ PctShrb2016Ws + PctGrs2016Ws)) %>%
  select(-any_of(unique(nlcd_vars$Metric.Name))) 

LU  <- select(d, starts_with('NLCD'))
cats <- apply(LU, 1, which.max)
cats[sapply(cats, function(x) ! length(x))] <- NA
d$NLCD_LUCat <- str_match(names(unlist(cats)), '^NLCD_Pct([A-Za-z]+)')[,2]

# Nitrogen deposition and application through agriculture:
# only keep the summarized inorganic N deposition (in Kg N/Ha)
d <- d %>% select(-SN_2008Ws, -NO3_2008Ws, -NH4_2008Ws) %>%
  rename(Inorg_N_WetDep_kgNhayr_2008 = InorgNWetDep_2008Ws)

# only keep the sum of all organic fertilizer (manure + bio fixation)
# and the sum of all inorganic fertilizer
d <- mutate(d, Inorg_N_fert_kgNhayr = FertWs,
            Org_N_fert_kgNhayr = CBNFWs + ManureWs) %>%
  select(-FertWs, -CBNFWs, -ManureWs)

# Combine Dam data
d <- mutate(d, Dam_densityperkm2 = (DamDensWs + NABD_DensWs)/2,
            Dam_total_vol_m3km2 = (DamNIDStorWs + NABD_NIDStorWs)/2,
            Dam_normal_vol_m3km2 = (DamNrmStorWs + NABD_NrmStorWs)/2) %>%
  select(-starts_with(c('DamN', 'DamD', 'NABD')))

# Waste point srcs
d <- mutate(d, Waste_point_srcs_perkm2 = SuperfundDensWs + NPDESDensWs + TRIDensWs) %>% 
  select(-SuperfundDensWs, -NPDESDensWs, -TRIDensWs)

# Combine data from the paired sites
w <- which(is.na(d[which(d$River == 'Polecat Creek'),] ))
d[which(d$River == 'Polecat Creek'), w] <- d[which(d$River == 'Big Springs'),w]

w <- which(is.na(d[which(d$River == 'Oakbrook Creek'),] ))
d[which(d$River == 'Oakbrook Creek'), w] <- d[which(d$River == 'Elijas Creek'),w]

lit_dat <- d %>%
  rename(Width = width_m, ws_area_km2 = NHD_TOTDASQKM_corr) %>%
  mutate(ratio_WA = Width/1000/sqrt(ws_area_km2))


# Calculate stream lengths below dams for each site using NHD
################################################################################
WGS84 = 4326
dams <- read_csv('data_ignored/nabd_dams.csv') %>% 
  sf::st_as_sf(coords = c('lon', 'lat'), crs = WGS84)
sites <- filter(lit_dat, !(River %in% c('Big Springs', 'Elijas Creek', 'La Choza')))
# functions for working with dams dataset:
COMID <- sites$COMID[1]
VPU <- sites$VPU[1]

a <- nhdplusTools::map_nhdplus(outlets = list(COMID),
                               flowline_only = TRUE)
flowline <- nhdplusTools::navigate_nldi(list(featureSource = 'comid',
                                             featureID = COMID),
                                        mode = 'upstreamTributaries',
                                        data_source = '')
sf::st_as_sf(sites[1,], coords = c('Longitude', 'Latitude'), crs = WGS84) %>%
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

lit_dat <- left_join(lit_dat, select(sites, nwis_gage,
                          total_flow_length, connected_flow_length), 
          by = 'nwis_gage') %>%
    mutate(drainage_density = total_flow_length/ws_area_km2,
           drainage_density_connected = connected_flow_length/ws_area_km2)

lit_dat <- lit_dat %>%
  mutate(width_to_area = Width/sqrt(ws_area_km2),
         Stream_PAR_sum = Stream_PAR_sum/365/1000, # convert units to match Bernhardt (mol/m2/d)
         MOD_ann_NPP = MOD_ann_NPP * 1000)         # kg C to g C /m2/y
           

write_csv(lit_dat, 'data_working/literature_streams_watershed_summary_data.csv')
# Pull out relevant variables
# first, fill in needed but missing values:
lit_dat <- read_csv('data_working/literature_streams_watershed_summary_data.csv')
# From US climate data.com: 
lit_dat$PrecipWs[lit_dat$River == 'SF Humboldt River'] <- 251.7 # mm
# From wikipedia:
lit_dat$ElevWs[lit_dat$River == 'SF Humboldt River'] <- 1544 # meters

lit_sum <- lit_dat %>%
  filter(!(River %in% c('Big Springs', 'Elijas Creek', 'La Choza'))) %>%
  # mutate(PAR_sum = PAR_sum * 0.001/365) %>%
  select(River, lat = Latitude, lon = Longitude, Type1, Type2, nwis_gage, 
         Width, PAR_kurt, width_to_area, PrecipWs, MOD_ann_NPP,
         RBI, Disch_cv, Disch_amp, Disch_ar1, drainage_density_connected, 
         ElevWs, Stream_PAR_sum, TmeanWs, max_interstorm) 

write_csv(lit_sum, 'data_working/literature_streams_data_for_PCA.csv')
