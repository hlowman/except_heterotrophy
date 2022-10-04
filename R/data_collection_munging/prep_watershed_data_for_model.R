# Clean up variables for annual site summaries and combine similar watershed characteristics
# A Carter
# March 2022

library(tidyverse)

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/annual_summary_data.csv')
ws_dat <- read_csv('data_356rivers/watershed_summary_data.csv')

# Condense categories and subset to useful columns ####
# Combine streamcat similar land use categories:
sc_vars <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')
data.frame(Metric.Name = colnames(dat)) %>%
  inner_join(sc_vars) %>% select(-Data.Location)

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
  mutate(NLCD_PctUrban = case_when(year < 2016 ~ PctUrbHi2011Ws + PctUrbMd2011Ws 
                                   + PctUrbLo2011Ws + PctUrbOp2011Ws,
                                   year >= 2016 ~ PctUrbHi2016Ws + PctUrbMd2016Ws 
                                   + PctUrbLo2016Ws + PctUrbOp2016Ws),
         NLCD_PctAgriculture = case_when(year < 2016 ~ PctCrop2011Ws + PctHay2011Ws,
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

# compare to IGBP classification:
d %>% select(starts_with('NLCD_LUC'), IGBP_LU_category) %>% data.frame()

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

# Other vars that could be grouped:
geo_vars <- sc_vars %>% filter(grepl('^Lithology', Data.Location)) %>% 
  select(-Data.Location) %>% data.frame()

# Variables that change over years (lu, imp surface) could be paired with the correct years
d <- d %>%
  mutate(Pct_impcov = case_when(year < 2004 ~ PctImp2001Ws,
                                year < 2006 ~ PctImp2004Ws,
                                year < 2008 ~ PctImp2006Ws,
                                year < 2011 ~ PctImp2008Ws,
                                year < 2013 ~ PctImp2011Ws,
                                year < 2016 ~ PctImp2013Ws,
                                year < 2019 ~ PctImp2016Ws,
                                year >= 2019 ~ PctImp2019Ws)) %>%
  select(-starts_with('PctImp'))

# The only remaining variable with many missing values is LAI, and by extension
# terrestrial NPP and Streamlight. Interestingly, all sites that are missing
# this are classified as urban:
filter(d, is.na(LAI_mean)) %>% 
  select(year, site_name, IGBP_LU_category, NLCD_LUCat) %>% data.frame()

urb_LAI <- filter(d, IGBP_LU_category == 'urban')

plot(density(d$LAI_mean, na.rm = T), ylim = c(0, .8))
lines(density(urb_LAI$LAI_mean, na.rm = T), lty = 2)
lines(density(filter(d, NLCD_LUCat == 'Urban')$LAI_mean, na.rm = T), lty = 3)  
plot(density(d$Stream_PAR_sum, na.rm = T), ylim = c(0, .2))
lines(density(urb_LAI$Stream_PAR_sum, na.rm = T), lty = 2)
lines(density(filter(d, NLCD_LUCat == 'Urban')$Stream_PAR_sum, na.rm = T), lty = 3)  
# it also doesn't look like there is a clear trend of LAI or streamPAR with urban

# I will try two different analyses, one with these sites taken out and another
# with these variables taken out.

# remove character columns and identification (not covariate) columns
d <- select(d, -long_name, -Source, -coord_datum, -alt, -alt_datum, -site_type, 
            -COMID, -VPU, -ndays, -nyears, -NHD_near_dist_m, -US_state, -region,
            -NHD_HYDROSEQ, -NHD_STARTFLAG, -NHD_REACHCODE, -NHD_AREASQKM,
            -NHD_TOTDASQKM,-NHD_LENGTHKM, -NHD_MAXELEVSMO, -NHD_MINELEVSMO, 
            -NHD_SLOPELENKM, -NHD_AREASQKM_corr, -areal_corr_factor) %>%
  rename(ws_area_km2 = NHD_TOTDASQKM_corr)

# Subset to high coverage variables and sites ####
# remove variables with low coverage then sites that don't have all variables
d <- d %>%
  select(-starts_with(c('HydroATLAS', 'struct.', 'dvqcoefs')), -HYRIV_ID,
         -Azimuth, -tree_height, -slope_calc)

# remove the two sites not in streamcat
w <- which(apply(d, 1, function(x) sum(is.na(x))) >50)
d <- d[-w,] 

glimpse(d)
glimpse(ws_dat)

d <- left_join(d, select(ws_dat, site_name, ends_with('flow_length'), 
                            starts_with('drainage_density')))

apply(d, 2, function(x) sum(is.na(x)))

write_csv(d, 'data_working/across_sites_model_data.csv')

