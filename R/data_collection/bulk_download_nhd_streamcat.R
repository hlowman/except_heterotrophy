#batch summarizing of watershed data for the continental USA
#Mike Vlah (vlahm13@gmail.com)

#use these tools to:
#1. acquire NHDPlusV2 COMIDs based on lat/long
#2. acquire NHDPlusV2 VPU IDs based on lat/long
#3. determine the position of a site along its NHD reach and calculate reach
#   proportion, then adjust linear and areal summary data accordingly
#4. use COMID and VPU to acquire NHDPlusV2 data for your sites
#5. use COMID to acquire StreamCat data for your sites

#see NHDPlusV2 docs (1) and StreamCat variable list (2) for help.
#1. ftp://ftp.horizon-systems.com/NHDplus/NHDPlusV21/Documentation/NHDPlusV2_User_Guide.pdf
#2. ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/Documentation/VariableList-QuickReference.html

#NOTE: these tools necessarily download a lot of large datasets and store them
#in memory. keep an eye on your usage.
library(readr)
library(stringr)
library(plyr)
library(dplyr)
# remotes::install_github("USGS-R/nhdplusTools")
library(nhdplusTools)
# remotes::install_github("jsta/nhdR")
library(nhdR)
library(RCurl)
library(sf)
# library(MODISTools)

# 1. setup and helper functions ####

#set your working directory to the location of your site data file.
#   site data must include latitude and longitude columns (decimal degrees)
# setwd('C:/Users/Alice Carter/git/except_heterotrophy/data_ignored/')
sites = readRDS('site_data.rds')
NAD83 = 4269 #EPSG code for coordinate reference system

comid_from_point = function(lat, long, CRS) {
  pt = sf::st_point(c(long, lat))
  ptc = sf::st_sfc(pt, crs=CRS)
  COMID = nhdplusTools::discover_nhdplus_id(ptc)
  if(! length(COMID)) COMID = NA
  return(COMID)
}

vpu_from_point = function(lat, long, CRS) {
  pt = sf::st_point(c(long, lat))
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

#this acquires nhdplusv2 data for a single site by COMID.
#it's just a thin wrapper around nhdR::nhd_plus_load
nhdplusv2_from_comid = function(VPU, COMID, component, DSN, quiet=FALSE) {
  
  if(! quiet){
    message(paste0('The nhdR package downloads NHDPlusV2 components to ',
                   nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
                   ' Fortunately, each component need only be downloaded once.'))
  }
  
  data = nhdR::nhd_plus_load(vpu=VPU, component=component,
                             dsn=DSN, approve_all_dl=TRUE)
  
  colnames(data)[colnames(data) == 'ComID'] = 'COMID'
  colnames(data)[colnames(data) == 'ReachCode'] = 'REACHCODE'
  data = data[data$COMID == COMID,]
  
  return(data)
}

#this calls nhdplusv2_from_comid repeatedly to get data for all your sites.
#the dataframe must include COMID and VPU columns
nhdplusv2_bulk = function(site_df, nhdplusv2_sets, quiet=FALSE){
  
  nhdplus_data = data.frame()
  if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')
  
  for(j in 1:nrow(site_df)){
    for(i in 1:length(setlist)){
      print(paste(j, nhdplusv2_sets[[i]]))
      
      if(i == 1 || initerr){
        row_base = try(nhdplusv2_from_comid(site_df$VPU[j],
                                            site_df$COMID[j], names(setlist[i]), setlist[[i]],
                                            quiet=quiet))
        if('try-error' %in% class(row_base) || nrow(row_base) > 1){
          initerr = TRUE
          row_base = data.frame(COMID=site_df$COMID[j])
        } else {
          initerr = FALSE
        }
      } else {
        row_ext = try(nhdplusv2_from_comid(site_df$VPU[j],
                                           site_df$COMID[j], names(setlist[i]), setlist[[i]],
                                           quiet=quiet))
        if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
          row_base = left_join(row_base, row_ext)
        }
      }
      
    }
    
    if(nrow(row_base) > 1){
      row_base = data.frame(COMID=site_df$COMID[j])
    }
    
    nhdplus_data = rbind.fill(nhdplus_data, row_base)
  }
  
  return(nhdplus_data)
}

query_streamcat_datasets = function(keyword=NULL){
  
  ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
                  'NHDPlusLandscapeAttributes/StreamCat/States/')
  
  url_list = getURL(ftpdir, dirlistonly=TRUE)
  url_list = strsplit(url_list, split='\n')[[1]]
  
  if(! is.null(keyword)){
    url_list = url_list[grep(keyword, url_list, ignore.case=TRUE)]
  }
  
  return(url_list)
}

#this function acquires streamcat data for a single site by NHDPlusV2 COMID.
streamcat_from_comid = function(USstate, COMID, dataset){
  
  ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
                  'NHDPlusLandscapeAttributes/StreamCat/States/')
  zip_name = paste0(dataset, '_', USstate, '.zip')
  
  csv_name = gsub('.zip', '.csv', zip_name)
  temp = tempfile()
  download.file(paste0(ftpdir, zip_name), temp)
  data = read.csv(unz(temp, csv_name), stringsAsFactors=FALSE)
  data = data[data$COMID == COMID,]
  
  return(data)
}

#this calls streamcat_from_comid repeatedly to get data for all your sites
#the dataframe must include COMID and region columns, where "region" refers to
#each state's 2-letter abbreviation.
streamcat_bulk = function(site_df, streamcat_sets){
  
  streamcat_data = data.frame()
  if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')
  
  for(j in 1:nrow(site_df)){
    for(i in 1:length(streamcat_sets)){
      print(paste(j, streamcat_sets[i]))
      
      if(i == 1 || initerr){
        row_base = try(streamcat_from_comid(site_df$region[j],
                                            site_df$COMID[j], streamcat_sets[i]))
        if('try-error' %in% class(row_base) || nrow(row_base) > 1){
          initerr = TRUE
          row_base = data.frame(COMID=site_df$COMID[j])
        } else {
          initerr = FALSE
        }
      } else {
        row_ext = try(streamcat_from_comid(site_df$region[j],
                                           site_df$COMID[j], streamcat_sets[i]))
        if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
          row_base = left_join(row_base, row_ext)
        }
      }
      
    }
    
    if(nrow(row_base) > 1){
      row_base = data.frame(COMID=site_df$COMID[j])
    }
    
    streamcat_data = rbind.fill(streamcat_data, row_base)
  }
  
  return(streamcat_data)
}


# 2. get NHDPlusV2 data ####

#COMID is the NHD identifier for any reach in the continental U.S.
#add COMIDs to your site table. If this doesn't work, update nhdplusTools
sites$COMID2 = unlist(mapply(comid_from_point, sites$lat,
                            sites$lon, NAD83))
# compare sites with mismatched COMIDS - it looks like the ones from my database are correct.
mismatch <- sites %>% 
  select(site_name, long_name, lat, lon, COMID, COMID2) %>% 
  mutate(diff = COMID - COMID2) %>% filter( diff != 0)

# add state abbreviations for streamcat:
states <- data.frame(US_state = state.name, region = state.abb)
sites <- mutate(sites, COMID = case_when(!is.na(COMID) ~ COMID,
                                         TRUE ~ COMID2)) %>%  
  select(-COMID2) %>%
  filter(!is.na(COMID)) %>%
  left_join(states, by = 'US_state')

sites %>% filter(is.na(region)) %>%
  select(long_name, US_state, region)

sites$region[is.na(sites$region)] <- c('MD', rep('AL', 7), 
                                       'IL', 'TX', 'OR', 'PR', 'PR')
#VPU == NHD vector processing unit. NHDPlusV2 data are downloaded per VPU.
#add VPUs to your site table and determine reach proportions.
a <- mapply(vpu_from_point, sites$lat, sites$lon, NAD83)
a[[409]] <- a[[409]][1] # one of these is giving two vpu's for some reason

sites$VPU = unlist(a, recursive = FALSE)

# reach_props = unlist(mapply(calc_reach_prop, sites$VPU, sites$COMID,
#                                        sites$lat, sites$lon, NAD83, quiet=TRUE))
sites$reach_proportion <- NA
for(i in 1:nrow(sites)){
  x <- try(calc_reach_prop(sites$VPU[i], sites$COMID[i], sites$lat[i], 
                           sites$lon[i], NAD83, quiet = TRUE))
  if(inherits(x, 'try-error')) x <- NA
  if(length(x)==0) x <- NA
  sites$reach_proportion[i] <- x
}
# if the reach proportion cannot be calculated, assume the entire reach is included

sites$reach_proportion[is.na(sites$reach_proportion)] <- 1
#construct list of DSN=component pairs to acquire. see NHDPlus docs for more.
setlist = list('NHDPlusAttributes'='PlusFlowlineVAA',
               'NHDPlusAttributes'='ElevSlope')
#'NHDSnapshot'='NHDFlowline', 'NHDPlusAttributes'='PlusFlowAR'

#retrieve NHDPlusV2 data
nhdplusv2_data = nhdplusv2_bulk(sites, setlist, quiet=TRUE)

#nhd variable names do not have consistent naming conventions. sometimes they're
#all caps; other times camel case. here's a crude way to deal with that.
colnames(nhdplusv2_data) = toupper(colnames(nhdplusv2_data))
nhdplusv2_data = nhdplusv2_data[, ! duplicated(colnames(nhdplusv2_data))]

#pick out the variables you want, then join them to your site data
# nhdplusv2_data = select(nhdplusv2_data, COMID, STREAMORDE, FROMMEAS, TOMEAS, SLOPE,
#                         REACHCODE, AREASQKM, TOTDASQKM, MAXELEVSMO, MINELEVSMO)
sites_nhd = sites %>% 
  select(site_name, lat, lon, coord_datum, region, COMID, VPU, reach_proportion) %>%
  left_join(nhdplusv2_data, by='COMID')
sites_nhd = sites_nhd[! duplicated(sites_nhd$site_name),]

#correct catchment area (AREASQKM) based on where each site falls within its reach.
#use this to correct watershed area (TOTDASQKM) and to determine an areal
#correction factor that can be multiplied with any areal summary data.
sites_nhd$AREASQKM_corr = round(sites_nhd$AREASQKM * sites_nhd$reach_proportion, 5)
sites_nhd$TOTDASQKM_corr = sites_nhd$TOTDASQKM - (sites_nhd$AREASQKM - sites_nhd$AREASQKM_corr)
sites_nhd$areal_corr_factor = sites_nhd$TOTDASQKM_corr / sites_nhd$TOTDASQKM

write_csv(sites_nhd, 'nhd_site_dat.csv')
sites_nhd <- read_csv( 'nhd_site_dat.csv')
# 3. get StreamCat data ####

#find out which streamcat datasets are available (case insensitive)
query_streamcat_datasets()
query_streamcat_datasets('ripbuf')

#construct vector of streamcat datasets to acquire (check variable list for deets)
setlist2 = c('Elevation', 'PRISM_1981_2010', 'NLCD2016', 'Runoff', 
             'ImperviousSurfaces2011', 'Dams', 'USCensus2010', 'EPA_FRS', 
             'Lithology', 'RoadDensity', 'RoadStreamCrossings', 'NABD',
             'AgriculturalNitrogen', 'STATSGO_Set2', 'NADP', 'GeoChemPhys1', 
             'GeoChemPhys2', 'GeoChemPhys3', 'BFI')

#save in chunks, this is a long process and a lot of data.
streamcat_data1 = streamcat_bulk(sites_nhd[1:100,], setlist2)
write_csv(streamcat_data1, 'streamcat_data1.csv')
streamcat_data1 = streamcat_bulk(sites_nhd[101:200,], setlist2)
write_csv(streamcat_data1, 'streamcat_data2.csv')
streamcat_data1 = streamcat_bulk(sites_nhd[201:300,], setlist2)
write_csv(streamcat_data1, 'streamcat_data3.csv')
streamcat_data1 = streamcat_bulk(sites_nhd[301:400,], setlist2)
write_csv(streamcat_data1, 'streamcat_data4.csv')
streamcat_data1 = streamcat_bulk(sites_nhd[401:500,], setlist2)
write_csv(streamcat_data1, 'streamcat_data5.csv')
streamcat_data1 = streamcat_bulk(sites_nhd[501:600,], setlist2)
write_csv(streamcat_data1, 'streamcat_data6.csv')

# compile chunks
streamcat_data <- data.frame()
for(i in 1:6){
  sc <- read_csv(paste0('streamcat_data', i, '.csv'))
  streamcat_data <- bind_rows(streamcat_data, sc)
}


#pick out the variables you want, then join them to your site data
streamcat_data = streamcat_data %>%
  # select(COMID, ElevWs, Precip8110Ws, Tmin8110Ws, Tmax8110Ws, Tmean8110Ws, 
  #        RunoffWs, matches('^Pct[a-zA-z]+2011Ws$'), PermWs, RckDepWs, OmWs, 
  #        WtDepWs, matches('^[a-zA-z0-9]_2008Ws$'), BFIWs, NWs, Al2O3Ws, CaOWs, 
  #        Fe2O3Ws, K2OWs, MgOWs, Na2OWs, P2O5Ws, SWs, SiO2Ws) %>%
  mutate(precip_runoff_ratio=PrecipWs / RunoffWs)

write_csv(streamcat_data, 'streamcat_data.csv')
sites_nhd = left_join(sites_nhd, streamcat_data, by='COMID')
sites_nhd = sites_nhd[! duplicated(sites_nhd$site_name),]

#save yer data
sites_nhd = arrange(sites_nhd, region, site_name)
write.csv(sites_nhd, 'watershed_summary_data.csv', row.names=FALSE)


# 4. get MODIS data (this section incomplete) ####
# VNP13A1
mt_bands("MOD13Q1")
subset1 = mt_subset(product = "MOD13Q1",
                    lat = 40,
                    lon = -110,
                    band = "250m_16_days_NDVI",
                    start = "2004-01-01",
                    end = "2004-02-01",
                    km_lr = 10,
                    km_ab = 10,
                    site_name = "testsite",
                    internal = TRUE,
                    progress = FALSE)

dfx = data.frame("site_name" = paste("test",1:2))
dfx$lat = 40
dfx$lon = -110

# test batch download
subsets = mt_batch_subset(dfx = dfx,
                          product = "MOD11A2",
                          band = "LST_Day_1km",
                          internal = TRUE,
                          start = "2004-01-01",
                          end = "2004-02-01")


