# batch summarizing of watershed data for the continental USA
# Modified script from Mike Vlah (vlahm13@gmail.com)

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

#   site data must include latitude and longitude columns (decimal degrees)
# dataframe generated from filter_powell_estimates scipt
sites = read_csv('data_working/Autotrophy_lit_table.csv')
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
  
  ftpdir = paste0('https://gaftp.epa.gov/epadatacommons/ORD/',
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
          if(row_ext[['WsPctFull']] < 90){
            readr::write_lines(paste(j, site_df$COMID[j], streamcat_sets[i], row_ext$WsPctFull),
                               file = 'data_ignored/streamcat/low_coverage_data.txt',
                               append = TRUE)
          }
          newcols <- setdiff(colnames(row_ext), colnames(row_base))
          row_ext <- select(row_ext, COMID, any_of(newcols))
          row_base = left_join(row_base, row_ext, by = 'COMID')
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
# only get COMIDS for the sites with missing ones because I have double checked all of the powell center ones
sites$COMID = unlist(mapply(comid_from_point, sites$Latitude,
                            sites$Longitude, NAD83))


# add state abbreviations for streamcat:
sites <- dplyr::rename(sites, region = US_State)

#VPU == NHD vector processing unit. NHDPlusV2 data are downloaded per VPU.
#add VPUs to your site table and determine reach proportions.
a <- mapply(vpu_from_point, sites$Latitude, sites$Longitude, NAD83)
a[sapply(a, function(x) !length(x))] <- NA
sites$VPU = unlist(a, recursive = FALSE)

# reach_props = unlist(mapply(calc_reach_prop, sites$VPU, sites$COMID,
#                                        sites$lat, sites$lon, NAD83, quiet=TRUE))
sites$reach_proportion <- NA
for(i in 1:nrow(sites)){
  x <- try(calc_reach_prop(sites$VPU[i], sites$COMID[i], sites$Latitude[i], 
                           sites$Longitude[i], NAD83, quiet = TRUE))
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
nhdplusv2_data = nhdplusv2_bulk(sites[!is.na(sites$COMID),], setlist, quiet=TRUE)

#nhd variable names do not have consistent naming conventions. sometimes they're
#all caps; other times camel case. here's a crude way to deal with that.
colnames(nhdplusv2_data) = toupper(colnames(nhdplusv2_data))
nhdplusv2_data = nhdplusv2_data[, ! duplicated(colnames(nhdplusv2_data))]

#pick out the variables you want, then join them to your site data
nhdplusv2_data = select(nhdplusv2_data, COMID, STREAMORDE, HYDROSEQ,  
                        STARTFLAG, REACHCODE, AREASQKM, TOTDASQKM, LENGTHKM, 
                        SLOPE, SLOPELENKM, TIDAL, MAXELEVSMO, MINELEVSMO)

colnames(nhdplusv2_data) = paste0("NHD_", colnames(nhdplusv2_data))
nhdplusv2_data <- dplyr::rename(nhdplusv2_data, COMID = NHD_COMID)
nhdplusv2_data <- nhdplusv2_data[!duplicated(nhdplusv2_data$COMID),]

sites_nhd = sites %>%
  left_join(nhdplusv2_data, by = 'COMID') 

#correct catchment area (AREASQKM) based on where each site falls within its reach.
#use this to correct watershed area (TOTDASQKM) and to determine an areal
#correction factor that can be multiplied with any areal summary data.
sites_nhd$NHD_AREASQKM_corr = round(sites_nhd$NHD_AREASQKM * sites_nhd$reach_proportion, 5)
sites_nhd$NHD_TOTDASQKM_corr = sites_nhd$NHD_TOTDASQKM - 
  (sites_nhd$NHD_AREASQKM - sites_nhd$NHD_AREASQKM_corr)
sites_nhd$areal_corr_factor = sites_nhd$NHD_TOTDASQKM_corr / sites_nhd$NHD_TOTDASQKM

# 3. get StreamCat data ####

#construct vector of streamcat datasets to acquire (check variable list for deets)

setlist2 = c('Elevation', 'PRISM_1981_2010', 'NLCD2011', 'NLCD2016', 'Runoff', 
             'ImperviousSurfaces2011', 'Dams', 'USCensus2010', 'EPA_FRS', #should this be just ImperviousSurfaces?
             'Lithology', 'RoadDensity', 'RoadStreamCrossings', 'NABD',
             'AgriculturalNitrogen', 'STATSGO_Set2', 'NADP', 'GeoChemPhys1',
             'GeoChemPhys2', 'GeoChemPhys3', 'BFI')

streamcat_data = streamcat_bulk(sites_nhd[!is.na(sites_nhd$COMID),], setlist2)

#pick out the variables you want, then join them to your site data
streamcat_filtered = streamcat_data %>%
  select(COMID, ElevWs, PrecipWs, TminWs, TmaxWs, TmeanWs, RunoffWs, 
         matches('[0-9]{4}Ws$'), ends_with(c('StorWs', 'DensWs')),
         matches('^Pct[a-zA-z]+Ws$'), CBNFWs, FertWs, ManureWs, WtDepWs, OmWs,
         PermWs, RckDepWs, BFIWs, HydrlCondWs, NWs, Al2O3Ws, CaOWs, Fe2O3Ws, 
         K2OWs, MgOWs, Na2OWs, P2O5Ws, SWs, SiO2Ws) %>%
  mutate(precip_runoff_ratio=PrecipWs / RunoffWs)

sites_nhd = sites_nhd %>%
  select(-ends_with(c('Cat', '2011Ws')), 
         -starts_with(c('NHD_CatPctFull', 'NHD_WsPctFull'))) %>%
  left_join(streamcat_filtered, by='COMID')
sites_nhd = sites_nhd[! duplicated(sites_nhd$River),]

#save yer data
write_csv(sites_nhd, 'data_working/literature_streams_watershed_summary_data_0.csv')
