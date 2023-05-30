# add baseflow to the daily datasets and calculate max time period between floods
# for each site year based on different flood thresholds.

library(tidyverse)
library(dataRetrieval)
library(moments)
# install.packages('hydrostats')
# library(hydrostats)
# the package EcoHydRology is no longer available on Cran, but this is an old version
# https://cran.r-project.org/src/contrib/Archive/EcoHydRology/
# install.packages('C:/Users/alice.carter/Downloads/EcoHydRology_0.4.12.1.tar.gz', 
#                  repos = NULL, type ='source', dependencies = TRUE)

# library(EcoHydRology)
# install.packages("FlowScreen")
library(FlowScreen)

sites <- read_csv('data_working/Autotrophy_lit_table.csv')

# Download USGS gage data to calculate discharge metrics:
discharge <- dataRetrieval::readNWISdv(
  siteNumbers = sites$nwis_gage[!is.na(sites$nwis_gage)],
  parameterCd = '00060',
  statCd = '00003')

group_by(discharge, site_no) %>%
  summarize(n = n()/365)
# after some exploration and reading about baseflow separation methods, 
# It seems that the Eckhardt method is recommended, so we'll use that.
add_baseflow_column <- function(dat){
  dat$bf <- FlowScreen::bf_eckhardt(dat$discharge_m3s, a = 0.975, 
                                    BFI = 0.8)
  dat$bfi <- dat$bf/dat$discharge_m3s
  plot(dat$Date, dat$discharge_m3s, type = 'l')
  lines(dat$Date, dat$bf, col = 2)
  
  return(dat)
}

discharge <- discharge %>%
  left_join(select(sites, River, site_no = nwis_gage)) %>%
  mutate(year = lubridate::year(Date),
         discharge_m3s = X_00060_00003 * (0.3048^3)) %>%
  select(-starts_with('X')) 

discharge <- split(discharge, discharge$River)
# names(discharge)
# add_baseflow_column(discharge[[9]])
# plot(discharge[[9]]$Date, discharge[[9]]$discharge_m3s, type = 'l')
# discharge[[9]]$Date[is.na(discharge[[9]]$discharge_m3s)]
## Big springs, Don't include the random data from 1923
discharge[[1]] <- filter(discharge[[1]], year > 1990)

## Colorado river, only use data from after the dam was installed
discharge[[3]] <- filter(discharge[[3]], year > 2000)

## Kinnickinnic river, remove data after 2021 due to data gap
discharge[[5]] <- filter(discharge[[5]], year < 2022)

## Klamath remove data before gap
discharge[[6]] <- filter(discharge[[6]], year > 1950)

## SF Humboldt remove data before gap and interpolate remaining missing day
discharge[[9]] <- filter(discharge[[9]], year < 2023)

discharge <- lapply(discharge, add_baseflow_column)
q_siteyears <- lapply(discharge, function(x) group_split(x, year))

q_siteyears_split <- unlist(q_siteyears, recursive = F, use.names = T)
names(q_siteyears_split) <- sapply(q_siteyears_split, 
                                   function(x) paste(x$site_no[1], 
                                                     x$year[1], sep = "_"))
q_siteyears_split <- Filter(function(x) nrow(x)>=274,
                            q_siteyears_split)

q_siteyears_split <- lapply(q_siteyears_split, as.data.frame)

calc_storm_gaps <- function(dat, threshold = 0.8){
  dat$storm = case_when(dat$bfi < threshold ~ 1,
                        TRUE ~ 0)    
  storms <- rle(dat$storm)
  storms <- data.frame(vals = storms$values,
                       lengths = storms$lengths)
  N = nrow(storms)
  
  # count interstorm intervals that wrap around from Dec -> Jan
  # if(storms$vals[1] == 0 & storms$vals[N] == 0){
  #    storms$lengths[1] = storms$lengths[1] + storms$lengths[N]
  #    storms <- storms[1:(N-1),]
  # }
  storms <- storms[2:N,]
  storms <- dplyr::filter(storms, vals == 0) 
  s <- data.frame(bfi_threshold = threshold,
                  med_interstorm = median(storms$lengths),
                  max_interstorm = max(storms$lengths))
  
  return(s)
  
}

plot_bf_separation<- function(dat){
  par(mfrow = c(2,1),
      mar = c(2,4,1,2), 
      oma = c(3,1,1,0))
  plot(dat$Date, dat$discharge_m3s, type = 'l', ylab = 'discharge', 
       main = dat$River[1])
  lines(dat$Date, dat$bf, col = 2)
  
  bf <- data.frame()
  for(threshold in seq(0.5, 0.9, by = 0.01)){
    d <- calc_storm_gaps(dat, threshold)
    bf <- bind_rows(bf, d)
  }
  
  plot(bf$bfi_threshold, bf$max_interstorm, type = 'l',
       ylim = range(c(bf$max_interstorm, bf$med_interstorm)),
       ylab = 'max interstorm days')
  lines(bf$bfi_threshold, bf$med_interstorm, lty = 2)
  mtext('baseflow index storm threshold', 1, 2.5)
}    

calc_RBI <- function(dat){
  q <- dat %>% select(Date, discharge_m3s)
  if(sum(is.na(q$discharge_m3s))/nrow(q) > 0.75){
    print('need at least 75% of days to have discharge to proceed')
    return(NA_real_)
  }
  q <- mutate(q, discharge_m3s = zoo::na.approx(discharge_m3s, na.rm = F)) %>%
    filter(!is.na(discharge_m3s))
  d <- diff(q$discharge_m3s)
  RBI <- sum(abs(d))/sum(q$discharge_m3s[2:nrow(q)])
  
  return(RBI)
  
}

# for(i in 1:length(lotic_siteyears_split)){
#   plot_bf_separation(lotic_siteyears_split[[923]])
# }

# based on looking through many of these plots, 0.75 seems like a reasonable
# BFI cutoff for what counts as a 'storm'

# calc storm intervals:
bf <- data.frame()
for(i in 1:length(q_siteyears_split)){
  dat <- q_siteyears_split[[i]]
  d <- calc_storm_gaps(dat, threshold = 0.75)
  d$RBI <- calc_RBI(dat)
  d$Disch_mean <- mean(dat$discharge_m3s)
  d$Disch_cv <- sd(dat$discharge_m3s)/mean(dat$discharge_m3s)
  d$site_no <- dat$site_no[1]
  d$Year = dat$year[1]
  
  bf <- bind_rows(bf, d)
}

bf <- relocate(bf, c(site_no, Year))

bf_sum <- bf %>% group_by(site_no) %>%
  summarize(across(.cols = any_of(c('med_interstorm', 'max_interstorm', 'RBI', 
                                    'Disch_mean', 'Disch_cv')), 
                   .fns = list(mean = ~mean(.), sd = ~sd(.), median = ~median(.)))) 
  


ggplot(bf, aes(Year, RBI)) + geom_line() +
  facet_wrap(.~site_no, scales = 'free_x')

sites_nhd <- read_csv('data_working/literature_streams_watershed_summary_data.csv')
# add width based on google earth
sites_nhd$width_m <- c(35, 5, NA, 40, 30, NA, 20, 110, 5, 20, NA, 5, 12, 35)
# add publication year
sites_nhd$year <- c(2016, 1972, 2011, 1957, 2003, 2003, 2012, 2017, 2008, 2020, 2017, 2017, 2020, 2008)
# sites_nhd <- select(sites_nhd, -max_interstorm, -RBI, -LAI, -PAR_inc, -PAR_surface)%>% 
sites_nhd <- sites_nhd %>%
          left_join(select(bf_sum, 
                 nwis_gage = site_no,
                 max_interstorm = max_interstorm_median,
                 RBI = RBI_median, 
                 Disch_mean = Disch_mean_median,
                 Disch_cv = Disch_cv_median),
          by = 'nwis_gage')

# add light data:
light <- read_csv('data_ignored/light/daily_modeled_light_lit_sites.csv') %>%
  mutate(year = lubridate::year(date),
         siteyear = paste(site, year, sep = '_')) 
full_ys <- light %>% 
  group_by(siteyear) %>%
  summarize(n = n()) %>%
  filter(n == 365)

light_years <- light %>% 
  filter(siteyear %in% full_ys$siteyear) %>%
  group_by(site, year) %>%
  summarize(LAI_mean = mean(LAI),
            PAR_mean = mean(PAR_inc),
            PAR_sum = sum(PAR_inc),
            Stream_PAR_sum = sum(PAR_surface),
            PAR_kurt = kurtosis(PAR_inc))
light <- light_years %>%
  group_by(site) %>%
  summarize(LAI = mean(LAI_mean),
            PAR_mean = mean(PAR_mean),
            PAR_sum = mean(PAR_sum),
            Stream_PAR_sum = mean(Stream_PAR_sum),
            PAR_kurt = mean(PAR_kurt)) %>%
  mutate(nwis_gage = substr(site, 6, 13)) %>%
  select(-site)
  
# add MOdis Annual NPP
NPP <- read_csv('data_ignored/light/MODIS/NPP/autotrophy-NPP-MOD17A3HGF-006-results.csv')
NPP <- NPP %>% select(ID, Date, NPP = MOD17A3HGF_006_Npp_500m,
               NPP_QC = MOD17A3HGF_006_Npp_QC_500m) %>%
    # filter(NPP_QC == 250) # this indicates the site is urban and NPP can't be estimated
    #                       # for our purposes, we will set it to zero. This is a highly
    #                       # developed site, so it is not a bad assumption.
    mutate(NPP = case_when(NPP_QC == 250 ~ 0,
                           TRUE ~ NPP)) %>%
    group_by(ID) %>%
    summarize(MOD_ann_NPP = mean(NPP)) %>%
    mutate(nwis_gage = substr(ID, 5, 12))

sites_nhd <- left_join(sites_nhd, NPP, by = 'nwis_gage') %>%
    select(-ID)

sites_nhd <- left_join(sites_nhd, light, by = 'nwis_gage')

write_csv(sites_nhd, 'data_working/literature_streams_watershed_summary_data.csv')
