# Compare/merge filtered Powell Center dataset to the 
# StreamPULSE dataset from Bernhardt et al 2022
# Add covariates from GRDO

# A Carter
# 3/7/2022

# setup
library(tidyverse)
library(purrr)

# Read in and filter Powell data: ####
pw_site_dat <- read.table('data_356rivers/site_data.tsv', sep = '\t', header = T)
pw <- read.table('data_ignored/daily_predictions.tsv', sep = '\t', header = T)

# Subset to high quality days, with Rhat < 1.05 and plausible GPP and ER values

sum(pw$GPP.Rhat >= 1.05 | pw$ER.Rhat >= 1.05 | pw$K600.Rhat >= 1.05)/nrow(pw) 
# poor convergence: 4.2% of days
sum(pw$GPP < 0 & pw$GPP.Rhat < 1.05)/sum(!is.na(pw$GPP)) # negative GPP 10% of days
sum(pw$ER > 0 & pw$ER.Rhat < 1.05)/sum(!is.na(pw$ER))  # positive ER 8.7% of days
# total lost days ~ 14%

HQdays <-  pw %>%
  rename(GPP_raw = GPP, ER_raw = ER) %>%
  mutate(GPP = case_when(GPP.Rhat >= 1.05 ~ NA_real_,
                         K600.Rhat >= 1.05 ~ NA_real_,
                         ER.Rhat >= 1.05 ~ NA_real_,
                         # GPP.upper < 0 ~ NA_real_,
                         GPP_raw < 0 ~ NA_real_,
                         # GPP_raw < -1 ~ NA_real_,
                         # ER_raw > 1 ~ NA_real_,
                         # GPP_raw < 0 ~ 0,
                         TRUE ~ GPP_raw),
         ER = case_when(ER.Rhat >= 1.05 ~ NA_real_,
                        K600.Rhat >= 1.05 ~ NA_real_,
                        GPP.Rhat >= 1.05 ~ NA_real_,
                        # ER.lower > 0 ~ NA_real_,
                        ER_raw > 0 ~ NA_real_,
                        # ER_raw > 1 ~ NA_real_,
                        # GPP_raw < -1 ~ NA_real_,
                        # ER_raw > 0 ~ 0,
                        TRUE ~ ER_raw),
         K600 = case_when(GPP.Rhat >= 1.05 ~ NA_real_,
                          ER.Rhat >= 1.05 ~ NA_real_,
                          K600.Rhat >= 1.05 ~ NA_real_),
         NEP = GPP + ER)

# png('figures/powell_dataset_filtering.png' )
#   plot(HQdays$GPP_raw, HQdays$ER_raw, pch = 20, col = 'grey')
#   points(HQdays$GPP, HQdays$ER, pch = 20)
#   points(HQdays$GPP_raw[!is.na(HQdays$GPP) & (HQdays$GPP_raw < 0)],
#          HQdays$ER_raw[!is.na(HQdays$GPP) & (HQdays$GPP_raw < 0)], pch = 20, col = 'brown3')
#   points(HQdays$GPP_raw[!is.na(HQdays$ER) & (HQdays$ER_raw > 0)],
#          HQdays$ER_raw[!is.na(HQdays$ER) & (HQdays$ER_raw > 0)], pch = 20, col = 'brown3')
#   legend('topright', c('raw', 'filtered (88% of days)', 'set to zero (2-6% of days)'), 
#          col = c('grey', 'black', 'brown3'), pch = 20, bty = 'n')
# dev.off()

sum(!is.na(HQdays$GPP))/sum(!is.na(HQdays$GPP_raw)) # 86% of days kept
sum(!is.na(HQdays$ER))/sum(!is.na(HQdays$ER_raw))   # 87% of days kept

# Comparison: if we set estimates that contain 0 in the CI to zero, we would 
# keep an additional:
nrow(HQdays[(HQdays$GPP.upper > 0) & (HQdays$GPP_raw < 0),])/
  sum(!is.na(HQdays$GPP_raw)) # 8.1% of days
nrow(HQdays[(HQdays$ER.lower < 0) & (HQdays$ER_raw > 0),])/
  sum(!is.na(HQdays$ER_raw))  # 3.6% of days 


# Data from Bernhardt et al 2022: ####
# https://figshare.com/articles/software/Code_and_RDS_data_for_Bernhardt_et_al_2022_PNAS_/19074140?backTo=/collections/Data_and_code_for_Bernhardt_et_al_2022_PNAS_/5812160
# Note, to get the lotic_standardized_full dataset 
sp <- readRDS('data_ignored/Bernhardt_2022/lotic_standardized_full.rds')
sp <- map_dfr(sp, bind_rows) 

# This dataset contains gapfilled data (following Bernhardt et al 2022) for all
# site years with at least 60% coverage
sp_filled <- readRDS('data_ignored/Bernhardt_2022/lotic_gap_filled.rds')
sp_filled <- map_dfr(sp_filled, bind_rows) %>%
  select(Site_ID, Date, ends_with('_filled')) %>%
  select(-ends_with('_C_filled'))
rownames(sp_filled) <- NULL

sp_sub <- left_join(sp, sp_filled, by = c('Site_ID', 'Date'))  %>%
  select(-starts_with('DO_'), -ends_with('raw'), -GPP, -ER, 
         -U_ID, -temp_water, -discharge) %>%
  rename(site_name = Site_ID, date = Date)

# subset columns:
pw_hq <- HQdays %>%
  select(-resolution, -ends_with('n_eff'), -day.length, -shortwave) %>%
  mutate(date = as.Date(date))

# Powell center estimates complete with bernhardt covatiates
sp_pw <- left_join(pw_hq, sp_sub, by = c('site_name', 'date'))
summary(sp_pw)

# Save joined database
saveRDS(sp_pw, 'data_356rivers/high_quality_daily_metabolism_with_SP_covariates.rds')
# zip('data_356rivers/high_quality_daily_metabolism_with_covariates.rds.zip', 
#     'data_ignored/high_quality_daily_metabolism_with_covariates.rds', 
#     flags = '-rj9X')


# Combine annual dataframes ####
# Dataframe containing summary from all site years that have:
# - minimum 60% coverage after filtering 
# - ER x K600 correlation < 0.6
# - maximum K600 < 100
sp_dat <- readRDS('data_ignored/Bernhardt_2022/lotic_site_info_filtered.rds')
GRDO <- read_csv('data_ignored/GRDO_GEE_HA_NHD.csv', guess_max = 10000)

# remove sites not in the powell center and select desired columns:
sp_dat <- sp_dat %>% tibble() %>%
  # filter(Source == 'USGS (Powell Center)') %>%
  select(site_name = Site_ID, COMID, VPU, Source, Lat, Lon, ndays, 
         nyears, StreamOrde, Azimuth, tree_height = TH, Width, ends_with('sum'),
         starts_with(c('ann', 'Wtemp', 'Disch', 'PAR', 'LAI')), MOD_ann_NPP)

# Remove sites that aren't streams and check for duplicates:
pw_site_dat %>% filter(site_type != "ST") %>% select(site_name, long_name, site_type, lat, lon)
# all sites are either streams, canals, ditches or tidal streams. For now I will keep them all

pw_site_dat[duplicated(pw_site_dat$site_name)|duplicated(pw_site_dat$site_name, fromLast = T),]
pw_site_dat[duplicated(pw_site_dat$nhdplus_id)|duplicated(pw_site_dat$nhdplus_id, fromLast = T),] %>%
  select(nhdplus_id, site_name, long_name, site_type, lat, lon) %>% filter(!is.na(nhdplus_id))

# a handful of these sites are bottom gages in deep rivers (the Klamath). Remove those.
pw_site_dat <- pw_site_dat %>% filter(!grepl('BOTTOM$', long_name))
duplicates <-
  pw_site_dat[duplicated(pw_site_dat$nhdplus_id)|duplicated(pw_site_dat$nhdplus_id, fromLast = T),] %>%
  select(nhdplus_id, site_name, long_name, site_type, lat, lon) %>% filter(!is.na(nhdplus_id))

# export this file to cross check lat longs with NHD segments manually:
write_csv(duplicates, 'data_ignored/duplicate_COMIDS.csv')

# several sites fall along Bayou Lake and City Park Lagoon in New Orleans,
# These show up as being in the NHD HR but not NHD and several share COMIDS.
# For this project, I think they should be excluded. (Checking above, they also
# contain no days of metabolism)
pw_site_dat <- pw_site_dat %>% 
  filter(!(nhdplus_id %in% c(18910456, 18910424)))

# The remaining duplicates are in fact on the same comid sections.

# Subset GRDO to contain hydroatlas and NHD data
GRDO_sub <- GRDO %>%
  rename(site_name = SiteID, IGBP_LU_category = LU_category) %>%
  select(c(3, 6:8, 13:65, 67:100, 145)) %>% glimpse()

df <- full_join(sp_dat, pw_site_dat, by = 'site_name') %>% 
  mutate(COMID = if_else(is.na(nhdplus_id), COMID, 
                         nhdplus_id),# the comids from the Powell data have been crossref'd with a map. Keep those
         lat = if_else(is.na(lat), Lat, lat),
         lon = if_else(is.na(lon), Lon, lon)) %>% 
  select(-nhdplus_id, -nwis_id, -Lat, -Lon) %>% 
  left_join(GRDO_sub, by = 'site_name') %>%
  relocate(long_name, 46:51, .after = site_name)

saveRDS(df, 'data_ignored/site_data.rds') 
