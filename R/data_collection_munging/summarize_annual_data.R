# Build dataframe with individual site year summaries paired with 
# watershed data for all good coverage siteyears

# package from Bernhardt et al 2022, 
# install by downloading zip file from figshare datarelease
install.packages('C:/Users/alice.carter/Downloads/BernhardtMetabolism_0.1.1.zip',
                 repos=NULL, type='source')

library(BernhardtMetabolism) 
library(stringr)
library(tidyverse)
# Summarize data from gapfilled years that have been filtered to have only 
# data with at least 60% coverage and uncorrelated KxER
lotic_gap_filled <- readRDS('data_ignored/Bernhardt_2022/lotic_gap_filled.rds')

lotic_siteyears <- lapply(lotic_gap_filled, function(x) group_split(x, Year))

lotic_siteyears_split <- unlist(lotic_siteyears, recursive = F, use.names = T)
names(lotic_siteyears_split) <- sapply(lotic_siteyears_split, 
         function(x) paste(x$Site_ID[1], x$Year[1], sep = "_"))

lotic_siteyears_split <- lapply(lotic_siteyears_split, as.data.frame)

metrics_compiled <- data.frame()


for(i in 1:length(lotic_siteyears_split)){
  newrow <- BernhardtMetabolism::calc_site_metrics(names(lotic_siteyears_split)[i], 
                                         lotic_siteyears_split)
  newrow$year <-
    str_match(names(lotic_siteyears_split)[i], '[A-Za-z0-9_]+?([0-9]{4})$')[2]

  metrics_compiled <- bind_rows(metrics_compiled, newrow)
}

write_csv(metrics_compiled, 
          'data_ignored/Bernhardt_2022/lotic_gap_filled_annual_summaries.csv')

metrics_compiled <- read_csv('data_ignored/Bernhardt_2022/lotic_gap_filled_annual_summaries.csv')
# Combine with Watershed Characteristics:
ws <- read_csv('data_356rivers/watershed_summary_data.csv')

sum_columns <- intersect(colnames(metrics_compiled), colnames(ws))
ws <- select(ws, -any_of(sum_columns))
yearly <- metrics_compiled %>%
  select(site_name = Site_ID, year, any_of(sum_columns)) %>%
  left_join(ws, by = 'site_name') %>%
  relocate(c('long_name', 'Source', 'lat', 'lon', 'coord_datum', 'alt', 
             'alt_datum', 'site_type', 'COMID', 'VPU'), .after = 'year') %>%
  select(-NHD_QE_MA, -NHD_VE_MA, -NHD_GNIS_NAME, -NHD_FTYPE, -NHD_FCODE) 


write_csv(yearly, 'data_working/annual_summary_data.csv')
