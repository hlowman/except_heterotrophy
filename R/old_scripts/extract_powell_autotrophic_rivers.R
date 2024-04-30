# Explore autotrophy in the cleaned and gap filled dataset from P Savoy 

# A Carter
# 2021-09

setwd('C:/Users/Alice Carter/git/except_heterotrophy/')
library(tidyverse)
library(lubridate)


# Filter Powell Center data set ####
dat <- read.table('data_ignored/daily_predictions.tsv', header = TRUE)

# Remove negative GPP, positive ER, and Rhat > 1.05
filtered <- dat %>%
  mutate(date = as.Date(date),
         year = year(date)) %>%
  filter(GPP >= 0, ER <= 0, 
         GPP.Rhat < 1.05, ER.Rhat < 1.05, K600.Rhat < 1.05)# %>%
  # mutate(GPP = case_when(GPP < 0 ~ 0, TRUE ~ GPP),
  #        ER = case_when(ER > 0 ~ 0, TRUE ~ ER)) 


# Compare filtered data set to the data currently being used in the SP paper 
#  (Bernhardt Savoy et al paper, data from mvlah 10/2021)
# 
sp_dat <- read_csv('../loticlentic_synthesis/data/phil_powell_data/lotic_gap_filled_unlisted.csv') %>%
  rename(site_name = sitecode, date = Date, year = Year)

c_years <- filtered %>%
  group_by(site_name, year) %>%
  filter(sum(!is.na(GPP))/365 >= 0.6,
         sum(!is.na(ER))/365 >= 0.6) %>%
  ungroup()

sp_only <- anti_join(sp_dat, c_years, by = c('site_name', 'date')) %>%
  filter(!is.na(GPP), !is.na(ER), grepl('nwis', site_name))
pw_only <- anti_join(c_years, sp_dat, by = c('site_name', 'date')) %>%
  filter(!is.na(GPP), !is.na(ER))

summary(pw_only)
# compare filtered and unfiltered data sets
png('figures/filtered_data_comparison.png', height = 480, width = 640)
  plot(dat$ER, dat$GPP, xlab = 'ER', ylab = 'GPP')
  points(filtered$ER, filtered$GPP, col = 2, pch = 20)
  points(pw_only$ER, pw_only$GPP, pch = 20, col = 4)
  legend('topright', 
         c('powell center', 'filtered data', 'data not included in Savoy data'),
         pch = c(1, 20, 20), col = c(1,2,4), bty = 'n')
dev.off()
nrow(dat)
nrow(filtered)
nrow(pw_only)

# it looks like the savoy dataset excludes an additional ~10% of the data, but
# they don't look like an autotrophically interesting subset, so likely it 
# wouldn't change our analysis.
# However, I do think it will still be useful to retain the info from the 
# gap filled savoy data set for calculating annual rates where possible. 
# To that end, I've created a summary dataframe below and added it to the 
# working data folder.
# summarize annual values ####
# 
# annual <- sp_dat %>%
#   group_by(site_name, year) %>%
#   summarize(disch_cv = sd(Disch_filled, na.rm = T)/mean(Disch_filled, na.rm = T),
#             across(c('GPP_filled', 'ER_filled', 'NEP_filled'), sum),
#             across(c('Wtemp_filled', 'Disch_filled', 'PAR_filled'),
#                    mean, na.rm = T),
#             gpp_er_cor = cor(GPP, -ER, use = 'complete.obs')) %>%
#   rename_with(function(x) sub('_filled', '', x), everything()) %>%
#   rename(watertemp_mean = Wtemp, disch_mean = Disch, PAR_mean = PAR) %>%
#   ungroup()
# 
# write_csv(annual, 'data_working/annual_gapfilled_values.csv')
annual <- read_csv('data_working/annual_gapfilled_values.csv')

# using only the years with at least 60% of data, which sites are autotrophic 
# at the annual timescale? 

auto_siteyears <- annual %>%
  filter(NEP > 0) %>%
  rename_with(.fn = function(x) paste0(x, '_annual'), 
              .cols = c(-site_name, -year)) %>%
  mutate(siteyear = paste(site_name, year, sep = '_'))

auto_sites <- c_years %>% 
  filter(paste(site_name, year, sep = '_') %in% auto_siteyears$siteyear) %>%
  left_join(auto_siteyears, by = c('site_name', 'year')) %>%
  select(-siteyear)

write_csv(auto_sites, 'data_working/autotrophic_siteyears_daily.csv')
write_csv(auto_siteyears, 'data_working/autotrophic_siteyears_annual.csv')

