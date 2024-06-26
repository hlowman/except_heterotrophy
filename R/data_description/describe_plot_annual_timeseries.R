# Patterns in Annual Autotrophy
# Alice Carter et al.

# Describe the patterns in metabolism across all siteyears with at least 60% coverage
# Dataset used is the gap-filled data from the Bernhardt data release

#### Setup ####

# Load packages.
library(tidyverse)
library(ggthemes)

# Load data.
dat <- read_csv('data_working/across_sites_model_data.csv')
sites <- read_tsv("data_356rivers/site_data.tsv")

#### Summaries ####

##### Autotrophic sites #####

# Descriptive Statistics and data distributions
# The dataset here uses gap-filled GPP and ER values, to match what is used in
# quantile regression modeling efforts for annual productivity estimates.
ann <- dat
# 921 site-years
unique(ann$site_name)
# 236 sites total

# Figure of all sites with autotrophic ones colored.
png('figures/distribution_annual_NEP.png', width = 600, height = 400)
    par(mfrow = c(1, 2), mar = c(5,2,3,1))
    plot(density(ann$ann_NEP_C), xlim = c(-6200, 2000), main = '', yaxt = 'n',
         xlab = expression(paste('Annual NEP (g',O[2],'/',m^2,'/d)' )), ylab = '')
    mtext('Density', 2, 0.7)
    dd = density(ann$ann_NEP_C, na.rm = T)
    ddo = order(dd$x)
    xdens = dd$x[ddo]
    ydens = dd$y[ddo]
    xdens_ut = xdens[xdens >= 0]
    ydens_ut = ydens[xdens >= 0]
    polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
            col='lightgreen', border='lightgreen')
    lines(dd)
    
    dd = density(log(-ann$ann_GPP_C/ann$ann_ER_C))
    plot(dd, main = '', xlab = 'Log(GPP/ER)' , ylab = '', yaxt = 'n')
    ddo = order(dd$x)
    xdens = dd$x[ddo]
    ydens = dd$y[ddo]
    xdens_ut = xdens[xdens >= 0]
    ydens_ut = ydens[xdens >= 0]
    polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
            col='lightgreen', border='lightgreen')
    lines(dd)
    par(mfrow = c(1,1), new = T)
    mtext(paste0('Autotrophic site years = ',
                 100 * round(length(which(ann$ann_NEP_C>0))/nrow(ann),3), 
                 '% of ', nrow(ann)), side = 3, line = 1, cex = 1.5)
dev.off()

# How many site-years are autotrophic at an annual scale?
ann_aut_sites <- ann %>% 
  mutate(siteyear = paste(site_name, year, sep = '_')) %>%
  filter(PR > 1) # 87, so ~9% of 921 site-years total

# How many UNIQUE sites?
length(unique(ann_aut_sites$site_name))
# 37 sites are autotrophic at the annual scale

# And what do the summaries of these look like?
# Start with the original, full dataset
ann_aut_sites_summ <- ann %>%
  # filter down to the 37 sites of choice
  filter(site_name %in% unique(ann_aut_sites$site_name)) %>%
  group_by(site_name, lat, lon) %>%
  summarize(med_PR = median(PR),
            max_PR = max(PR)) %>%
  ungroup()

# Percentiles of autotrophy
quantile(ann_aut_sites_summ$med_PR, c(0.10, 0.90))
# 10th = 0.84
# 90th = 1.35

quantile(ann_aut_sites_summ$max_PR, c(0.10, 0.90))
# 10th = 1.04
# 90th = 1.48

# How many sites are autotrophic?
aut_sites <- ann %>% 
  select(site_name, ER = ann_ER_C, GPP = ann_GPP_C, PR, NEP = ann_NEP_C) %>%
  group_by(site_name) %>%
  # Using Alice C's workflow to calculate median PR/NEP across all site-years
  summarize(across(ER:NEP, median, na.rm = T)) %>%
  ungroup() %>%
  filter(PR > 1) # 15, so ~6% of 236 sites total

# make a table of sites that are at least occasionally autotrophic
aut_sites2 <- ann %>% 
  select(site_name, lat, lon, ER = ann_ER_C, GPP = ann_GPP_C, 
         PR, NEP = ann_NEP_C, NLCD_LUCat, IGBP_LU_category, Pct_impcov, PrecipWs,
         Dam_densityperkm2, Dam_normal_vol_m3km2, Waste_point_srcs_perkm2) %>%
  group_by(site_name, NLCD_LUCat, IGBP_LU_category) %>%
  # Using Alice C's workflow to calculate median PR/NEP across all site-years
  summarize(NEP_max = max(NEP, na.rm = T),
            PR_max = max(PR, na.rm = T),
            across(ER:NEP, median, na.rm = T)) %>%
  ungroup() %>%
  arrange(desc(PR)) %>%
  left_join(select(sites, site_name, long_name)) %>%
  filter(PR_max > 1) # 15, so ~6% of 236 sites total

write_csv(aut_sites2, 'data_working/autotrophic_sites.csv')

# Print site names
site_names <- sites %>%
  select(site_name, long_name)

aut_sites <- left_join(aut_sites, site_names)

unique(aut_sites$long_name)

# also printing sites that *aren't* autotrophic overall for categorization in manuscript
always_aut <- unique(aut_sites$site_name)
sometimes_aut <- unique(ann_aut_sites$site_name)
occasionally_aut <- sites %>%
  select(site_name, long_name) %>% # select only site info
  filter(site_name %in% sometimes_aut) %>% # keep sites that have at least 1 year of aut
  filter(!site_name %in% always_aut) # remove sites that are predominantly autotrophic on annual scale

# Need to re-create dataset with site-years of data but only for truly
# autotrophic sites.
aut_site_names <- unique(aut_sites$site_name)

aut_sites_allyrs <- ann %>% 
  filter(site_name %in% aut_site_names)

# GPP Summary:
# minimum GPP
min(aut_sites_allyrs$ann_GPP_C) # 282

# maximum GPP
max(aut_sites_allyrs$ann_GPP_C) # 3,755

# mean GPP
mean(aut_sites_allyrs$ann_GPP_C) # 732

# sd GPP
sd(aut_sites_allyrs$ann_GPP_C) # 712

# NEP Summary:
# minimum NEP
min(aut_sites_allyrs$ann_NEP_C) # -254

# maximum NEP
max(aut_sites_allyrs$ann_NEP_C) # 1,199

# mean NEP
mean(aut_sites_allyrs$ann_NEP_C) # 109

# sd NEP
sd(aut_sites_allyrs$ann_NEP_C) # 202

##### Heterotrophic sites #####

# How many site-years are heterotrophic?
ann_het_sites <- ann %>% 
  mutate(siteyear = paste(site_name, year, sep = '_')) %>%
  filter(PR <= 1 ) # 834

# How many UNIQUE sites?
length(unique(ann_het_sites$site_name))
# 228 sites are heterotrophic at the annual scale (including occasionally het. sites)

# How many sites are heterotrophic?
het_sites <- ann %>% 
  select(site_name, ER = ann_ER_C, GPP = ann_GPP_C, PR, NEP = ann_NEP_C) %>%
  group_by(site_name) %>%
  # Using Alice C's workflow to calculate median PR/NEP across all site-years
  summarize(across(everything(), median, na.rm = T)) %>%
  ungroup() %>%
  filter(PR <= 1) # 221, so ~94% of 236 sites total

# Print site names
het_sites <- left_join(het_sites, site_names)
unique(het_sites$long_name) # some long names missing

# Need to re-create dataset with site-years of data but only for truly
# heterotrophic sites.
het_site_names <- unique(het_sites$site_name)

het_sites_allyrs <- ann %>% 
  filter(site_name %in% het_site_names)

# min GPP at heterotrophic sites
min(het_sites_allyrs$ann_GPP_C) # 10

# max GPP at heterotrophic sites
max(het_sites_allyrs$ann_GPP_C) # 3,650

# mean GPP at heterotrophic sites
mean(het_sites_allyrs$ann_GPP_C) # 264

# sd GPP at heterotrophic sites
sd(het_sites_allyrs$ann_GPP_C) # 314

# NEP Summary:
# minimum NEP
min(het_sites_allyrs$ann_NEP_C) # -1,770

# maximum NEP
max(het_sites_allyrs$ann_NEP_C) # 308

# mean NEP
mean(het_sites_allyrs$ann_NEP_C) # -276

# sd NEP
sd(het_sites_allyrs$ann_NEP_C) # 279

#### Other Figures ####

# Create new dataframe only for autotrophic sites and create new column
# with cumulative GPP on a daily basis.
aut <- dat %>% 
  mutate(siteyear = paste(site_name, year, sep = '_')) %>%
  filter(siteyear %in% aut_sites$siteyear) %>%
  group_by(siteyear) %>%
  arrange(DOY) %>%
  mutate(GPP_cum = cumsum(GPP_filled))

# Plot cumulative GPP curves for autotrophic sites.
# Pecos River is far and away the greatest.
dev.off() 
plot(aut$DOY, aut$GPP_cum, type = 'n',  xlab = 'Day of Year',
     ylab = '', ylim = c(0, 15000))
mtext(expression(paste('GPP (g',O[2],m^{-2},d^{-1},')')), 2, 2.2)
for(i in 1:nrow(aut_sites)){
  d <- aut %>% filter(siteyear == aut_sites$siteyear[i]) %>%
    arrange(DOY)
  lines(d$DOY, d$GPP_cum, col = alpha('forestgreen', 0.5))
  
}

# Distribution of cumulative(?) water temperature, discharge, PAR, GPP, ER
png('figures/distribution_annual_GPP_ER_Q_PAR_T.png', width = 800, height = 400)
aut_sites %>%
  mutate(log10GPP=log10(GPP),
         log10ER=log10(abs(ER)),
         log10NEP=log10(NEP)) %>%
  select(-GPP,-ER,-NEP, -siteyear) %>%
  pivot_longer(-(1:2))%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density() +
  ggtitle("What is the distribution of annual values of GPP, ER, NEP, Q, PAR, and W temp for all site-years?")+
  labs(x="Coefficient of variation (sd/mean * 100)")+
  ggthemes::theme_few()+
  theme(legend.position="none") +
  facet_wrap(~name, scales="free")
dev.off()

# How variable is water temp, discharge, PAR, GPP, ER at these sites?
cv <- function(x) {
  sd(x)/mean(x)*100
}

# Calculate cv for each site-year
cv_aut <- aut %>%
  group_by(site_name, year) %>%
  summarize(across(any_of(ends_with('filled')), cv)) %>%
  rename_with(~str_c("cv_", .), GPP_filled:PAR_filled) 

png('figures/distribution_cv_GPP_ER_Q_PAR_T.png', width = 600, height = 400)
cv_aut %>%
  pivot_longer(-c(1:2)) %>%
  filter(!name=="cv_NEP_filled") %>% #some extremely high basically uninterpretable values; ignore for now
  ggplot(aes(x=abs(value), fill=name)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  ggtitle("What are the ranges of CV of Q, ER, GPP, PAR, and Wtemp among all site-years?")+
  labs(x="Coefficient of variation (sd/mean * 100)")+
  ggthemes::theme_few()+
  theme(legend.position="none")
dev.off()

# End of script.
