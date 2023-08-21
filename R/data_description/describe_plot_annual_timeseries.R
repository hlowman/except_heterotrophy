# Patterns in Annual Autotrophy
# Alice Carter et al.

# Describe the patterns in metabolism across all siteyears with at least 60% coverage
# Dataset used is the gap-filled data from the Bernhardt data release

# Load packages.
library(tidyverse)

# Load data.
dat <- readRDS('data_356rivers/high_quality_daily_metabolism_with_SP_covariates.rds')

# Descriptive Statistics and data distributions

ann <- dat %>%
  group_by(site_name, year) %>%
  summarize(across(any_of(ends_with('filled')), sum, na.rm = T)) %>%
  rename_with(.cols = ends_with('_filled'), .fn = ~sub('_filled', '', .) )

png('figures/distribution_annual_NEP.png', width = 600, height = 400)
par(mfrow = c(1, 2), mar = c(5,2,3,1))
plot(density(ann$NEP), xlim = c(-6200, 2000), main = '', yaxt = 'n',
     xlab = expression(paste('Annual NEP (g',O[2],'/',m^2,'/d)' )), ylab = '')
mtext('Density', 2, 0.7)
dd = density(ann$NEP, na.rm = T)
ddo = order(dd$x)
xdens = dd$x[ddo]
ydens = dd$y[ddo]
xdens_ut = xdens[xdens >= 0]
ydens_ut = ydens[xdens >= 0]
polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
        col='lightgreen', border='lightgreen')
lines(dd)

dd = density(log(-ann$GPP/ann$ER))
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
             100 * round(length(which(ann$NEP>0))/nrow(ann),2), 
             '% of ', nrow(ann)), side = 3, line = 1, cex = 1.5)
dev.off()

aut_sites <- ann %>% 
  mutate(siteyear = paste(site_name, year, sep = '_')) %>%
  filter(NEP > 0 )
aut <- dat %>% 
  mutate(siteyear = paste(site_name, year, sep = '_')) %>%
  filter(siteyear %in% aut_sites$siteyear) %>%
  group_by(siteyear) %>%
  arrange(DOY) %>%
  mutate(GPP_cum = cumsum(GPP_filled))

dev.off() 
plot(aut$DOY, aut$GPP_cum, type = 'n',  xlab = 'Day of Year',
     ylab = '', ylim = c(0, 15000))
mtext(expression(paste('GPP (g',O[2],m^{-2},d^{-1},')')), 2, 2.2)
for(i in 1:nrow(aut_sites)){
  d <- aut %>% filter(siteyear == aut_sites$siteyear[i]) %>%
    arrange(DOY)
  lines(d$DOY, d$GPP_cum, col = alpha('forestgreen', 0.5))
  
}

# How many UNIQUE sites?
length(unique(aut_sites$site_name))

#Distribution of cumulative(?) water temperature, discharge, PAR, GPP, ER
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


#How variable is water temp, discharge, PAR, GPP, ER at these sites?
cv <- function(x) {
  sd(x)/mean(x)*100
}

#Calculate cv for each site-year
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

