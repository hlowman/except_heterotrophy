# Simulate Respiration data based on:
#   GPP
#   temperature
#   discharge
#   terrestrial carbon input (approx from LAI)

library(tidyverse)
library(lubridate)
setwd('C:/Users/alice.carter/git/except_heterotrophy/')
source('R/data_simulation/simulation_functions.R')
dat <- readRDS('data_356rivers/high_quality_daily_metabolism_with_SP_covariates.rds')

# add site names and subset to autotrophic sites
autosites <- read_csv('data_working/autotrophic_siteyears_daily.csv') %>%
    select(site_name) %>% pull() %>% unique()
md <- read_tsv('data_356rivers/site_data.tsv') %>%
    select(site_name, long_name)
dat <- left_join(dat, md, by = 'site_name') %>%
    filter(site_name %in% autosites) %>%
    mutate(year = year(date))

glimpse(dat)
unique(dat$long_name)

# Simulation for Brandywide Creek nwis_01481500 ####
site <- 'nwis_01481500'
dd <- filter(dat, site_name == site, year !=2007) %>%
    select(site_name, long_name, year, date, GPP_obs = GPP, ER_obs = ER,
           temp.water, discharge, depth, LAI = LAI_proc,
           Stream_PAR = Stream_PAR_sum) %>%
    mutate(GPP_interp = zoo::na.approx(GPP_obs, na.rm = F),
           Q = discharge/max(discharge))

# Simulate litterfall (based on total litterfall/year = 100* max(LAI) gc/m2/y,
# this is arbitrary and should be changed at some point)
dd$litter <- calc_litter_from_LAI(dd)
ggplot(dd, aes(date, LAI)) +
    geom_line(col = 'forestgreen') + geom_line(aes(y = litter)) +
    facet_wrap(.~year, scale = 'free_x') +
    ylim(0,10)


# Calculate antecedent GPP based on desired number of previous days (nweights)
# and the weights (w) on each. sum(w) must equal 1.

nweights <- 5
w <- c(0, 0, 0, .5, .5) # here, antecedent P is a function of 4 and 5 days ago
dd$P_ante <- calc_antecedent_GPP(dd$GPP_interp, nweights, w)
ggplot(dd, aes(date, GPP_obs)) +
    geom_line(col = 'forestgreen') + geom_line(aes(y = P_ante)) +
    facet_wrap(.~year, scale = 'free_x')

AR_f = 0.5# fraction of autotrophic respiration
dd$AR <- -AR_f * dd$GPP_interp

beta_p = 0.8 # the coefficient on antecedent productivity
dd$HR_alg = -beta_p * dd$P_ante
# Calculate underlying terrestrial C and detrital respiration
# temperature dependent respiration rate coefficient
dd$K_er <- calc_rate_coef(dd$temp.water, K_20 = 0.01)

ndays <- nrow(dd)
C0 = 100 # initial organic carbon in stream (g/m2)
beta_s = 0.8# percent C lost to largest storm
sigma_proc = 0.01 # process error (lognormal, % error in each timestep)
sigma_obs = 0.1   # observation error (normally distributed)

dd$C <- dd$HR <- NA_real_ # initialize columns for Carbon and heterotrophic respiration
C_hat <- rep(C0, ndays)
dd$C[1] <- exp(rnorm(1, log(C0), sigma_proc))
dd$HR[1] <- -dd$K_er[1] * dd$C[1]

for(i in 2:ndays){
    C_hat[i] = (dd$C[i-1] + dd$litter[i] + dd$HR[i-1]) * (1 - beta_s * dd$Q[i])
    dd$C[i] = exp(rnorm(1, log(C_hat[i]), sigma_proc))
    dd$HR[i] = -dd$K_er[i] * dd$C[i]
}

dd$ER_latent = dd$AR + dd$HR_alg + dd$HR
dd$ER_sim = rnorm(ndays, dd$ER_latent, sigma_obs)

dd %>% filter(year == 2013) %>%
    mutate(Q = log(discharge)) %>%
    select(date, GPP = GPP_obs, ER_sim, C, temp.water, Q) %>%
    pivot_longer(-date, names_to = 'variable', values_to = 'value') %>%
    ggplot(aes(date, value)) +
        geom_line() +
        facet_wrap(.~variable, scale = 'free', ncol = 1)

write_csv(dd, 'data_working/simulated_ER_nwis_01481500.csv')
# Simulate ER function: ####
# coming soon
