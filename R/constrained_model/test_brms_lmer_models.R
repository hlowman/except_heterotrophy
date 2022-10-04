# Test hierarchical model for across sites 
# pooled by site as there was insufficient evidence of variation across siteyears
# A Carter

library(tidyverse)
library(brms)
library(lme4)

dat <- read_csv('data_working/across_sites_model_data.csv')

dat$PR = -dat$ann_GPP_C/dat$ann_ER_C
dat$NEP = dat$ann_GPP_C + dat$ann_ER_C

dd <- dat %>%
    filter(!is.na(PR),
           !is.na(drainage_density_connected),
           !is.infinite(drainage_density_connected),
           !is.na(Stream_PAR_sum), 
           !is.na(MOD_ann_NPP)) %>%
    select(site_name, PR, NEP, drainage_density_connected, Stream_PAR_sum,
           MOD_ann_NPP, Disch_cv) %>%
    mutate(log_PR = log(PR),
           across(-site_name, ~scale(.)[, 1]))

b <- brm(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected + 
      (1|site_name), data = dd)

b <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected + 
      (1|site_name), data = dd)

l <- lmer(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + 
            drainage_density_connected + (1|site_name), data = dd)

l <- lmer(NEP ~ Stream_PAR_sum + MOD_ann_NPP + 
            drainage_density_connected + (1|site_name), data = dd)

summary(l)


sum(is.na(dat$MOD_ann_NPP))
