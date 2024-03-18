# compare different metrics for quantifying the influence of dams
library(tidyverse)
library(corrplot)

dat <- read_csv('data_working/across_sites_model_data.csv')


# filter out all of the dam and connectivity metrics:

dams <- dat %>% select(site_name, width_to_area, drainage_density_connected, 
               drainage_density, Dam_densityperkm2, Dam_total_vol_m3km2,
               Dam_normal_vol_m3km2, Dam_influence)%>%
  mutate(Dam_influence = case_when(Dam_influence == 0 ~ NA,
                                   Dam_influence == 0.5 ~ 0,
                                   TRUE ~ Dam_influence), 
         connectivity_ratio = drainage_density_connected/drainage_density) %>%
  group_by(site_name) %>%
  summarize(across(everything(), mean, na.rm = T)) %>% 
  select(-site_name)


summary(dams)


dd <- select(dams, -drainage_density, -connectivity_ratio)
d <- cor(dd, use = 'complete.obs')
d1 <- cor.mtest(dd, conf.level = 0.95)

pal1<-colorRampPalette(c('coral2', 'white', 'steelblue'))

png('figures/dam_covariate_correlations.png', width = 4, height = 4,
    res = 300, units = 'in')
par(oma = c(.5,0,0,0))
corrplot(d, type = "lower", order = "AOE", addCoef.col = "black", 
         tl.col = "black", tl.srt = 90, tl.cex = .7,
         method='color', diag=FALSE, cl.cex = .6,
         cl.length = 11, col=pal1(20),
         pch.cex = 1, pch.col='grey20')

dev.off()

