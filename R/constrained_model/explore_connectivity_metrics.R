# Test out different derived metrics of connectivity 
# These will be used as drivers in the researcher constrained model

# A Carter
#8/2022

library(tidyverse)

dat <- read_csv('data_working/across_sites_model_data.csv')
mdat <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')

colnames(dat)
# hydrologic connectivity ####
dat <- dat %>%
  mutate(ratio_WA = Width/1000/sqrt(ws_area_km2),
         ratio_PR = log(-ann_GPP_C/ann_ER_C))

wa <- dat %>%
  filter(log(ratio_WA) < -4)%>%
ggplot( aes(log(ratio_WA), Stream_PAR_sum, col = ratio_PR))+
  geom_point(size = 2) +
  scale_color_continuous(type = 'viridis')+
  theme_bw()+
  xlab('Width to Area Ratio')+
  ylab('Light') +
  labs(col = 'log(P/R)')


npp <- ggplot(dat, aes(MOD_ann_NPP, Stream_PAR_sum, col = ratio_PR))+
  geom_point(size = 2) +
  scale_color_continuous(type = 'viridis')+
  theme_bw()+
  xlab('Terrestrial NPP')+
  ylab('Light') +
  labs(col = 'log(P/R)')


pc<-ggplot(dat, aes(PrecipWs, Stream_PAR_sum, col = ratio_PR))+
  geom_point(size = 2) +
  scale_color_continuous(type = 'viridis')+
  theme_bw()+
  xlab('Annual Precipitation')+
  ylab('Light') +
  labs(col = 'log(P/R)')


dam <- ggplot(dat, aes(-log(Dam_densityperkm2), Stream_PAR_sum, col = ratio_PR))+
  geom_point(size = 2) +
  scale_color_continuous(type = 'viridis')+
  theme_bw()+
  xlab('-Dam density')+
  ylab('Light') +
  labs(col = 'log(P/R)')

png('figures/connectivity_light_plots.png',
    width = 600)
  ggpubr::ggarrange(wa, npp, pc, dam, nrow = 2, ncol = 2, common.legend = TRUE,
                    legend = 'right')
dev.off()
ggplot(dat, aes(Pct_impcov, ratio_PR))+
  geom_point()
