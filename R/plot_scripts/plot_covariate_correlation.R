

library(tidyverse)
library(corrplot)

dat <- read_csv('data_working/across_sites_model_data.csv')
vars <- dat %>%
  group_by(site_name) %>%
  summarize(across(everything(), mean, na.rm = T)) %>%
  select(ann_GPP_C, ann_ER_C, ann_NEP_C, Stream_PAR_sum, Disch_cv,
         Disch_amp, RBI, max_interstorm,
         width_to_area, PrecipWs, MOD_ann_NPP, PAR_kurt, TmeanWs, 
         ElevWs, Width, drainage_density_connected)

res <- cor(vars, use='complete.obs')
res1 <- cor.mtest(vars, conf.level=0.95)


# dd <- vars %>%
#   select(slope, depth, 'log(Q)' = discharge, 'water temp' = watertemp,
#          GPP, ER, ':O[2]' = O2.mgL, 'log(DOC)' = DOC.mgl,
#          ':log(NO[3])' = NO3.mgl)
dd <- select(vars, -drainage_density_connected)
d <- cor(dd, use = 'complete.obs')
d1 <- cor.mtest(dd, conf.level = 0.95)

pal1<-colorRampPalette(c('coral2', 'white', 'steelblue'))

png('figures/covariate_correlations.png', width = 4, height = 4,
     res = 300, units = 'in')
par(oma = c(.5,0,0,0))
corrplot(d, type = "lower", order = "original", p.mat = d1$p,
         tl.col = "black", tl.srt = 90, tl.cex = .7,
         method='color', diag=FALSE, cl.cex = .6,
         cl.length = 11, col=pal1(20),
         sig.level = c(0.01, 0.05, 0.1), insig = 'label_sig',
         pch.cex = 1, pch.col='grey20')
mtext('Correlation Coefficient', 1, 4.5, cex = .7, adj = .57)

dev.off()

