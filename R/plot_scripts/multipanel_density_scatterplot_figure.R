################################################################################
# Try a multiplot version
# Pair-wise Scatter Plots
library(tidyverse)
library(GGally)


across_sites_model_data <- read_csv('data_working/across_sites_model_data.csv')
lit_sites <- read_csv('data_working/literature_streams_data_for_PCA.csv')

dat <- lit_sites %>%
  mutate(trophic = 'Literature Sites',
         # Stream_PAR_sum = Stream_PAR_sum/365/1000, # convert units to match
         # MOD_ann_NPP = MOD_ann_NPP * 1000,
         site_name = paste('nwis', nwis_gage, sep = '_')) %>%
  select(-nwis_gage) %>%
  bind_rows(mutate(across_sites_model_data,
                   trophic = case_when(PR > 1 ~ 'Autotrophic',
                                       PR < 1 ~ 'Heterotrophic')))

dat_panels <- dat %>%
  select(trophic, Type = Type1, Stream_PAR_sum, max_interstorm, MOD_ann_NPP) %>%
  mutate(max_interstorm = log10(max_interstorm)) %>%
  mutate(Type = case_when(is.na(Type) ~ 'Study',
                          TRUE ~ Type),
         alpha = case_when(trophic == 'Literature Sites' ~ 1,
                           trophic == 'Autotrophic' ~ 0.8,
                           trophic == 'Heterotrophic' ~ 0.8),
         trophic = factor(trophic, levels = c('Literature Sites', 'Autotrophic',
                                              'Heterotrophic'))) %>%
  
  arrange(desc(trophic))

dat_lit <- filter(dat_panels, trophic == 'Literature')

col_pal <- c(grey.colors(3)[1],
             terrain.colors(10, alpha = 0.8)[2],
             terrain.colors(10, alpha = 0.6)[8])
fill_pal <- c(grey.colors(3, alpha = 0.4)[1],
              terrain.colors(10, alpha = 0.5)[2],
              terrain.colors(10, alpha = 0.5)[8])
png('figures/correlation_densities_plot.png',
    width = 6.5, height = 6.5, units = 'in', res = 300)

    ggpairs(dat_panels, columns = 3:5, 
            mapping = ggplot2::aes(color = trophic),
            upper = 'blank',
            lower = list(continuous = wrap('points')),
            columnLabels = c('Light', 'Max interstorm', 'Terrestrial NPP'),
            axisLabels = 'show',
            # switch = 'both'#, legend = c(1,3)
            ) +
        scale_color_manual(values = col_pal)+
        scale_fill_manual(values = fill_pal)+
        # geom_point(data = dat_lit, aes(pch = Type), col = col_pal[3])+
        theme_classic()+
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, 'line'))
    
dev.off()


# Other versions
dat_panels <- dat %>%
  select(trophic, Stream_PAR_sum, max_interstorm, Disch_cv, RBI, 
         MOD_ann_NPP, PrecipWs, width_to_area, drainage_density_connected) %>%
  mutate(max_interstorm = log10(max_interstorm)) %>%
  mutate(alpha = case_when(trophic == 'Literature Sites' ~ 1,
                           trophic == 'Autotrophic' ~ 0.8,
                           trophic == 'Heterotrophic' ~ 0.8),
         trophic = factor(trophic, levels = c('Literature Sites', 'Autotrophic',
                                              'Heterotrophic'))) %>%
  
  arrange(desc(trophic))

png('figures/correlation_densities_plot_all_vars.png',
    width = 10.5, height = 10.5, units = 'in', res = 300)
  ggpairs(dat_panels, columns = 2:9, 
          ggplot2::aes(color = trophic),
          upper = 'blank') +
    scale_color_manual(values = col_pal)+
    scale_fill_manual(values = fill_pal)+
    theme(panel.border = element_rect(fill = NULL))+
    theme_classic()
dev.off()
