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

col_pal <- c(terrain.colors(10)[2],
             terrain.colors(10)[8],
             'darkgreen')
fill_pal <- c(
  # grey.colors(3, alpha = 0.4)[1],
              terrain.colors(10, alpha = 0.5)[2],
              terrain.colors(10, alpha = 0.5)[8],
              alpha('darkgreen', alpha = 0.4))
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

dd <- select(dat,
       site_name, River, Type1, Type2, trophic, lat, lon, year, 
       ann_GPP_C, ann_ER_C, ann_NEP_C, PR,
       Stream_PAR_sum, max_interstorm, Disch_cv, RBI, 
       width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected,
       TmeanWs, PAR_kurt, ElevWs, Disch_amp, Width) 

dd %>%
  pivot_longer(cols = c('Stream_PAR_sum', 'max_interstorm', 'Disch_cv', 'RBI', 
       'width_to_area', 'PrecipWs', 'MOD_ann_NPP', 'drainage_density_connected',
       'TmeanWs', 'PAR_kurt', 'ElevWs', 'Disch_amp', 'Width')) %>%
ggplot( aes(value, group = trophic, fill = trophic))+
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual('Category', values = fill_pal) +
  facet_wrap(.~name, scales = 'free') +
  theme_bw()

# dd %>%
#   pivot_longer(cols = c('Stream_PAR_sum', 'max_interstorm', 'Disch_cv', 'RBI',
#        'width_to_area', 'PrecipWs', 'MOD_ann_NPP', 'drainage_density_connected',
#        'TmeanWs', 'PAR_kurt', 'ElevWs', 'Disch_amp', 'Width')) %>%
# ggplot( aes(value, y = trophic, fill = trophic))+
#   geom_density_ridges(adjust=1.5, alpha=.4) +
#   theme_ridges() + 
#   
#   scale_fill_manual('Category', values = fill_pal) +
#   facet_wrap(.~name, scales = 'free') +
#   theme_bw()


plot_comb_dens <- function(dat, var, name, log = FALSE){
  dat <- data.frame(dat)
  add_dens_plot <- function(d, med, fill, line_col = 'black', lit = FALSE){
    x <- c(d$x, rev(d$x))
    y <- c(d$y, rep(0, length(d$y)))
    y_med <- y[which.min(abs(x - med))]
    if(lit)  polygon(x,y, col = fill, density = 30, border = FALSE)
    if(!lit) polygon(x,y, col = fill, border = FALSE)
    segments(med, 0, med, y_med, col = line_col)
    lines(d$x, d$y, col = line_col)
  }
  
  d_aut <- density(dat[dat$trophic == 'Autotrophic', var], na.rm = T)
  d_het <- density(dat[dat$trophic == 'Heterotrophic', var], na.rm = T)
  d_lit <- density(dat[dat$trophic == 'Literature Sites', var], na.rm = T)
  med_aut <- median(dat[dat$trophic == 'Autotrophic', var], na.rm = T)
  med_het <- median(dat[dat$trophic == 'Heterotrophic', var], na.rm = T)
  med_lit <- median(dat[dat$trophic == 'Literature Sites', var], na.rm = T)
  
  if(log) {
  d_aut <- density(log10(dat[dat$trophic == 'Autotrophic', var]), na.rm = T)
  d_het <- density(log10(dat[dat$trophic == 'Heterotrophic', var]), na.rm = T)
  d_lit <- density(log10(dat[dat$trophic == 'Literature Sites', var]), na.rm = T)
  med_aut <- median(log10(dat[dat$trophic == 'Autotrophic', var]), na.rm = T)
  med_het <- median(log10(dat[dat$trophic == 'Heterotrophic', var]), na.rm = T)
  med_lit <- median(log10(dat[dat$trophic == 'Literature Sites', var]), na.rm = T)
  }
  x_range <- range(d_aut$x, d_het$x, d_lit$x)
  y_range <- range(d_aut$y, d_het$y, d_lit$y)
  y_range[2] <- y_range[2]*1.3
  plot(1, xlim = x_range, ylim = y_range, type = 'n',
       ylab = 'Density', xlab = var, xaxt = 'n')
  if(log){
    axis(1, at = seq(ceiling(x_range[1]), floor(x_range[2])),
         labels = 10 ^ seq(ceiling(x_range[1]), floor(x_range[2])))
  } else {
    axis(1)
  }
  add_dens_plot(d_lit, med_lit, fill_pal[3], col_pal[3], lit = TRUE)
  add_dens_plot(d_aut, med_aut, fill_pal[1], col_pal[1])
  add_dens_plot(d_het, med_het, fill_pal[2], col_pal[2])
  mtext(name, line = -1.3, adj = 0.05)
}

col_pal <- c('forestgreen',
             'sienna4',
             'darkgreen')
png(file = 'figures/covariate_densities.png', width = 8.5, height = 6.75, type = 'cairo',
    res = 300, units = 'in')
# par(mfrow = c(4,1),
#     mar = c(0.2,0,0.2,0),
#     oma = c(3,3,2,2))
par(mfrow = c(4,3),
    mar = c(2.5,2.1,1,0.1),
    oma = c(2.5,3.5,1,3.5))

# plot.new()
box_coords <- par("usr")  # Get the user coordinates of the plotting area
# rect(box_coords[1], box_coords[3], box_coords[2], box_coords[4], border = 'red', lwd = 2)
plot_comb_dens(dat, 'Stream_PAR_sum', 'Stream PAR')
rect(65.5, -.16, 68, 0.48, col = "#FBA59F", lwd = 2, xpd = NA)
rect(-4.2, -.16, 68, 0.48, border = "black", lwd = 2, xpd = NA)
rect(65.5, -.8, 68, -.16, col =  "#D0ED96", lwd = 2, xpd = NA)
rect(-4.2, -.8, 68, -.16, border =  "black", lwd = 2, xpd = NA)
rect(65.5, -1.44, 68, -0.8, col = "#75E5E5", lwd = 2, xpd = NA)
rect(-4.2, -1.44, 68, -0.8, border = "black", lwd = 2, xpd = NA)
rect(65.5, -2.08, 68, -1.44, col = 'grey', lwd = 2, xpd = NA)
rect(-4.2, -2.08, 68, -1.44, border = 'black', lwd = 2, xpd = NA)
plot_comb_dens(dat, 'PAR_kurt', 'PAR Kurtosis')
plot_comb_dens(dat, 'Width', 'Width', log = TRUE)
mtext('Light', side = 4, line = 1.3, las = 0, cex = 1.2, adj = 0.3)

# disturbance
plot_comb_dens(dat, 'max_interstorm', 'Max Interstorm', log = TRUE)
plot_comb_dens(dat, 'Disch_cv', 'Discharge CV')
plot_comb_dens(dat, 'RBI', 'Flashiness (RBI)') 
mtext('Disturbance', side = 4, line = 1.3, las = 0, cex = 1.2, adj = 1)

# connectivity
plot_comb_dens(dat, 'MOD_ann_NPP', 'Terrestrial NPP')
plot_comb_dens(dat, 'PrecipWs', 'Precipitation')
plot_comb_dens(dat, 'width_to_area', 'Width:Area', log = TRUE)
mtext('Connectivity', side = 4, line = 1.3, las = 0, cex = 1.2, adj = 1)


# other
plot_comb_dens(dat, 'ElevWs', 'Elevation', log = TRUE)
plot_comb_dens(dat, 'TmeanWs', 'Air Temperature')
plot_comb_dens(dat, 'Disch_amp', 'Discharge Range')
mtext('Other', side = 4, line = 1.3, las = 0, cex = 1.2, adj = 0.35)

par(mfrow = c(1,1), new = T)
mtext('Covariate Density', side = 2, line = 4.4, cex = 1.2, adj = 0.3)
mtext('Covariate Range', side = 1, line = 3.8, cex = 1.2, adj = 0.5)

dev.off()



# Light

plot_comb_dens(dat, 'Stream_PAR_sum')
plot_comb_dens(dat, 'PAR_kurt')
plot_comb_dens(dat, 'Width', log = TRUE)


# disturbance
plot_comb_dens(dat, 'max_interstorm', log = TRUE)
plot_comb_dens(dat, 'Disch_cv')
plot_comb_dens(dat, 'RBI') 

# connectivity
plot_comb_dens(dat, 'MOD_ann_NPP')
plot_comb_dens(dat, 'PrecipWs')
plot_comb_dens(dat, 'width_to_area', log = TRUE)
plot_comb_dens(dat, 'drainage_density_connected', 'Connected drainage density')

# other
plot_comb_dens(dat, 'ElevWs', log = TRUE)
plot_comb_dens(dat, 'TmeanWs')
plot_comb_dens(dat, 'Disch_amp')




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
