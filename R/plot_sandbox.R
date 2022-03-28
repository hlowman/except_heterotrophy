
setwd('C:/Users/Alice Carter/git/except_heterotrophy/')
library(tidyverse)
library(lubridate)
# library(here)

# read in autotrophic site years:
daily <- read_csv('data_working/autotrophic_siteyears_daily.csv')
annual <- read_csv('data_working/annual_gapfilled_values.csv')

# read in site information
siteinfo <- read_csv('data_working/across_sites_model_data.csv')

# plot NEP vs GPP
png('figures/NEPvsGPP_ERcorrelation.png',  height = 480, width = 640)
annual %>% filter(NEP > 0) %>%
  mutate(gpp_er_cor = case_when(gpp_er_cor < 0 ~ 0,
                                TRUE ~ gpp_er_cor)) %>%
  ggplot(aes(GPP, NEP, col = gpp_er_cor)) +
    geom_point(size = 3) +
    scale_color_gradientn(colors = viridis::magma(5)) +
    scale_x_log10(limits = c(300, 13000))+ 
    scale_y_log10(limits = c(10, 6000))+
    labs(col = "GPP x ER\ncorrelation") +
    xlab("GPP (g O2/m2/y)") + ylab("NEP (g O2/m2/y)") +
    geom_abline(slope = 1, intercept = -0.252, linetype = 'dashed') +
    theme_minimal()
dev.off()

png('figures/NEPvsGPP_cvQ.png',  height = 480, width = 640)
annual %>% filter(NEP > 0) %>%
  ggplot(aes(GPP, NEP, col = disch_cv)) +
    geom_point(size = 3) +
    scale_color_gradientn(colors = viridis::magma(4)) +
    scale_x_log10(limits = c(300, 13000))+ 
    scale_y_log10(limits = c(10, 6000))+
    xlab("GPP (g O2/m2/y)") + ylab("NEP (g O2/m2/y)") +
    geom_abline(slope = 1, intercept = -0.252, linetype = 'dashed') +
    theme_minimal()
dev.off()


auto_dm <- readRDS("data_working/autotrophic_event_durations_magnitudes.rds")

png('figures/short_autoevents_magnitude_seasonally.png',
    height = 480, width = 640 )
auto_dm %>%
  mutate(month = lubridate::month(start_date),
         event_duration = as.numeric(event_duration)) %>%
  filter(event_duration <= 3) %>%
  ggplot(aes(factor(month), NEP/event_duration)) +
  geom_boxplot()+
  ggtitle('Magnitude of 1-3 day autotrophic events by month') +
  scale_y_log10() + ylab('average NEP (g O2/m2/d)') + xlab('month')
dev.off()
png('figures/autotrophic_events_magnitude_by_duration.png',
    height = 480, width = 640 )
auto_dm %>%
  ggplot(aes(factor(event_duration), NEP/as.numeric(event_duration))) +
  geom_boxplot()+
  scale_y_log10() + 
  ggtitle('Magnitude of autotrophic events by duration') +
  ylab('average NEP (g O2/m2/d)') + xlab('event duration (days)')
dev.off()

ggpubr::ggarrange(a,b, ncol = 2)
auto_dm %>%
  mutate(month = lubridate::month(start_date),
         event_duration = as.numeric(event_duration),
         NEP_daily = NEP/event_duration,
         GPP_daily = GPP/event_duration) %>%
  ggplot(aes(event_duration, GPP_daily)) +
  geom_point() +
  geom_point(aes(y = NEP_daily), col = 'grey75')+
  theme_minimal()
 boxplot()+
  scale_y_log10()

# plot of autotrophic site-years by stream order

# select only autotrophic site-years
annual_aut <- annual %>%
   filter(NEP > 0)

# pull out stream order data
site_so <- siteinfo %>%
  select(site_name, lat, lon, NHD_STREAMORDE) %>%
  distinct(site_name, NHD_STREAMORDE, .keep_all = TRUE)

aut_streams <- inner_join(site_so, annual_aut, by = (c("site_name"))) %>%
  mutate(order = factor(NHD_STREAMORDE))

(fig1 <- ggplot(aut_streams, aes(x = order, y = NEP)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Stream Order",
       y = "Net Ecosystem Productivity") +
  theme_minimal())

(fig2 <- ggplot(aut_streams, aes(x = lat, y = NEP)) +
    geom_point() +
    scale_y_log10() +
    labs(x = "Latitude",
         y = "Net Ecosystem Productivity") +
    theme_minimal())

(fig3 <- ggplot(aut_streams, aes(x = lon, y = NEP)) +
    geom_point() +
    scale_y_log10() +
    labs(x = "Longitude",
         y = "Net Ecosystem Productivity") +
    theme_minimal())

# End of script.
