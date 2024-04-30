library(lubridate)
library(broom)
library(scales)
library(furrr)
library(hydrostats)
library(patchwork)
library(xts)
library(viridis)
library(zoo)
library(readxl)
library(changepoint)
library(imputeTS)
library(tidyverse)
library(tseries)


df_met <- readRDS("Data/Loire_DO/metab_extremelyconstrainedK_gppconstrained_all_discharge_bins") %>%
  mutate(GPP = ifelse(GPP < 0, NA, GPP),
         ER = ifelse(ER > 0, NA, ER),
         year = year(date),
         month = month(date)) %>%
  filter(!(year %in% c(2000,2001,2007,2008)),
         between(month, 5, 9))

d_p1 <- df_met %>%
  filter(year == 1995)  %>%
  na.trim() %>%
  na_interpolation(.$GPP, option = "stine") %>%
  na_interpolation(.$ER, option = "stine") %>%
  as.data.frame()
plot(d_p1$GPP)
plot(d_p1$ER)
d_p2 <- df_met %>%
  filter(year == 2015)  %>%
  na.trim() %>%
  na_interpolation(.$GPP, option = "stine") %>%
  na_interpolation(.$ER, option = "stine") %>%
  as.data.frame()
plot(d_p2$GPP)
plot(d_p2$ER)
# Test for stationarity
adf.test(d_p1$GPP)
kpss.test(d_p1$GPP)

diff_p1_min <- diff(-d_p1$ER,1)
diff_p1_max <- diff(d_p1$GPP,1)
plot(diff_p1_max)
plot(diff_p1_min)
diff_p2_min <- diff(-d_p2$ER,1)
diff_p2_max <- diff(d_p2$GPP,1)
plot(diff_p2_max)
plot(diff_p2_min)
adf.test(diff_p1_min)
kpss.test(diff_p1_min)
acf(d_p1$GPP)
acf(diff_p1_max)
ccf(d_p1$GPP, -d_p1$ER)
ccf(d_p2$GPP, -d_p2$ER)
ccf(diff_p1_max, diff_p1_min)
ccf(diff_p2_max, diff_p2_min)

p1_ccf <- ccf(diff_p1_max, diff_p1_min)
p2_ccf <- ccf(diff_p2_max, diff_p2_min)
df_ccf <- data.frame(period = "Period 1", Lag = p1_ccf$lag, CCF = p1_ccf$acf) %>%
  bind_rows(data.frame(period = "Period 2", Lag = p2_ccf$lag, CCF = p2_ccf$acf))

ggplot(df_ccf) +
    geom_col(aes(x = Lag, y = CCF), width = 0.2) +
    facet_wrap(~period) +
    scale_x_continuous(limits = c(-7,7),
                       breaks = seq(-7, 7, 1)) +
    scale_y_continuous(breaks = seq(-0.2, 0.5, 0.1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0.05,
               color = "blue",
               linetype = "dashed") +
    geom_hline(yintercept = -0.05,
               color = "blue",
               linetype = "dashed") +
    ggtitle("Daily max-min cross correlation")


  ggsave(filename = "Figures/max_min_ccf.tiff",
         device = "tiff",
         width = 12,
         height = 8,
         units = "in")