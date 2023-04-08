
library(dataRetrieval)
library(tidyverse)

syc <- '09510200'
ss <- '02239500'
cf <- '12340500'

syc_flow <- readNWISdata(sites = syc,
                         service = "dv",
                         parameterCd = "00060",
                         startDate = "2017-10-01",
                         endDate = "2018-09-30")
syc_flow$index = seq(1:nrow(syc_flow))
ss_flow <- readNWISdata(sites = ss,
                        service = "dv",
                        parameterCd = "00060",
                        startDate = "2006-10-01",
                        endDate = "2007-09-30")
ss_flow$index = seq(1:nrow(ss_flow))
cf_flow <- readNWISdata(sites = cf,
                        service = "dv",
                        parameterCd = "00060",
                        startDate = "2020-10-01",
                        endDate = "2021-09-30")
cf_flow$index = seq(1:nrow(cf_flow))

plot(cf_flow$X_00060_00003)

flow <- bind_rows(ss_flow, syc_flow, cf_flow) %>%
    rename(Q = X_00060_00003)
sss = 2
syc <- syc_flow %>%
    ggplot(aes(index,X_00060_00003)) +
    geom_line(size = sss) +
    # facet_wrap(site_no~., ncol = 1 , scales = 'free_y') +
    scale_y_log10() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
          axis.text = element_blank()) + ylab("") + xlab("")

cf <- cf_flow %>%
    ggplot(aes(index,X_00060_00003)) +
    geom_line(size = sss) +
    # facet_wrap(site_no~., ncol = 1 , scales = 'free_y') +
    scale_y_log10() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
          axis.text = element_blank()) + ylab("") + xlab("")
ss <- ss_flow %>%
    ggplot(aes(index,X_00060_00003)) +
    geom_line(size = sss) +
    # facet_wrap(site_no~., ncol = 1 , scales = 'free_y') +
    scale_y_log10() +
    theme_minimal() +
    ylim(1, 500) +
    theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
          axis.text = element_blank()) + ylab("") + xlab("")
ss2 <- ss_flow %>%
    ggplot(aes(index,X_00060_00003)) +
    geom_line(size = sss, col = 'white') +
    # facet_wrap(site_no~., ncol = 1 , scales = 'free_y') +
    scale_y_log10() +
    theme_minimal() +
    ylim(1, 500) +
    theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
          axis.text = element_blank()) + ylab("") + xlab("")

ggpubr::ggarrange(ss, ss2, ss2, cf, syc, ncol = 1)
