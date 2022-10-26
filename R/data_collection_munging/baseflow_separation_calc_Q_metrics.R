# add baseflow to the daily datasets and calculate max time period between floods
# for each site year based on different flood thresholds.

library(tidyverse)
# install.packages('hydrostats')
# library(hydrostats)
# the package EcoHydRology is no longer available on Cran, but this is an old version
# https://cran.r-project.org/src/contrib/Archive/EcoHydRology/
# install.packages('C:/Users/alice.carter/Downloads/EcoHydRology_0.4.12.1.tar.gz', 
#                  repos = NULL, type ='source', dependencies = TRUE)

# library(EcoHydRology)
install.packages("FlowScreen")
library(FlowScreen)


lotic_gap_filled <- readRDS('data_ignored/Bernhardt_2022/lotic_gap_filled.rds')

# after some exploration and reading about baseflow separation methods, 
# It seems that the Eckhardt method is recommended, so we'll use that.
add_baseflow_column <- function(dat){
    dat$bf <- FlowScreen::bf_eckhardt(dat$Disch_filled, a = 0.975, 
                                  BFI = 0.8)
    dat$bfi <- dat$bf/dat$Disch_filled
    
    return(dat)
    # plot(dat$Date, dat$Disch_filled, type = 'l')
    # lines(dat$Date, dat$ek, col = 2)
    # lines(dat$Date, dat$bf, col = 2, lty = 2)
}
# dat <- lotic_gap_filled[[2]]
# add_baseflow_column <- function(dat){
#     dd <- rename(dat, Q = Disch_filled)
#     dd <- hydrostats::baseflows(dd, a = 0.925, ts = 'daily',
#                                 n.reflected = 30)
#     dat <- left_join(dat, select(dd, -Q))
#     plot(dat$Date, dat$discharge, type = 'l')
#     lines(dat$Date, dd$bf, col = 2)
# }
# add_baseflow_column <- function(dat){
#     # dd <- rename(dat, Q = Disch_filled)
#     dd <- EcoHydRology::BaseflowSeparation(streamflow = dat$Disch_filled,
#                                           filter_parameter = 0.975)
#     dat <- bind_cols(dat, dd)
#     plot(dat$Date, dat$Disch_filled, type = 'l')
#     lines(dat$Date, dd$bt, col = 2)
#     lines(dat$Date, dat$bf, col = 2, lty = 2)
# }

lotic_gap_filled <- lapply(lotic_gap_filled, add_baseflow_column)
lotic_siteyears <- lapply(lotic_gap_filled, function(x) group_split(x, Year))

lotic_siteyears_split <- unlist(lotic_siteyears, recursive = F, use.names = T)
names(lotic_siteyears_split) <- sapply(lotic_siteyears_split, 
                                       function(x) paste(x$Site_ID[1], x$Year[1], sep = "_"))

lotic_siteyears_split <- lapply(lotic_siteyears_split, as.data.frame)

calc_storm_gaps <- function(dat, threshold = 0.8){
    dat$storm = case_when(dat$bfi < threshold ~ 1,
                          TRUE ~ 0)    
    storms <- rle(dat$storm)
    storms <- data.frame(vals = storms$values,
                         lengths = storms$lengths)
    N = nrow(storms)
    
    # count interstorm intervals that wrap around from Dec -> Jan
    # if(storms$vals[1] == 0 & storms$vals[N] == 0){
    #    storms$lengths[1] = storms$lengths[1] + storms$lengths[N]
    #    storms <- storms[1:(N-1),]
    # }
    storms <- storms[2:N,]
    storms <- dplyr::filter(storms, vals == 0) 
    s <- data.frame(bfi_threshold = threshold,
                    med_interstorm = median(storms$lengths),
                    max_interstorm = max(storms$lengths))
    
    return(s)
    
}

plot_bf_separation<- function(dat){
    par(mfrow = c(2,1),
        mar = c(2,4,1,2), 
        oma = c(3,1,1,0))
    plot(dat$Date, dat$Disch_filled, type = 'l', ylab = 'discharge', 
         main = dat$Site_ID[1])
    lines(dat$Date, dat$bf, col = 2)
    
    bf <- data.frame()
    for(threshold in seq(0.5, 0.9, by = 0.01)){
        d <- calc_storm_gaps(dat, threshold)
        bf <- bind_rows(bf, d)
    }
    
    plot(bf$bfi_threshold, bf$max_interstorm, type = 'l',
         ylim = range(c(bf$max_interstorm, bf$med_interstorm)),
         ylab = 'max interstorm days')
    lines(bf$bfi_threshold, bf$med_interstorm, lty = 2)
    mtext('baseflow index storm threshold', 1, 2.5)
}    

# for(i in 1:length(lotic_siteyears_split)){
#   plot_bf_separation(lotic_siteyears_split[[923]])
# }

# based on looking through many of these plots, 0.75 seems like a reasonable
# BFI cutoff for what counts as a 'storm'

# calc storm intervals:
bf <- data.frame()
for(i in 1:length(lotic_siteyears_split)){
    dat <- lotic_siteyears_split[[i]]
    d <- calc_storm_gaps(dat, threshold = 0.75)
    d$Site_ID <- dat$Site_ID[1]
    d$Year = dat$Year[1]
  
    bf <- bind_rows(bf, d)
}

bf <- relocate(bf, c(Site_ID, Year))

write_csv(bf, 'data_ignored/annual_interstorm_intervals.csv')
