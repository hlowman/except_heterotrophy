# Annual autotrophy in streams
# May 26, 2022
# J. Blaszczak, C. Barbosa, M. Desiervo

## load more packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "ggmap","maps","mapdata",
         "ggsn","wesanderson"), require, character.only=T)


## Import Data
NEP_info <- read.csv("")

## Generate map of NEP
NEP_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                             maptype='toner'))+
    geom_point(data = NEP_info, aes(x = lon, y = lat, 
                                  fill=NEP, size=NEP), shape=21)+
    theme(legend.position = "right")+
    labs(x="Longitude", y="Latitude")+
    scale_fill_gradient("Mean annual NEP (units C)",
                        low = "blue", high = "red")+
    scale_size_continuous("Mean annual NEP (units C)")

## Generate map of P:R
PR_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = NEP_info, aes(x = lon, y = lat, 
                                fill=PtoR, size=PtoR), shape=21)+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")+
  scale_fill_gradient("P:R",
                      low = "blue", high = "red")+
  scale_size_continuous("P:R")







