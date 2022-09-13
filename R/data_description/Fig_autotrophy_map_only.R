# Annual autotrophy in streams
# May 26, 2022
# J. Blaszczak, C. Barbosa, M. Desiervo

## load more packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "ggmap","maps","mapdata",
         "ggsn","wesanderson"), require, character.only=T)

library(dplyr)

## Import Data
NEP_info <- read.csv("data_working/riversiteannual.csv")
NEP_info_2<-NEP_info %>% mutate(logPtoP=log(PtoR))


## Generate map of NEP
NEP_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                             maptype='toner'))+
    geom_point(data = NEP_info_2, aes(x = lon, y = lat, 
                                  fill=NEP, size=NEP), shape=21)+
    theme(legend.position = "right")+
    labs(x="Longitude", y="Latitude")+
    scale_fill_gradient("Mean annual NEP (units C)",
                        low = "blue", high = "red")+
    scale_size_continuous("Mean annual NEP (units C)")

## Generate map of P:R
PR_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = NEP_info_2, aes(x = lon, y = lat, 
                                fill=PtoR, size=PtoR), shape=21)+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")+
  scale_fill_gradient("P:R",
                      low = "green", high = "brown")+
  scale_size_continuous("P:R")


PR_map_fig
NEP_map_fig

NEP_quantile <- quantile(NEP_info_2[,5], probs = seq(0, 1, 0.25), na.rm = FALSE,
         names = TRUE)


#library(ggExtra)
#ggMarginal(PR_map_fig, data=NEP_info_2)

NEP_hist <- hist(NEP_info_2[,8],main="logNEP",xlab="",xlim=c(-4,1),col="brown", ylim = c(0,100))



