#### attempt at making Fig 1 for heterotrophy paper##

library(here)
library(plyr)
library(dplyr)

riversiteyears<-read.csv("data_working/across_sites_model_data.csv")

riversiteyears_2<-riversiteyears %>% dplyr::select(site_name, year, lat, lon, ann_GPP_C, ann_ER_C) %>% mutate(NEP=(abs(ann_GPP_C) - abs(ann_ER_C))) %>% mutate(PtoR = abs(ann_GPP_C) / abs(ann_ER_C))

riversite<- riversiteyears_2 %>%  group_by(site_name, lat, lon)  %>% summarise(NEP=mean(NEP), PtoR=mean(PtoR), nyears = n())

riversite<-write.csv("data_working/riversiteannual.csv")

                                                                               