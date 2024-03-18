## Matt Dunkle made this code May 29 2023 ##
## Melissa D. edited on June 27 2023 ##
## Alice. edited on Dec 2023 ##

#data = across_sites_model_data#

library(tidyverse);
library(vegan);
#library(ggvegan);
library(BiodiversityR);
library(ggforce);
library(ggrepel); 
library(here)


across_sites_model_data <- read.csv(here::here("data_working/across_sites_model_data.csv"), 
                                    header=TRUE, fileEncoding = "UTF-8-BOM")

literature_data <- read.csv(here::here("data_working/literature_streams_data_for_PCA.csv"), 
                            header=TRUE, fileEncoding = "UTF-8-BOM")

literature_data2 <- literature_data %>%  mutate(Class = "Literature")

auto_class <- read_csv('data_working/Autotrophic_site_classifications.csv') %>%
  select(-Site, -lat, -lon)

autotrophic = across_sites_model_data %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(max_NEP = max(ann_NEP_C, na.rm = T),
            across(where(is.numeric), \(x) median(x, na.rm = T))) %>% 
  ungroup() %>% 
  select(site_name, lat, lon, ann_GPP_C, ann_ER_C, ann_NEP_C, max_NEP, PR, 
         Stream_PAR_sum, max_interstorm, MOD_ann_NPP, PAR_kurt, ElevWs, Width) %>% 
  # select(site_name, year, lat, lon, ann_GPP_C, ann_ER_C, ann_NEP_C,  PR,PAR_sum, Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs, Width, Stream_PAR_sum) %>% 
  filter(complete.cases(across(ann_GPP_C:Width))) %>% as_tibble() %>% 
  mutate(Class = case_when(PR > 1 ~ "Autotrophic", 
                           max_NEP > 0 ~ "Sometimes Autotrophic",
                           PR < 1 ~ "Heterotrophic")) %>%
  left_join(auto_class, by = 'site_name')



autotrophic_data = autotrophic %>% 
  select( Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, 
          drainage_density_connected, PAR_kurt, ElevWs, Width,Stream_PAR_sum)

autotrophic_data2 = autotrophic %>% 
  select( site_name, year, Disch_cv, max_interstorm, RBI, width_to_area, 
          PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs,
          Width,Stream_PAR_sum)


litauto<-dplyr::bind_rows(autotrophic, literature_data2)

## additional dataframes #

autotrophic_QR = autotrophic %>% 
  select( site_name, Stream_PAR_sum, max_interstorm, Disch_cv, RBI,
          width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected)

autotrophic_sparse = autotrophic %>% 
  select( site_name, year, max_interstorm, MOD_ann_NPP, PAR_kurt, ElevWs, Width)


#auto_hellinger = disttransform(autotrophic_data, method = "log")


pca1 = rda(autotrophic_data, data = autotrophic, scale=T)

plot1 = ordiplot(pca1, choices = c(1,2))

scores<-scores(pca1, display=c("sites"))

sites.long1 <- sites.long(plot1, env.data=autotrophic[,1:8]) #extract data for plotting

head(sites.long1)


spec.envfit <- envfit(plot1, env=autotrophic[,10:19])
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)


species.long1 <- species.long(plot1)
species.long1 #all data




rda_raw <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  
  geom_point(data=sites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=axis1, y=axis2, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1,0.1))+
  
 geom_text_repel(data = sites.long1 %>% filter(PR>1.5), aes(x=axis1, y = axis2, label = site_name))+
  
  geom_mark_ellipse(data=sites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                               PR>1~"Autotrophic")),
                    aes(x=axis1, y=axis2, linetype = Class),
                    expand=0, size=0.2, show.legend=T, alpha =.1) +
  geom_segment(data=species.long1# %>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
                                             # axis2 = case_when(labels == "ElevWs"~1,T~axis2),
                                             # weight = case_when(labels %in% c("ElevWs","Width")~2,T~1))
               ,
               aes(xend=axis1, y = 0,x=0, yend=axis2),
               colour="red",  arrow=arrow(), show.legend = F) +

  geom_text_repel(data=species.long1#%>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
                                     #          axis2 = case_when(labels == "ElevWs"~1,T~axis2))
                  ,
                  aes(x=axis1, y=axis2, label=labels),
                  colour="black", check.overlap = T, show.legend = F)+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))+
  theme_bw()+theme(panel.grid = element_blank())

rda_raw


## MD version ##

autotrophic_2<-autotrophic_data2[,3:13]

pca_all = rda(autotrophic_2,  scale=T)  ## 52 % var 2 axes ##

# with just the sparse predictors ##

autotrophic_sparse2<-autotrophic_sparse[,3:7]

pca_sparse=rda(autotrophic_sparse2, scale=T) ## 61.37 % 2 axes  var ##

plot3<-plot(pca_sparse)

sparsesites.long1 <- sites.long(plot3, env.data=autotrophic[,1:8]) #extract data for plotting

sparsespecies.long1 <- species.long(plot3)

rda_sparse <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  
  geom_point(data=sparsesites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                              PR>1~"Autotrophic")), 
             aes(x=axis1, y=axis2, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1,0.1))+
  
  #geom_text_repel(data = sites.long1 %>% filter(PR>1.5), aes(x=axis1, y = axis2, label = site_name))+
  
  geom_mark_ellipse(data=sparsesites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                                     PR>1~"Autotrophic")),
                    aes(x=axis1, y=axis2, linetype = Class),
                    expand=0, size=0.2, show.legend=T, alpha =.1) +
  geom_segment(data=sparsespecies.long1# %>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
               # axis2 = case_when(labels == "ElevWs"~1,T~axis2),
               # weight = case_when(labels %in% c("ElevWs","Width")~2,T~1))
               ,
               aes(xend=axis1, y = 0,x=0, yend=axis2),
               colour="red",  arrow=arrow(), show.legend = F) +
  
  geom_text_repel(data=sparsespecies.long1#%>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
                  #          axis2 = case_when(labels == "ElevWs"~1,T~axis2))
                  ,
                  aes(x=axis1, y=axis2, label=labels),
                  colour="black", check.overlap = T, show.legend = F)+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))+
  theme_bw()+theme(panel.grid = element_blank())+ggtitle("PCA with predictors from sparse model")




## with just the QR predictors ###

autotrophic_QR2<-autotrophic_QR[,3:10]

scaled_auto_QR<-scale(autotrophic_QR2)

pca_qr = rda(autotrophic_QR2,  scale=T)  ## 61.4 % var 2 axes ##

#

qrscores<-scores(pca_qr, display=c("sites"))

plot2<-plot(pca_qr)

qrsites.long1 <- sites.long(plot2, env.data=autotrophic[,1:8]) #extract data for plotting

qrspecies.long1 <- species.long(plot2)



rda_QR <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  
  geom_point(data=qrsites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=axis1, y=axis2, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1,0.1))+
  
  #geom_text_repel(data = sites.long1 %>% filter(PR>1.5), aes(x=axis1, y = axis2, label = site_name))+
  
  geom_mark_ellipse(data=qrsites.long1 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                                   PR>1~"Autotrophic")),
                    aes(x=axis1, y=axis2, linetype = Class),
                    expand=0, size=0.2, show.legend=T, alpha =.1) +
  geom_segment(data=qrspecies.long1# %>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
               # axis2 = case_when(labels == "ElevWs"~1,T~axis2),
               # weight = case_when(labels %in% c("ElevWs","Width")~2,T~1))
               ,
               aes(xend=axis1, y = 0,x=0, yend=axis2),
               colour="red",  arrow=arrow(), show.legend = F) +
  
  geom_text_repel(data=qrspecies.long1#%>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
                  #          axis2 = case_when(labels == "ElevWs"~1,T~axis2))
                  ,
                  aes(x=axis1, y=axis2, label=labels),
                  colour="black", check.overlap = T, show.legend = F)+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))+
  theme_bw()+theme(panel.grid = element_blank())+ggtitle("PCA with predictors from QR")



##

## WITH BOTH THE LITERATURE AND OUR AUTOROPHIC SITES

head(litauto)

litauto2<-litauto %>% select(site_name,year, ann_GPP_C , PR, Class, River, Type1, Type2, nwis_gage, Disch_cv:Stream_PAR_sum)

litauto3<-litauto2[,10:20]

pca_litauto = rda(litauto3,  scale=T)  ## 48 % var 2 axes ##

summary(pca_litauto) %>% selec

plotlitautopca<-plot(pca_litauto)

litautositeslong <- sites.long(plotlitautopca, env.data=litauto2[,1:9]) #extract data for plotting

litautospecieslong <- species.long(plotlitautopca)

##need to split up the sites dataframes into lit/auto again for plotting purposes###

litsiteslong2<-litautositeslong %>% filter(nwis_gage>0) %>% select(River, Type1, Type2, nwis_gage, axis1, axis2, labels)

autositeslong2<-litautositeslong %>% filter(year>0) %>% select(site_name, year, ann_GPP_C, PR, Class, axis1, axis2, labels)


rda_litauto <-  ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  
  
  geom_point(data=autositeslong2 %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                              PR>1~"Autotrophic")), 
             aes(x=axis1, y=axis2, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1,0.1))+
  
  geom_mark_ellipse(data=autositeslong2%>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                                     PR>1~"Autotrophic")),
                    aes(x=axis1, y=axis2, linetype = Class),
                    expand=0, size=0.2, show.legend=T, alpha =.1) +
  geom_segment(data=litautospecieslong# %>% mutate(axis1 = case_when(labels == "Width"~-1,T~axis1),
               # axis2 = case_when(labels == "ElevWs"~1,T~axis2),
               # weight = case_when(labels %in% c("ElevWs","Width")~2,T~1))
               ,
               aes(xend=axis1, y = 0,x=0, yend=axis2),
               colour="red3",  arrow=arrow(), show.legend = F) +
  
  geom_text(data=litautospecieslong %>% mutate(labels = case_when(labels =="Disch_cv"~ "                                      Discharge (cv)",
                                                                  labels =="max_interstorm"~ "                    Interstorm (max)",
                                                                  labels =="RBI"~ "Flashiness (RBI)",
                                                                  labels =="width_to_area"~ "Width/Area",
                                                                  labels =="PrecipWs"~ "Precipitation",
                                                                  labels =="MOD_ann_NPP"~ "Annual NPP",
                                                                  labels =="drainage_density_connected"~ "Drainage Density                                         ",
                                                                  labels =="PAR_kurt"~ "PAR Kurtosis",
                                                                  labels =="ElevWs"~ "                Watershed Elev.",
                                                                  labels =="Width"~ "Width",
                                                                  labels =="Stream_PAR_sum"~ "        PAR (stream)",
                                                                  T~NA
                                                                  )),
                  aes(x=axis1/1.5, y=axis2/1.5, label=labels, 
                      angle =case_when(labels =="Width/Area"~ -65,T~(axis2/axis1)*atan(
                          # slope
                          1 *
                            # aspect ratio of plot:
                            unit(0, 'npc') %>% grid::convertY('native', valueOnly = T) /
                            unit(1, 'npc') %>% grid::convertX('native', valueOnly = T) /
                            # ratio of y-range to x-range of plot:
                            ( 2 / 2 )
                        ) * 180 / pi)), position = position_nudge(y = c(-.085,#Disch_cv
                                                                       -.085,
                                                                       -.085,
                                                                       -.085,#width_to_area
                                                                       .085,
                                                                       -.085,
                                                                       .085,#drain_dens
                                                                       .085,
                                                                       -.085,
                                                                       -.085,
                                                                       .085),
                                                                  x = c(0, #Disch_cv
                                                                        0,
                                                                        0,
                                                                        0,#width_to_area
                                                                        0,
                                                                        0,
                                                                        0,#drain_dens
                                                                        0,
                                                                        0,
                                                                        0,
                                                                        0)),
                  colour="red3", check.overlap = T, show.legend = F)+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))+
    geom_point(data=litsiteslong2, aes(x=axis1, y=axis2), color="blue", 
               size=2, show.legend = F, shape = 19) +
    geom_text(data = litsiteslong2, 
              mapping = aes(x=axis1, y=axis2, label = River), 
              inherit.aes = FALSE, color="blue", position = position_nudge(y = .085, set.seed(5))) +
  theme_bw()+theme(panel.grid = element_blank())+
    labs(x = "pca1 (32.34%)", y = "pca2 (16.88%)")+coord_cartesian()

save_plot("Heterotrophy_Figure6_PCA.png", rda_litauto, dpi = 600, base_width = 8, base_height = 6)

ggplot(litauto, aes(x = Disch_cv, fill = Class)) + 
  geom_histogram()

