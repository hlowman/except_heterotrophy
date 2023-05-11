data = across_sites_model_data
library(tidyverse);library(vegan);library(ggvegan);library(BiodiversityR);library(ggforce);library(ggrepel)


autotrophic = across_sites_model_data %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(across(where(is.numeric), median, na.rm = T)) %>% ungroup() %>% 
  mutate(code = paste(site_name, year, sep = "_")) %>% 
  select(site_name, year, lat, lon, ann_GPP_C, ann_ER_C, ann_NEP_C,  PR,PAR_sum, Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs, Width, Stream_PAR_sum) %>% 
  filter(complete.cases(across(ann_GPP_C:Width))) %>% as_tibble() %>% 
  mutate(Class = case_when(PR<1~"Heterotrophic",
                           PR>1~"Autotrophic"))



autotrophic_data = autotrophic %>% select( Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs, Width,Stream_PAR_sum)
auto_hellinger = disttransform(autotrophic_data, method = "log")


pca1 = rda(autotrophic_data, data = autotrophic, scale=T)

plot1 = ordiplot(pca1, choices = c(1,2))


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



watershed_summary_data %>% filter(site_name == "nwis_13173600") %>% view()      
