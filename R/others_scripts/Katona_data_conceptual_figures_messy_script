###### working with Modelscape Autotrophy data
###### in effort to create data-driven conceptual figure

### add in watershed_summary_data
watershed_summary_data<-read.csv(file.choose(), header=T)

### add in across_sites_model_data
### this file contains the data used in quantile regression model
across_sites_model_data <- read.csv(file.choose(), header=T)

### high_quality_daily_metabolism_with_SP_covariates(1)

### merge these together
AutoDat <- merge(watershed_summary_data, `high_quality_daily_metabolism_with_SP_covariates(1)`, by = c("site_name"))

### conceptual figure brainstorm:

###  figure that has 3 axes representing the hypothesized driver categories with 
###     the  data for autotrophic sites plotted on top of this
###  What about comparing this to a PCA using more of the drivers? We could 
###     do this with the drivers selected by the sparse model. Maybe a two 
###     panel figure, one with the three hypothesized drivers, and one with a 
###     PCA of the sparse drivers could be a good parallel to the two models.

### hypothesized drivers: light, disturbance, and watershed connectivity 

### what is median max_interstorm?
median(across_sites_model_data$max_interstorm)
## 35 days
## min = 7; max = 365

### move to ggplot2
library(ggplot2)

### Stream PAR ~ width-to-area
plot1 <- ggplot(across_sites_model_data, aes(x = Stream_PAR_sum, y = width_to_area, color = max_interstorm))+
                  geom_point(size = 2, alpha = 0.45)
plot1 + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("stream width:area"))+
  labs(color = " maximum \n interstorm \n interval")

### Stream PAR ~ watershed precipitation 
plot2 <- ggplot(across_sites_model_data, aes(x = Stream_PAR_sum, y = PrecipWs, color = max_interstorm))+
  geom_point(size = 2, alpha = 0.45)
plot2 + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("watershed precipitation"~("mm / year")))+
  labs(color = " maximum \n interstorm \n interval")

### Stream PAR ~ modeled annual terrestrial NPP
plot3 <- ggplot(across_sites_model_data, aes(x = Stream_PAR_sum, y = MOD_ann_NPP, color = max_interstorm))+
  geom_point(size = 2, alpha = 0.45)
plot3 + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("terrestrial NPP"))+
  labs(color = " maximum \n interstorm \n interval")

##### make similar figures with the autotrophic sites only

Autotrophic <- subset(across_sites_model_data, PR >1.0)

### Stream PAR ~ width-to-area
plot1a <- ggplot(Autotrophic, aes(x = Stream_PAR_sum, y = width_to_area, color = max_interstorm))+
  geom_point(size = 2, alpha = 0.45)
plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("stream width:area"))+
  labs(color = " maximum \n interstorm \n interval")

### Stream PAR ~ watershed precipitation 
plot2a <- ggplot(Autotrophic, aes(x = Stream_PAR_sum, y = PrecipWs, color = max_interstorm))+
  geom_point(size = 2, alpha = 0.45)
plot2a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("watershed precipitation"~("mm / year")))+
  labs(color = " maximum \n interstorm \n interval")

### Stream PAR ~ modeled annual terrestrial NPP
plot3a <- ggplot(Autotrophic, aes(x = Stream_PAR_sum, y = MOD_ann_NPP, color = max_interstorm))+
  geom_point(size = 2, alpha = 0.45)
plot3a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_minimal()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("terrestrial NPP"))+
  labs(color = " maximum \n interstorm \n interval")



################################################################################
################################################################################
########## format according to Matt Dunkle's code for PCA: #####################
################################################################################
################################################################################

data = across_sites_model_data
library(tidyverse);library(vegan);library(ggvegan);library(BiodiversityR);library(ggforce);library(ggrepel)


autotrophic = across_sites_model_data %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(across(where(is.numeric), median, na.rm = T)) %>% ungroup() %>% 
  mutate(code = paste(site_name, year, sep = "_")) %>% 
  select(site_name, year, lat, lon, ann_GPP_C, ann_ER_C, ann_NEP_C,  PR,PAR_sum, Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs, Width) %>% 
  filter(complete.cases(across(ann_GPP_C:Width))) %>% as_tibble() %>% 
  mutate(Class = case_when(PR<1~"Heterotrophic",
                           PR>1~"Autotrophic"))

autotrophic_data = autotrophic %>% select( Disch_cv, max_interstorm, RBI, width_to_area, PrecipWs, MOD_ann_NPP, drainage_density_connected, PAR_kurt, ElevWs, Width)

##### Stream PAR ~ width-to-area
plot1a <- ggplot()+
  geom_point(data=autotrophic %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=PAR_sum, y=width_to_area, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
   geom_point(size = 2, alpha = 0.45)+
#plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_bw()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("stream width:area"))+
  labs(color = "P:R")+
  scale_alpha_manual(values = c(1,0.1))+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))
plot1a

### Stream PAR ~ watershed precipitation 
plot2a <- ggplot()+
  geom_point(data=autotrophic %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=PAR_sum, y=PrecipWs, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  geom_point(size = 2, alpha = 0.45)+
  #plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_bw()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("precipitation (mm/year)"))+
  labs(color = "P:R")+
  scale_alpha_manual(values = c(1,0.1))+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))
plot2a

### Stream PAR ~ modeled annual terrestrial NPP
plot3a <- ggplot()+
  geom_point(data=autotrophic %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=PAR_sum, y=MOD_ann_NPP, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  geom_point(size = 2, alpha = 0.45)+
  #plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_bw()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("terrestrial NPP"))+
  labs(color = "P:R")+
  scale_alpha_manual(values = c(1,0.1))+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))
plot3a

### now, let's add max interstorm interval to this
plot4a <- ggplot()+
  geom_point(data=autotrophic %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=PAR_sum, y=max_interstorm, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  geom_point(size = 2, alpha = 0.45)+
  #plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_bw()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("maximum interstorm interval"))+
  labs(color = "P:R")+
  scale_alpha_manual(values = c(1,0.1))+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))
plot4a


### now, let's add elevation!
plot5a <- ggplot()+
  geom_point(data=autotrophic %>%  mutate(Class = case_when(PR<1~"Heterotrophic",
                                                            PR>1~"Autotrophic")), 
             aes(x=PAR_sum, y=ElevWs, colour=PR, shape = Class, alpha = Class), 
             size=5, show.legend = T) +
  geom_point(size = 2, alpha = 0.45)+
  #plot1a + scale_color_gradientn(colours = c("darkred", "red", "darkorange", "orange", "yellow", "white", "lightblue", "blue", "darkblue", "green"), limits = c (7,350))+
  theme_bw()+
  xlab(expression("stream PAR"~("µmol / m² / s")))+
  ylab(expression("elevation (MAS)"))+
  labs(color = "P:R")+
  scale_alpha_manual(values = c(1,0.1))+
  scale_colour_stepsn(colors = terrain.colors(10, rev=T))
plot5a


