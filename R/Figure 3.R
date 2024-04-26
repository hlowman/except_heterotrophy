library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(cowplot)
library(here)

mod_tableNEP<-read.csv('data_working/constrained_quantile_regression_results_NEP.csv') %>%
  mutate(Metric = 'NEP')
mod_tablePR<-read.csv('data_working/constrained_quantile_regression_results_PR.csv') %>%
  mutate(Metric = 'PR')


mod_tableNEP <- bind_rows(mod_tableNEP, mod_tablePR) %>% 
  mutate(modeltype= case_when((!is.na(D_mean) & !is.na(C_mean) ~ "Light + Dist. + Conn."), 
                              (!is.na(C_mean) & is.na(D_mean) ~ "Light + Conn."),
                              (!is.na(D_mean) & is.na(C_mean) ~ "Light + Dist."),
                              (is.na(D_mean) & is.na(C_mean)~ "Light")))


mod_tableNEP$modeltype <- factor(mod_tableNEP$modeltype,
                                 levels=c("Light", "Light + Dist.", 
                                          "Light + Conn.", 
                                          "Light + Dist. + Conn."))


mod_tableNEP <- mod_tableNEP %>%  
  arrange(Metric, modeltype) %>% 
  mutate(modelorder=c(seq(1:20), seq(1:20)))

mod_tableNEPcoeflong <- melt(mod_tableNEP,
                  id.vars=c("mod", "modelorder", "modeltype", 
                            "D_var", "C_var", "AIC", "Metric"),
                  measure.vars=c("L_mean", "D_mean", "C_mean" ),
                  variable.name="covariate",
                  value.name="coef")


mod_tableNEPerrorlong <- melt(mod_tableNEP,
                             id.vars=c("mod", "modelorder", "modeltype", 
                                       "D_var", "C_var", "AIC", "Metric"),
                             measure.vars=c("L_se", "D_se", "C_se"),
                             variable.name="covariate2",
                             value.name="error")

mod_tablelong<-merge(mod_tableNEPcoeflong,mod_tableNEPerrorlong) 

mod_tablelong2<-mod_tablelong %>% drop_na(error)


mod_tablelong2$covariate <- recode(mod_tablelong2$covariate, C_mean = 'Connectivity', 
                          D_mean = 'Disturbance',
                          L_mean = 'Light')

mod_tablelong2$covariate <- factor(mod_tablelong2$covariate,
                                   levels=c("Light", "Disturbance", "Connectivity"))


#### making the figure ####

fig3 <- mod_tablelong2 %>% 
  as_tibble() %>% 
  mutate(dAIC = as.numeric(round(AIC-min(AIC), 1))) %>% 
  mutate(modeltype = as.factor(modeltype)) %>% 
  # mutate(Metric = case_when(is.na(coef)==T~"PR",T~"NEP")) %>%  
  ungroup() %>% 
  ggplot() + 
  geom_rect(aes(fill = modeltype,alpha = modeltype),xmin = -Inf,xmax = Inf,
                       ymin = -Inf,ymax = Inf, show.legend = F) +
  scale_alpha_manual(values = c(.25,.05,.05,.01))+
  geom_point(aes(x = coef ,  y = as.factor(desc(modelorder)),
                 size = -dAIC, shape =Metric), show.legend = T)+
  geom_errorbarh(aes(xmin = coef-error, xmax = coef+error, 
                     y = as.factor(desc(modelorder)),
                     height = 0, color=modeltype), show.legend = F)+
  facet_grid(modeltype~covariate, scales = "free_y", drop = T,
             switch = "y", space = "free")+
  geom_vline(xintercept=0, linetype = "dashed")+
  
  geom_point(aes(x = coef , y = as.factor(desc(modelorder)),  size = -(dAIC), 
                 color=modeltype), show.legend = F)+
  scale_shape_manual(values = c(19,1))+
  scale_size_continuous(range = c(.5,8), guide = T)+
  theme(legend.position = "bottom")+
  theme(panel.spacing.x  =unit(0.04, "lines"), panel.spacing.y = unit(0,"lines"),
        panel.border = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y=element_blank())+
  theme_bw()+theme(panel.grid=element_blank(),strip.background = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.text.y.left  = element_text(angle = 0, hjust = 1),
        legend.position = "bottom")+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+ 
  guides(size = guide_legend(expression(paste(delta, 'AIC')), 
                             direction = "horizontal", nrow = 1, reverse = T),
    shape = guide_legend(override.aes = list(size = 5), title = "Metric"), 
    color = guide_none(), fill = guide_none(), alpha = guide_none())

  
  
fig3


## Updated Figure 3 ##
#Figure 3 changes 11/02/2023 (MD),  12/5/2023 (CT), and 1/8/24 (CT)

fig3.2 <- mod_tablelong2 %>% as_tibble() %>% 
  
  filter(Metric == 'NEP') %>%
  mutate(dAIC = as.numeric(round(min(AIC)-AIC, 1))) %>% 
  mutate(modeltype = as.factor(modeltype)) %>% 
  # mutate(Metric = case_when(is.na(coef)==TRUE ~ "PR",
  #                           TRUE ~ "NEP")) %>% 
  ungroup() %>% 
  mutate(L_var = "Light") %>% 
  mutate(best = case_when( dAIC >= -2 ~ 'Best',
                           dAIC < -2 ~ 'Other' )) %>%
  
  filter(is.na(coef) == F) %>% 
  mutate(se = case_when(covariate == "Light" & covariate2 == "L_se"|
                          covariate == "Disturbance" & covariate2 == "D_se"|
                          covariate == "Connectivity" & covariate2 == "C_se"~ error, T~NA)) %>% 
  mutate(D_var = case_when(D_var == "Disch_cv"~ "Discharge (cv)",
                           D_var == "log_RBI" ~ "Flashiness (RBI)",
                           D_var == "max_interstorm" ~ "Max Interstorm"),
         C_var = case_when(C_var == "drainage_density_connected"~"Drainage Density",
                           C_var == "width_to_area"~"Width:Area",
                           C_var == "MOD_ann_NPP"~ "Terrestrial NPP",
                           C_var == "PrecipWs" ~ "Precipitation")) %>% 
  
  ggplot() + 
  geom_rect(aes(fill = modeltype, alpha = modeltype),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, show.legend = F) +
  scale_alpha_manual(values = c(.15,.03,.03,.005))+
  
  geom_errorbarh(aes(xmin = coef-se, xmax = coef+se, 
                     y = as.factor(modelorder), height= 0, color = best), 
                 show.legend = F)+
  
  geom_point(aes(x = coef ,  y = as.factor(modelorder),
                 color = best,
                 shape =case_when(covariate =="Light"~L_var,
                                  covariate == "Connectivity"~C_var,
                                  covariate == "Disturbance" ~ D_var)), 
             show.legend = T, size = 4)+
  
  facet_grid((modeltype)~covariate, scales = "free_y", drop = T, 
             switch = "y", space = "free")+
  geom_vline(xintercept=0, linetype = "dashed")+
  
  scale_shape_manual(values = c(15,13,17,8,19,7,10,9,10))+  
  scale_color_manual(values = c("black","grey70"))+
  
  theme(legend.position = "bottom")+
  theme(panel.spacing.x  =unit(0.04, "lines"), panel.spacing.y = unit(0,"lines"),
        panel.border = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y=element_blank())+
  theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 30))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.text.y.left  = element_text(angle = 90),
        legend.position = "bottom")+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+ 
  guides(size = guide_legend(expression(paste(delta, 'AIC')), 
                             direction = "vertical", nrow = 1, reverse = T),
         shape = guide_legend(override.aes = list(size = 5), title = "Metric"), 
         color = guide_none(), fill = guide_none(), alpha = guide_none())+
  theme(text = element_text(size=16))

fig3.2 

cowplot::save_plot(here("figures/Figure3_NEP.png"), 
                   fig3.2, base_width = 9, base_height = 10, dpi = 600)
cowplot::save_plot(here("figures/Figure3_PR.png"), 
                   fig3.2, base_width = 9, base_height = 10, dpi = 600)


## Updated Figure 3 ##
#Figure 3 changes 4/22/2023 (AC)


fig3.3 <- mod_tablelong2 %>% as_tibble() %>% 
  
  # filter(Metric == 'NEP') %>%
  filter(Metric == 'PR') %>%
  mutate(dAIC = as.numeric(round(min(AIC)-AIC, 1))) %>% 
  mutate(modeltype = as.factor(modeltype)) %>% 
  ungroup() %>% 
  mutate(L_var = "Light") %>% 
  mutate(best = case_when( dAIC >= -2 ~ 'Best',
                           dAIC < -2 ~ 'Other' )) %>%
  
  filter(is.na(coef) == F) %>% 
  mutate(se = case_when(covariate == "Light" & covariate2 == "L_se"|
                          covariate == "Disturbance" & covariate2 == "D_se"|
                          covariate == "Connectivity" & covariate2 == "C_se"~ error, T~NA)) %>% 
  mutate(D_var = case_when(D_var == "Disch_cv"~ "Discharge (cv)",
                           D_var == "log_RBI" ~ "Flashiness (RBI)",
                           D_var == "max_interstorm" ~ "Max Interstorm"),
         C_var = case_when(C_var == "drainage_density_connected"~"Drainage Density",
                           C_var == "width_to_area"~"Width:Area",
                           C_var == "MOD_ann_NPP"~ "Terrestrial NPP",
                           C_var == "PrecipWs" ~ "Precipitation"),
         shape =case_when(covariate =="Light"~L_var,
                          covariate == "Connectivity"~C_var,
                          covariate == "Disturbance" ~ D_var),
         shape = factor(shape, levels = c('Light', 'Precipitation', 
                                          'Max Interstorm', 'Width:Area', 
                                          'Flashiness (RBI)', 'Terrestrial NPP', 
                                          'Discharge (cv)', 'Drainage Density'
                                          ))) %>% 
  
  ggplot() + 
  geom_rect(fill = NA, xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, show.legend = F) +
  geom_rect(
    data=data.frame(modeltype = factor('Light + Dist. + Conn.')),
    aes(xmin = -Inf, xmax = Inf, ymin = 18.5, ymax = Inf),
    # aes(xmin = -Inf, xmax = Inf, ymin = 17.5, ymax = Inf),
    alpha = 0.5, fill = 'grey') +
  geom_rect(
    data=data.frame(modeltype = factor('Light')),
    aes(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5),
    alpha = 0.5, fill = 'white') +
  geom_rect(
    data=data.frame(modeltype = factor('Light + Dist.')),
    aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 4.5),
    alpha = 0.5, fill = 'white') +
  geom_rect(
    data=data.frame(modeltype = factor('Light + Conn.')),
    aes(xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 8.5),
    alpha = 0.5, fill = 'white') +
  geom_errorbarh(aes(xmin = coef-se, xmax = coef+se, 
                     y = (modelorder), height= 0, color = shape), 
                 show.legend = F)+
  geom_point(aes(x = coef ,  y = (modelorder),
                 color = shape,
                 shape = shape), 
             show.legend = T, size = 4)+
  
  facet_grid((modeltype) ~ covariate, scales = "free_y", drop = T, 
             switch = "y", space = "free")+
  geom_vline(xintercept=0, linetype = "dashed")+
  # scale_y_continuous(expand = c(0.2, 0.1)) +
  
  scale_shape_manual(name = "", 
                     values = c(20, 17, 9, 0, 15,7, 25,8),
                     labels = c('Light', 'Precipitation', 
                                'Max Interstorm', 'Width:Area', 
                                'Flashiness (RBI)', 'Terrestrial NPP', 
                                'Discharge (cv)', 'Drainage Density'
                     ))+  
  scale_color_manual(name = "", 
                     values = c("#F8766D", "#00BFC4", "#7CAE00", "#00BFC4",
                                "#7CAE00", "#00BFC4", "#7CAE00", "#00BFC4"),
                     labels = c('Light', 'Precipitation', 
                                'Max Interstorm', 'Width:Area', 
                                'Flashiness (RBI)', 'Terrestrial NPP', 
                                'Discharge (cv)', 'Drainage Density'))+  
  theme(legend.position = "bottom")+
  theme(panel.spacing.x  =unit(0.04, "lines"), panel.spacing.y = unit(0,"lines"),
        panel.border = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y=element_blank())+
  theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 30))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.text.y.left  = element_text(angle = 90),
        legend.position = "bottom")+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+ 
  theme(text = element_text(size=16))


cowplot::save_plot(here("figures/Figure3_NEP.png"), 
                   fig3.3, base_width = 9, base_height = 10, dpi = 600)
cowplot::save_plot(here("figures/Figure3_PR.png"), 
                   fig3.3, base_width = 9, base_height = 10, dpi = 600)



#Figure 4 Attempt by MD
library(ggrepel)
coefs = read.csv("data_working/sparse_quantile_regression_results_NEP_PR.csv")
fig_4 = coefs %>% 
  mutate(color = case_when(lab %in% c("Width", "Stream Surface PAR ", 
                                      "PAR kurtosis") ~"Light",
                           lab %in% c("Max interstorm", "Annual Discharge Amp") ~"Disturbance",
                           lab == "Terrestrial NPP"~"Connectivity",
                           TRUE ~ "Other")) %>% 
  mutate(lab = case_when(lab == "Stream Surface PAR "~ "Stream Surface Light",
                           lab == "PAR kurtosis" ~ "Light Kurtosis",
                           lab == "Max interstorm" ~ "Max Interstorm",
                           lab == "Annual Discharge Amp" ~ "Discharge Range",
                           TRUE ~ lab)) %>% 
  mutate(color = factor(color, levels = c("Light","Disturbance","Connectivity",NULL))) %>% 
  
  rename(PR = value_pr, NEP = value_nep, PR_se = se_pr, NEP_se = se_nep, label = lab) %>% 
  pivot_longer(cols = c(PR,NEP), names_to = "Metric") %>% 
  mutate(se = case_when(Metric == "PR"~PR_se,
                        Metric == "NEP"~NEP_se)) %>% 
  select(-c(PR_se,NEP_se)) %>% 
  group_by(Metric, color) %>% arrange(-value) %>%
  mutate(num = row_number()) %>% 
  mutate(order = case_when(label == "Width"~2,
                           label == "Stream Surface Light"~1,
                           label == "Light Kurtosis" ~3,
                           label == "Discharge Range"~4,
                           label == "Max Interstorm" ~ 5,
                           label == "Terrestrial NPP" ~ 6,
                           label == "Elevation"~7,
                           label == "Mean Air Temp" ~ 8)) %>% 
 
  ggplot(aes(y=-order, x = value))+
  geom_point(alpha = .8, size = 2, show.legend = T, aes(color = color, shape = Metric))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_errorbarh(aes(xmin = value-se, xmax = value+se, color = color), alpha = .8, height = 0, show.legend = F)+
  scale_shape_manual(values = c(19,1), guide = guide_legend(order = 1))+
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4","grey"),
                     label = c("Light","Disturbance","Connectivity","Other"),
                     guide = guide_legend(order = 2))+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+
  guides(
    shape = guide_legend(order = 1), 
    color = guide_legend(order = 2))+
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    legend.box = 'vertical',
    legend.key=element_blank(),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(-0.2,0,0,0, unit="cm"),
    panel.border = element_rect(fill = "transparent"),
    strip.background = element_blank())+
  lims(x = c(-.4,.4))+
  geom_text(aes(label = label, color = color),
            min.segment.length = .001,
            position = position_nudge(y = -.2),
            show.legend = F, check_overlap = T)+
  theme(panel.spacing = unit(0,"lines"))
fig_4

#Merging Figures 3-4


cowplot::save_plot(here("figures/Figure4_sparse_NEP_PR.png"), 
                   fig_4, base_width = 4, base_height = 5, dpi = 600)
                          

Fig_34 = ((fig3+theme(legend.position = "bottom")|fig_4+guides(color = guide_none(),shape = guide_none())+theme(legend.position = "none")))+plot_annotation(tag_levels = "A")+ 
  plot_layout(widths = c(5,3), heights = c(6,.5), guides = "collect")& theme(legend.position = 'bottom')
Fig_34
# save_plot("figures/Fig_34.png",Fig_34, dpi = 300, base_width = 9, base_height = 5)





