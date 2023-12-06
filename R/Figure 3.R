library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(cowplot)


mod_tableNEP<-read.csv('data_working/constrained_quantile_regression_results_NEP.csv')





mod_tableNEP<-mod_tableNEP %>% mutate(modeltype= case_when((!is.na(D_mean) & !is.na(C_mean) ~ "Light + Dist. + Conn."), 
                                                           (!is.na(C_mean) & is.na(D_mean) ~ "Light + Conn."),
                                                           (!is.na(D_mean) & is.na(C_mean) ~ "Light + Dist."),
                                                           (is.na(D_mean) & is.na(C_mean)~ "Light")))


mod_tableNEP$modeltype<-factor(mod_tableNEP$modeltype,levels=c("Light", "Light + Dist.", "Light + Conn.", "Light + Dist. + Conn."))


mod_tableNEP<-mod_tableNEP %>% arrange(modeltype) %>% mutate(modelorder=seq(1:20))

mod_tableNEPcoeflong <- melt(mod_tableNEP,
                  id.vars=c("mod", "modelorder", "modeltype", "D_var", "C_var", "AIC"),
                  measure.vars=c("L_mean", "D_mean", "C_mean" ),
                  variable.name="covariate",
                  value.name="coef")


mod_tableNEPerrorlong <- melt(mod_tableNEP,
                             id.vars=c("mod", "modelorder", "modeltype", "D_var", "C_var", "AIC"),
                             measure.vars=c("L_se", "D_se", "C_se"),
                             variable.name="covariate2",
                             value.name="error")

mod_tablelong<-merge(mod_tableNEPcoeflong,mod_tableNEPerrorlong) 

mod_tablelong2<-mod_tablelong %>% drop_na(error)


mod_tablelong2$covariate <- recode(mod_tablelong2$covariate, C_mean = 'Connectivity', 
                          D_mean = 'Disturbance',
                          L_mean = 'Light')

mod_tablelong2$covariate<-factor(mod_tablelong2$covariate,levels=c("Light", "Disturbance", "Connectivity"))


#### making the figure ####

fig3<-mod_tablelong2 %>% as_tibble() %>% mutate(dAIC = as.numeric(round(AIC-min(AIC), 1))) %>% 
  mutate(modeltype = as.factor(modeltype)) %>% 
  mutate(Metric = case_when(is.na(coef)==T~"PR",T~"NEP")) %>% ungroup() %>% 
  ggplot() + 
  geom_rect(aes(fill = modeltype,alpha = modeltype),xmin = -Inf,xmax = Inf,
                       ymin = -Inf,ymax = Inf, show.legend = F) +
  scale_alpha_manual(values = c(.25,.05,.05,.01))+
  geom_point(aes(x = coef ,  y = as.factor(desc(modelorder)),size = -dAIC, shape =Metric), show.legend = T)+
  geom_errorbarh(aes(xmin = coef-error, xmax = coef+error, y = as.factor(desc(modelorder)),height = 0, color=modeltype), show.legend = F)+
  facet_grid(modeltype~covariate, scales = "free_y", drop = T, switch = "y", space = "free")+
  geom_vline(xintercept=0, linetype = "dashed")+
  
  geom_point(aes(x = coef , y = as.factor(desc(modelorder)),  size = -(dAIC), color=modeltype), show.legend = F)+
  scale_shape_manual(values = c(19,1))+
  scale_size_continuous(range = c(.5,8), guide = T)+
  theme(legend.position = "bottom")+
  theme(panel.spacing.x  =unit(0.04, "lines"), panel.spacing.y = unit(0,"lines"),
        panel.border = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y=element_blank())+
  theme_bw()+theme(panel.grid=element_blank(),strip.background = element_blank())+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y = element_blank(),strip.text.y.left  = element_text(angle = 0, hjust = 1),
        legend.position = "bottom")+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+ 
  guides(size = guide_legend(expression(paste(delta, 'AIC')), direction = "horizontal", nrow = 1, reverse = T),
    shape = guide_legend(override.aes = list(size = 5), title = "Metric"), color = guide_none(), fill = guide_none(), alpha = guide_none())

  
  
fig3


## Updated Figure 3 ##
#Figure 3 changes 11/02/2023 (MD) and 12/5/2023 (CT)
fig3.2 <- rbind(mod_tablelong2 %>% as_tibble() %>% mutate(dAIC = as.numeric(round(AIC-min(AIC), 1))) %>% 
                  mutate(modeltype = as.factor(modeltype)) %>% 
                  mutate(Metric = case_when(is.na(coef)==T~"PR",T~"NEP")) %>% ungroup() %>% 
                  mutate(L_var = "Light") %>% 
                  # selecting the lowest 27% of AIC values overall: NOTE: Grouping by modeltype did not work if there were <4 models in a modeltype (when prop was = 0.25)
                  slice_min(order_by=AIC, prop = .27) %>% mutate(best = "Best"),    # changed prop to 0.27 and 0.73 to keep both models with AICd = 9.4 in the "Best" category
                
                mod_tablelong2 %>% as_tibble() %>% mutate(dAIC = as.numeric(round(AIC-min(AIC), 1))) %>% 
                  mutate(modeltype = as.factor(modeltype)) %>% 
                  mutate(Metric = case_when(is.na(coef)==T~"PR",T~"NEP")) %>% ungroup() %>% 
                  mutate(L_var = "Light") %>% 
                  # selecting the highest 73% of AIC values overall. See above re: grouping by modeltype
                  slice_max(order_by=AIC, prop = .73) %>% mutate(best = "Other")) %>% # changed this prop to 0.73 to keep both models with AICd = 9.4 in the "Best" category
  
  filter(is.na(coef) == F) %>% 
  mutate(se = case_when(covariate == "Light" & covariate2 == "L_se"|
                          covariate == "Disturbance" & covariate2 == "D_se"|
                          covariate == "Connectivity" & covariate2 == "C_se"~ error, T~NA)) %>% 
  mutate(D_var = case_when(D_var == "Disch_cv"~ "Discharge (cv)",
                           D_var == "log_RBI" ~ "Flashiness (RBI)",
                           D_var == "max_interstorm" ~ "Max Interstorm."),
         C_var = case_when(C_var == "drainage_density_connected"~"Drainage Density",
                           C_var == "width_to_area"~"Width:Area",
                           C_var == "MOD_ann_NPP"~ "Terrestrial NPP",
                           C_var == "PrecipWs" ~ "Precipitation")) %>% 
  
  ggplot() + 
  geom_rect(aes(fill = modeltype,alpha = modeltype),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, show.legend = F) +
  scale_alpha_manual(values = c(.25,.05,.05,.01))+
  
  
  geom_errorbarh(aes(xmin = coef-se, xmax = coef+se, y = as.factor(desc(modelorder)), height= 0, color = best), show.legend = F)+
  
  
  geom_point(aes(x = coef ,  y = as.factor(desc(modelorder)),
                 color = best,
                 shape =case_when(covariate =="Light"~L_var,
                                  covariate == "Connectivity"~C_var,
                                  covariate == "Disturbance" ~ D_var)), show.legend = T, size = 4)+
  
  
  
  facet_grid((modeltype)~covariate, scales = "free_y", drop = T, switch = "y", space = "free")+
  geom_vline(xintercept=0, linetype = "dashed")+
  
  scale_shape_manual(values = c(15,13,17,8,19,7,10,9,10))+  
  scale_color_manual(values = c("black","grey70"))+
  
  theme(legend.position = "bottom")+
  theme(panel.spacing.x  =unit(0.04, "lines"), panel.spacing.y = unit(0,"lines"),
        panel.border = element_blank())+
  theme(axis.text.y = element_blank(), axis.title.y=element_blank())+
  theme_bw()+theme(panel.grid=element_blank(),strip.background = element_blank(),panel.spacing = unit(1.5, "lines"),plot.margin = margin(t = 10, r = 20, b = 10, l = 30))+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y = element_blank(),strip.text.y.left  = element_text(angle = 90),
        legend.position = "bottom")+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+ 
  guides(size = guide_legend(expression(paste(delta, 'AIC')), direction = "vertical", nrow = 1, reverse = T),
         shape = guide_legend(override.aes = list(size = 5), title = "Metric"), color = guide_none(), fill = guide_none(), alpha = guide_none())+
  theme(text = element_text(size=16))

fig3.2 

cowplot::save_plot("Figure3_UpdatedDec5_SinglePanel_WithCovariates.png", fig3.2, base_width = 9, base_height = 10, dpi = 600)



#Figure 4 Attempt by MD
library(ggrepel)
coefs = read.csv("data_working/sparse_quantile_regression_results_NEP_PR.csv")
fig_4 = coefs %>% mutate(color = case_when(lab %in% c("Width","Stream Surface PAR ","PAR kurtosis")~"Light",
                                           lab %in% c("Max interstorm", "Annual Discharge Amp")~"Disturbance",
                                           lab == "Elevation"~"Connectivity",
                                           T~"Other")) %>% 
  mutate(color = factor(color, levels = c("Light","Disturbance","Connectivity",NULL))) %>% 
  
  rename(PR = value_pr, NEP = value_nep, PR_se = se_pr, NEP_se = se_nep, label = lab) %>% 
  pivot_longer(cols = c(PR,NEP), names_to = "Metric") %>% 
  mutate(se = case_when(Metric == "PR"~PR_se,
                        Metric == "NEP"~NEP_se)) %>% 
  select(-c(PR_se,NEP_se)) %>% 
  group_by(Metric, color) %>% arrange(-value) %>%
  mutate(num = row_number()) %>% 
  mutate(order = case_when(label == "Width"~2,
                           label == "Stream Surface PAR "~1,
                           label == "PAR kurtosis" ~3,
                           label == "Annual Discharge Amp"~4,
                           label == "Max interstorm" ~ 5,
                           label == "Elevation"~6,
                           label == "Terrestrial NPP" ~ 7,
                           label == "Mean Air Temp" ~ 8)) %>% 
  mutate(col = case_when(col == "black"~"Positive Expected",
                         col == "red"~"Negative Expected",
                         col == "grey"~"No Hyp.")) %>% 
  
  
  
  
  ggplot(aes(y=-order, x = value))+
  geom_point(alpha = .8, size = 2, show.legend = T, aes(color = color, shape = Metric))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_errorbarh(aes(xmin = value-se, xmax = value+se, color = color), alpha = .8, height = 0, show.legend = F)+
  scale_shape_manual(values = c(19,1))+
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4","grey"),
                     label = c("","","",""))+
  labs(x = expression(paste( 'Coefficient estimate (', beta, ')')))+
  theme(#axis.text.x =element_text(angle = 45, hjust =1), 
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank())+
  lims(x = c(-.4,.4))+
  theme(panel.border = element_rect(fill = "transparent"))+
  #facet_wrap(Metric~., scales = "free",drop = T)+
  geom_text(aes(label = label, color = color),
            min.segment.length = .001,
            position = position_nudge(y = -.2),
            show.legend = F, check_overlap = T)+
  theme(panel.spacing = unit(0,"lines"))
fig_4

#Merging Figures 3-4


                          

Fig_34 = ((fig3+theme(legend.position = "bottom")|fig_4+guides(color = guide_none(),shape = guide_none())+theme(legend.position = "none")))+plot_annotation(tag_levels = "A")+ 
  plot_layout(widths = c(5,3), heights = c(6,.5), guides = "collect")& theme(legend.position = 'bottom')
Fig_34
# save_plot("figures/Fig_34.png",Fig_34, dpi = 300, base_width = 9, base_height = 5)





