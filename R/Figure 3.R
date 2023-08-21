

mod_tableNEP<-read.csv('data_working/constrained_quantile_regression_results_NEP.csv')





mod_tableNEP<-mod_tableNEP %>% mutate(modeltype= case_when((!is.na(D_mean) & !is.na(C_mean) ~ "Light + Dist. + Conn."), 
                                                           (!is.na(C_mean) & is.na(D_mean) ~ "Light + Conn."),
                                                           (!is.na(D_mean) & is.na(C_mean) ~ "Light + Dist."),
                                                           (is.na(D_mean) & is.na(C_mean)~ "Light")))


mod_tableNEP$modeltype<-factor(mod_tableNEP$modeltype,levels=c("Light", "Light + Dist.", "Light + Conn.", "Light + Dist. + Conn."))


mod_tableNEP<-mod_tableNEP %>% arrange(modeltype) %>% mutate(modelorder=seq(1:20))

mod_tableNEPcoeflong <- melt(mod_tableNEP,
                  id.vars=c("mod", "modelorder", "modeltype", "D_var", "C_var"),
                  measure.vars=c("L_mean", "D_mean", "C_mean" ),
                  variable.name="covariate",
                  value.name="coef")


mod_tableNEPerrorlong <- melt(mod_tableNEP,
                             id.vars=c("mod", "modelorder", "modeltype", "D_var", "C_var"),
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

fig3<-ggplot(data = mod_tablelong2) +
  geom_point(aes(x = coef , y = as.factor(desc(modelorder)))) +
  facet_wrap(~covariate)+
  geom_vline(xintercept=0)


                          