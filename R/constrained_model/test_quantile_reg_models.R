# Test quantile regression model for across sites 
# pooled by site as there was insufficient evidence of variation across siteyears
# A Carter

library(tidyverse)
library(quantreg)
library(Qtools)
library(corrplot)

dat <- read_csv('data_working/across_sites_model_data.csv')

dd <- dat %>%
length(unique(dd$site_name))
# filter dataframe to remove NA's in the relevant variables:
dd <- dat %>%
    select(site_name, ER = ann_ER_C, GPP = ann_GPP_C, PR, NEP = ann_NEP_C, 
           drainage_density_connected, drainage_density, Stream_PAR_sum,
           MOD_ann_NPP, PrecipWs,# precip_runoff_ratio, 
           Disch_cv, Disch_ar1, RBI, max_interstorm, ws_area_km2, 
           width_to_area) %>%
    group_by(site_name) %>%
    summarize(across(everything(), mean, na.rm = T)) %>%
    filter(!is.na(PR),
           !is.na(drainage_density_connected),
           !is.infinite(drainage_density_connected),
           !is.na(Stream_PAR_sum),
           !is.na(MOD_ann_NPP)) %>%
    mutate(ER = - ER, 
           log_PR = log(PR),
           log_RBI = log(RBI),
           # across(c(NEP, PR, log_PR), ~./max(., na.rm = T)),
           across(-c(site_name), ~scale(.)[,1]))

# look at correlation of predictors with GPP and ER, and with eachother
cor(select(dd, where(is.numeric)), use = 'complete.obs')%>%
  as.data.frame() %>%
  select(GPP, ER) %>%
  arrange(desc(abs(GPP)))

cc <- cor(select(dd, where(is.numeric)), use = 'complete.obs')
cc_sig <- cor.mtest(select(dd, where(is.numeric)), conf.level=0.95)

corrplot(cc, type = "lower", order = "original", p.mat = cc_sig$p,
         tl.col = "black", tl.srt = 90, tl.cex = .7,
         method='color', diag=FALSE, cl.cex = .6,
         cl.length = 11,
         sig.level = c(0.05), insig = 'label_sig',
         pch.cex = 1, pch.col='grey20')

# Example of 'why quantile regression' ####
dat %>%
  rename(GPP = ann_GPP_C, ER = ann_ER_C, NEP = ann_NEP_C) %>%
  pivot_longer(cols = any_of(c('GPP', 'NEP', 'ER')), 
               names_to = 'metab', values_to = 'gO2m2d') %>%
  pivot_longer(cols = any_of(c('RBI', 'width_to_area', 
                               'drainage_density_connected')), 
               names_to = 'covariate', values_to = 'value') %>%
  mutate(metab = factor(metab, levels = c('GPP', 'ER', 'NEP')))%>%
  ggplot(aes(value, gO2m2d)) +
  geom_point() +
  facet_grid(metab~covariate, scales = 'free')+
  geom_hline(yintercept = 0)+
  theme_classic()+ 
  ylab(expression(Metabolism~(gC~m^2~yr^-1)))+
  xlab('Covariate Value')+
  theme(panel.border = element_rect(fill = NA))


# test quantile regression:####
summary(lm(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv, dd))
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv,
                     tau = c(0.1, 0.5, 0.9),
                     data = dd)
summary(qmod, se = 'boot')
anova(qmod, test = 'Wald', joint = FALSE)
eps <- 0.05
kt <- KhmaladzeTest(NEP ~  Disch_cv, data = dd, 
              taus = seq(0.05, 0.95, by = 0.01),
              trim = c(eps, 1 - eps))

KhmaladzeFormat(kt, epsilon = eps)
# 
# plot(qmod)
# plot(kt)
# build table of different model fits
add_mod <- function(qmod, mod_table, D = TRUE, C = TRUE, data, y = 'NEP'){
    ss <- summary(qmod, se = 'boot')
    qm <- data.frame(mod = as.character(qmod$formula)[3],
                     tau = ss$tau,
                     L_mean = ss$coefficients[2,1],
                     L_se = ss$coefficients[2,2],
                     L_pval = ss$coefficients[2,4])
    if(D & C){
    qm <- qm %>%
      mutate(D_mean = ss$coefficients[3,1],
             D_se = ss$coefficients[3,2],
             D_pval = ss$coefficients[3,4],
             C_mean = ss$coefficients[4,1],
             C_se = ss$coefficients[4,2],
             C_pval = ss$coefficients[4,4],
             C_var = row.names(ss$coefficients)[4],
             D_var = row.names(ss$coefficients)[3])
    }
    if(D & !C){
    qm <- qm %>% 
      mutate(D_mean = ss$coefficients[3,1],
             D_se = ss$coefficients[3,2],
             D_pval = ss$coefficients[3,4],
             D_var = row.names(ss$coefficients)[3])
      
    }
    if(!D & C){
    qm <- qm %>%
      mutate(C_mean = ss$coefficients[3,1],
             C_se = ss$coefficients[3,2],
             C_pval = ss$coefficients[3,4],
             C_var = row.names(ss$coefficients)[3])
      
    }
    
    qm$AIC = AIC(qmod)
    qm$R1 = calc_R1(qmod, tau, data)
    mod_table <- bind_rows(mod_table, qm)

    return(mod_table)
}

calc_R1 <- function(qmod, tau, data, y = 'NEP'){
  rho <- function(u,tau=.5) u*(tau - (u < 0))
  qmod0 <- quantreg::rq(paste0(y,' ~ 1'),
                       tau = tau, data = data)
  
  V_tilde <- sum(rho(qmod0$resid, qmod0$tau))
  V_hat <- sum(rho(qmod$resid, qmod$tau))
  
  R1 = 1 - V_hat/V_tilde
  
  return(R1)
}

mod_table <- data.frame()
tau = 0.95
# dd$NEP2 <- dd$NEP
# dd$NEP <- dd$PR
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, C = FALSE, D = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + log_RBI,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, C = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, C = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, C = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + MOD_ann_NPP,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, D = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, D = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + drainage_density_connected,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, D = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + width_to_area,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, D = FALSE, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + log_RBI + MOD_ann_NPP,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + log_RBI + PrecipWs,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + log_RBI + drainage_density_connected,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + log_RBI + width_to_area,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + MOD_ann_NPP,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + PrecipWs,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + drainage_density_connected,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + width_to_area,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + MOD_ann_NPP,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + PrecipWs,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + drainage_density_connected,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + width_to_area,
                     tau = tau, data = dd)
mod_table <- add_mod(qmod, mod_table, data = dd)

mod_table <- mod_table %>% arrange(C_var) %>%
  mutate(C_col = case_when(C_var == 'width_to_area' ~ 'black',
                           TRUE ~ 'brown3'),
         D_col = case_when(D_var == 'log_RBI' ~ 'brown3',
                           D_var == 'Disch_cv'~'brown3',
                           TRUE ~ 'black'),
         C_pch = as.numeric(factor(C_var)) + 14,
         D_pch = as.numeric(factor(D_var)))
mod_table <- arrange(mod_table, -AIC )
mod_table
write_csv(mod_table, 'data_working/constrained_quantile_regression_results_NEP.csv')
write_csv(mod_table, 'data_working/constrained_quantile_regression_results_PR.csv')

png('figures/quantile_regression_coefficients_grouped_NEP.png', width = 7, height = 4,
    units = 'in', res = 300)
    par(mfrow = c(1,4),
        mar = c(4,0.3,5,0.3),
        oma = c(1,1,2,1))
    
    plot(mod_table$L_mean, seq(1:nrow(mod_table)), xlim = c(-0.5, 1), pch = 19,
         xlab = expression(paste('Light (', beta[1], ')')), yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$L_mean - mod_table$L_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$L_mean + mod_table$L_se, y1 = seq(1:nrow(mod_table)))
    abline(v = 0)
    mtext('Coefficient', line = 3.1, adj = 0)
    mtext('Estimates', line = 1.9, adj = 0)
    plot(mod_table$D_mean, seq(1:nrow(mod_table)), xlim = c(-.5, 1), 
         pch = mod_table$D_pch, col = mod_table$D_col, 
         xlab = expression(paste('Distrubance (', beta[2], ')')),
         yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$D_mean - mod_table$D_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$D_mean + mod_table$D_se, y1 = seq(1:nrow(mod_table)),
             col = mod_table$D_col)
    abline(v = 0)
    legend(x = -0.1, y = nrow(mod_table) + 6, 
           legend = c('max interstorm', 'CV discharge', 'RBI'),
           pch = c(3, 2, 1), col = c('black', 'brown3', 'brown3'),
           bty = 'n', xpd = TRUE)
    plot(mod_table$C_mean, seq(1:nrow(mod_table)), xlim = c(-.75, 0.75), 
         pch = mod_table$C_pch, col = mod_table$C_col, 
         xlab = expression(paste('Connectivity (', beta[3], ')')),
         yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$C_mean - mod_table$C_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$C_mean + mod_table$C_se, y1 = seq(1:nrow(mod_table)), 
             col = mod_table$C_col)
    abline(v = 0)
    legend(x = -0.5, y = nrow(mod_table) + 6, 
           legend = c('width to area', 'precipitation',
                      'terrestrial NPP', 'con. drainage dens'),
           pch = c(18, 17, 16, 15), col = c('black', rep('brown3', 3)),
           bty = 'n', xpd = TRUE)
    plot(rep(0, nrow(mod_table)), seq(1:nrow(mod_table)), ann = FALSE, type = 'n',
         axes = FALSE, xlim = c(-1,2))
    text(rep(0, nrow(mod_table)), seq(1:nrow(mod_table)),
         labels = paste(round(mod_table$AIC, 0), round(mod_table$R1, 2), sep = ",   "))
    mtext('Model', cex = 0.6, adj = 0.29, line = 1)
    mtext('AIC, R1', cex = 0.6, adj = 0.29)
dev.off()
png('figures/quantile_regression_coefficients_grouped_PR.png', width = 7, height = 4,
    units = 'in', res = 300)
    par(mfrow = c(1,4),
        mar = c(4,0.3,5,0.3),
        oma = c(1,1,2,1))
    
    plot(mod_table$L_mean, seq(1:nrow(mod_table)), xlim = c(-0.4, 1.4), pch = 19,
         xlab = expression(paste('Light (', beta[1], ')')), yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$L_mean - mod_table$L_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$L_mean + mod_table$L_se, y1 = seq(1:nrow(mod_table)))
    abline(v = 0)
    mtext('Coefficient', line = 3.1, adj = 0)
    mtext('Estimates', line = 1.9, adj = 0)
    plot(mod_table$D_mean, seq(1:nrow(mod_table)), xlim = c(-.9, 0.9), 
         pch = mod_table$D_pch, col = mod_table$D_col, 
         xlab = expression(paste('Distrubance (', beta[2], ')')),
         yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$D_mean - mod_table$D_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$D_mean + mod_table$D_se, y1 = seq(1:nrow(mod_table)),
             col = mod_table$D_col)
    abline(v = 0)
    legend(x = -0.1, y = nrow(mod_table) + 6, 
           legend = c('max interstorm', 'CV discharge', 'RBI'),
           pch = c(3, 2, 1), col = c('black', 'brown3', 'brown3'),
           bty = 'n', xpd = TRUE)
    plot(mod_table$C_mean, seq(1:nrow(mod_table)), xlim = c(-.9, 0.9), 
         pch = mod_table$C_pch, col = mod_table$C_col, 
         xlab = expression(paste('Connectivity (', beta[3], ')')),
         yaxt = 'n', bty = 'n')
    segments(x0 = mod_table$C_mean - mod_table$C_se, y0 = seq(1:nrow(mod_table)),
             x1 = mod_table$C_mean + mod_table$C_se, y1 = seq(1:nrow(mod_table)), 
             col = mod_table$C_col)
    abline(v = 0)
    legend(x = -0.5, y = nrow(mod_table) + 6, 
           legend = c('width to area', 'precipitation',
                      'terrestrial NPP', 'con. drainage dens'),
           pch = c(18, 17, 16, 15), col = c('black', rep('brown3', 3)),
           bty = 'n', xpd = TRUE)
    plot(rep(0, nrow(mod_table)), seq(1:nrow(mod_table)), ann = FALSE, type = 'n',
         axes = FALSE, xlim = c(-1,2))
    text(rep(0, nrow(mod_table)), seq(1:nrow(mod_table)),
         labels = paste(round(mod_table$AIC, 0), signif(mod_table$R1, 2), sep = ",   "))
    mtext('Model', cex = 0.6, adj = 0.29, line = 1)
    mtext('AIC, R1', cex = 0.6, adj = 0.29)
dev.off()


#### model fits:
mod_table

qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + width_to_area,
                     tau = tau, data = dd)

dd$NEP_pred <- predict(qmod, dd[,-5])
png('figures/quantile_model_predictions.png',
    width = 7, height = 3, units = 'in', res = 300)

dd %>%
  mutate(width_to_area = case_when(width_to_area > 4 ~ 3.57,
                                   TRUE ~ width_to_area),
         max_interstorm = case_when(max_interstorm > 10 ~ 2.2,
                                    TRUE ~ max_interstorm)) %>%
  select(Light = Stream_PAR_sum, `Width to Area` = width_to_area, 
         `Max interstorm interval` = max_interstorm) %>%
  mutate(across(.fns = ~ (. - min(., na.rm = T))/
                  (max(., na.rm = T)-min(., na.rm = T)))) %>% 
  bind_cols(select(dd, NEP, NEP_pred)) %>% 
  pivot_longer(cols = any_of(c('Light', 'Width to Area', 
                               'Max interstorm interval')),
               values_to = 'Value', names_to = 'covariate') %>%
  ggplot(aes(NEP, NEP_pred, col = Value))+
  geom_abline(slope = 1, intercept = 0, lty = 2, col = 'grey50')+
  geom_point() +
  scale_color_continuous(type = 'viridis') +
  ylab(expression(Predicted~95^th~'quantile\nNEP'~(g~O[2]~m^2~d^-1)))+
  xlab(expression('NEP (g '~ O[2]~ m^2~ d^-1* ')'))+
  scale_x_continuous(breaks=c(-4,-2,0,2,4), limits=c(-4,4))+
  scale_y_continuous(breaks=c(-4,-2,0,2,4), limits=c(-4,4))+
  facet_wrap(.~covariate, ncol = 3)+
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        legend.position = "bottom")
dev.off()
# individual variable quantile regressions:####
par(mfrow = c(3,3),
    mar = c(2,2,1,1),
    oma = c(1,3,1,1))

qmod <- quantreg::rq(NEP ~ Stream_PAR_sum,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$Stream_PAR_sum, dd$NEP, ylab = 'NEP', xlab = 'Light', pch = 20)
summary(qmod)
abline(-1.956, -0.0474, col = 'brown3', lty = 2)
abline(1.1, 0.34, col = 'brown3', lty = 2)
# abline(-1.935, -0.276, col = 'brown3', lty = 2)
# abline(1.08, 0.278, col = 'brown3', lty = 2)
mtext('light', 3, -1.5, adj = .1)
qmod <- quantreg::rq(NEP ~ Disch_cv,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$Disch_cv, dd$NEP, ylab = 'NEP', xlab = 'CV discharge', pch = 20)
summary(qmod)
abline(-1.68, 0.792, col = 'brown3', lty = 2)
abline(1.018, -0.1677, col = 'brown3', lty = 2)
# abline(-1.755, 0.5296, col = 'brown3', lty = 2)
# abline(1.096, -0.145, col = 'brown3', lty = 2)
mtext('CV Q', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ max_interstorm,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$max_interstorm, dd$NEP, ylab = 'NEP', xlab = 'Max interstorm', pch = 20)
summary(qmod)
abline(-1.87, -1.04, col = 'brown3', lty = 2)
abline(1.15, 0.684, col = 'brown3', lty = 2)
# abline(-2.029, -0.8756, col = 'brown3', lty = 2)
# abline(1.12, 0.277, col = 'brown3', lty = 2)
mtext('max interstorm', 3, -1.5, adj = .1)
qmod <- quantreg::rq(NEP ~ MOD_ann_NPP,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$MOD_ann_NPP, dd$NEP, ylab = 'NEP', xlab = 'terrestrial NPP', pch = 20)
summary(qmod)
abline(-1.885, -0.177, col = 'brown3', lty = 2)
abline(1.149, -0.126, col = 'brown3', lty = 2)
# abline(-1.97, -0.1207, col = 'brown3', lty = 2)
# abline(1.108, -0.091, col = 'brown3', lty = 2)
mtext('NEP',2, 2)
mtext('terr NPP', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ PrecipWs,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$PrecipWs, dd$NEP, ylab = 'NEP', xlab = 'precipitation', pch = 20)
summary(qmod)
abline(-1.99, -0.02, col = 'brown3', lty = 2)
abline(1.096, -0.208, col = 'brown3', lty = 2)
# abline(-2.077, 0.353, col = 'brown3', lty = 2)
# abline(1.102, -0.055, col = 'brown3', lty = 2)
mtext('precip', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ drainage_density_connected,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$drainage_density_connected, dd$NEP, ylab = 'NEP', xlab = 'connected drainage density', pch = 20)
summary(qmod)
abline(-1.505, 0.532, col = 'brown3', lty = 2)
abline(1.164, -0.222, col = 'brown3', lty = 2)
# abline(-1.678, 0.806, col = 'brown3', lty = 2)
# abline(1.118, -0.0037, col = 'brown3', lty = 2)
mtext('DD con', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ width_to_area,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$width_to_area, dd$NEP, ylab = 'NEP', xlab = 'width to area', pch = 20)
summary(qmod)
abline(-1.95, 0.57, col = 'brown3', lty = 2)
abline(0.91, .117, col = 'brown3', lty = 2)

# abline(-1.816, 1.83, col = 'brown3', lty = 2)
# abline(1.124, 0.16, col = 'brown3', lty = 2)
mtext('W/A', 3, -1.5, adj = .9)

