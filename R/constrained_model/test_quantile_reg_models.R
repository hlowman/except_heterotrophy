# Test quantile regression model for across sites 
# pooled by site as there was insufficient evidence of variation across siteyears
# A Carter

library(tidyverse)
library(quantreg)
library(corrplot)

dat <- read_csv('data_working/across_sites_model_data.csv')

dat$PR = -dat$ann_GPP_C/dat$ann_ER_C
dat$NEP = dat$ann_GPP_C + dat$ann_ER_C

dd <- dat %>%
    filter(!is.na(PR),
           !is.na(drainage_density_connected),
           !is.infinite(drainage_density_connected),
           !is.na(Stream_PAR_sum), 
           !is.na(MOD_ann_NPP)) %>%
    mutate(width_to_area = Width/sqrt(ws_area_km2))%>%
    select(site_name, ER = ann_ER_C, GPP = ann_GPP_C, PR, NEP, 
           drainage_density_connected, drainage_density, Stream_PAR_sum,
           MOD_ann_NPP, PrecipWs,# precip_runoff_ratio, 
           Disch_cv, max_interstorm, ws_area_km2, 
           width_to_area) %>%
    mutate(ER = - ER, 
           log_PR = log(PR),
           # across(c(NEP, PR, log_PR), ~./max(., na.rm = T)),
           across(-c(site_name), ~scale(.)[,1]))

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

# test quantile regression:
summary(lm(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv, dd))
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv,
                     tau = c(0.05, 0.5, 0.95),
                     data = dd)
summary(qmod, se = 'boot')
anova(qmod, test = 'Wald', joint = FALSE)

plot(qmod)

# build table of different model fits
add_mod <- function(qmod, mod_table, D = TRUE, C = TRUE){
    ss <- summary(qmod, se = 'boot')
    if(D & C){
    qm <- data.frame(mod = as.character(qmod$formula)[3],
                     tau = ss$tau,
                     L_mean = ss$coefficients[2,1],
                     L_se = ss$coefficients[2,2],
                     D_mean = ss$coefficients[3,1],
                     D_se = ss$coefficients[3,2],
                     C_mean = ss$coefficients[4,1],
                     C_se = ss$coefficients[4,2],
                     C_var = row.names(ss$coefficients)[4],
                     D_var = row.names(ss$coefficients)[3])
    }
    if(D & !C){
    qm <- data.frame(mod = as.character(qmod$formula)[3],
                     tau = ss$tau,
                     L_mean = ss$coefficients[2,1],
                     L_se = ss$coefficients[2,2],
                     D_mean = ss$coefficients[3,1],
                     D_se = ss$coefficients[3,2],
                     D_var = row.names(ss$coefficients)[3])
      
    }
    if(!D & C){
    qm <- data.frame(mod = as.character(qmod$formula)[3],
                     tau = ss$tau,
                     L_mean = ss$coefficients[2,1],
                     L_se = ss$coefficients[2,2],
                     C_mean = ss$coefficients[3,1],
                     C_se = ss$coefficients[3,2],
                     C_var = row.names(ss$coefficients)[3])
      
    }
      
    mod_table <- bind_rows(mod_table, qm)

    return(mod_table)
}

mod_table <- data.frame()

# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, C = FALSE)
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, C = FALSE)
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + MOD_ann_NPP,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, D = FALSE)
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, D = FALSE)
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + drainage_density_connected,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, D = FALSE)
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + width_to_area,
#                      tau = 0.95, data = dd)
# mod_table <- add_mod(qmod, mod_table, D = FALSE)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + MOD_ann_NPP,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + PrecipWs,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + drainage_density_connected,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + Disch_cv + width_to_area,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + MOD_ann_NPP,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + PrecipWs,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + drainage_density_connected,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)
qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + max_interstorm + width_to_area,
                     tau = 0.95, data = dd)
mod_table <- add_mod(qmod, mod_table)

mod_table <- mod_table %>% arrange(C_var) %>%
  mutate(C_col = case_when(C_var == 'width_to_area' ~ 'black',
                           TRUE ~ 'brown3'),
         D_col = case_when(D_var == 'Disch_cv' ~ 'brown3',
                           TRUE ~ 'black'),
         C_pch = as.numeric(factor(C_var)) + 14,
         D_pch = as.numeric(factor(D_var)))

par(mfrow = c(1,3),
    mar = c(4,0.3,5,0.3),
    oma = c(1,1,1,1))

plot(mod_table$L_mean, seq(1:nrow(mod_table)), xlim = c(-.2, .4), pch = 19,
     xlab = 'Light', yaxt = 'n', bty = 'n')
segments(x0 = mod_table$L_mean - mod_table$L_se, y0 = seq(1:nrow(mod_table)),
         x1 = mod_table$L_mean + mod_table$L_se, y1 = seq(1:nrow(mod_table)))
abline(v = 0)
mtext('Coefficient', line = 3.1, adj = 0)
mtext('Estimates', line = 1.9, adj = 0)
plot(mod_table$D_mean, seq(1:nrow(mod_table)), xlim = c(-.2, .4), 
     pch = mod_table$D_pch, col = mod_table$D_col, xlab = 'Distrubance',
     yaxt = 'n', bty = 'n')
segments(x0 = mod_table$D_mean - mod_table$D_se, y0 = seq(1:nrow(mod_table)),
         x1 = mod_table$D_mean + mod_table$D_se, y1 = seq(1:nrow(mod_table)),
         col = mod_table$D_col)
abline(v = 0)
legend(x = -0.2, y = nrow(mod_table) + 2, 
       legend = c('max interstorm', 'CV discharge'),
       pch = c(2, 1), col = c('black', 'brown3'),
       bty = 'n', xpd = TRUE)
plot(mod_table$C_mean, seq(1:nrow(mod_table)), xlim = c(-.2, .4), 
     pch = mod_table$C_pch, col = mod_table$C_col, xlab = 'Connectivity',
     yaxt = 'n', bty = 'n')
segments(x0 = mod_table$C_mean - mod_table$C_se, y0 = seq(1:nrow(mod_table)),
         x1 = mod_table$C_mean + mod_table$C_se, y1 = seq(1:nrow(mod_table)), 
         col = mod_table$C_col)
abline(v = 0)
legend(x = -.2, y = nrow(mod_table) + 2.5, 
       legend = c('width to area', 'precipitation',
                  'terrestrial NPP', 'con. drainage dens'),
       pch = c(18, 17, 16, 15), col = c('black', rep('brown3', 3)),
       bty = 'n', xpd = TRUE)

# individual variable quantile regressions:
par(mfrow = c(3,3), 
    mar = c(2,2,1,1),
    oma = c(1,3,1,1))

qmod <- quantreg::rq(NEP ~ Stream_PAR_sum,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$Stream_PAR_sum, dd$NEP, ylab = 'NEP', xlab = 'Light', pch = 20)
summary(qmod)
abline(-1.935, -0.276, col = 'brown3', lty = 2)
abline(1.08, 0.278, col = 'brown3', lty = 2)
mtext('light', 3, -1.5, adj = .1)
qmod <- quantreg::rq(NEP ~ Disch_cv,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$Disch_cv, dd$NEP, ylab = 'NEP', xlab = 'CV discharge', pch = 20)
summary(qmod)
abline(-1.755, 0.5296, col = 'brown3', lty = 2)
abline(1.096, -0.145, col = 'brown3', lty = 2)
mtext('CV Q', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ max_interstorm,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$max_interstorm, dd$NEP, ylab = 'NEP', xlab = 'Max interstorm', pch = 20)
summary(qmod)
abline(-2.029, -0.8756, col = 'brown3', lty = 2)
abline(1.12, 0.277, col = 'brown3', lty = 2)
mtext('max interstorm', 3, -1.5, adj = .1)
qmod <- quantreg::rq(NEP ~ MOD_ann_NPP,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$MOD_ann_NPP, dd$NEP, ylab = 'NEP', xlab = 'terrestrial NPP', pch = 20)
summary(qmod)
abline(-1.97, -0.1207, col = 'brown3', lty = 2)
abline(1.108, -0.091, col = 'brown3', lty = 2)
mtext('NEP',2, 2)
mtext('terr NPP', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ PrecipWs,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$PrecipWs, dd$NEP, ylab = 'NEP', xlab = 'precipitation', pch = 20)
summary(qmod)
abline(-2.077, 0.353, col = 'brown3', lty = 2)
abline(1.102, -0.055, col = 'brown3', lty = 2)
mtext('precip', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ drainage_density_connected,
                     tau = c(0.05, 0.95), data = dd)
plot(dat$drainage_density_connected, dat$NEP, ylab = 'NEP', xlab = 'connected drainage density', pch = 20)
summary(qmod)
abline(-1.678, 0.806, col = 'brown3', lty = 2)
abline(1.118, -0.0037, col = 'brown3', lty = 2)
mtext('DD con', 3, -1.5, adj = .9)
qmod <- quantreg::rq(NEP ~ width_to_area,
                     tau = c(0.05, 0.95), data = dd)
plot(dd$width_to_area, dd$NEP, ylab = 'NEP', xlab = 'width to area', pch = 20)
summary(qmod)
abline(-1.816, 1.83, col = 'brown3', lty = 2)
abline(1.124, 0.16, col = 'brown3', lty = 2)
mtext('W/A', 3, -1.5, adj = .9)

entify(dat$drainage_density_connected, dat$ann_ER_C, labels = dat$site_name)
plot(dd$drainage_density_connected, dd$NEP)
plot(dd$width_to_area, dd$NEP, xlim = c(-0.51, 2))
plot(dd$max_interstorm, dd$GPP)
mods <- data.frame()
add_model_to_table <- function(b, mods){
  bb <- data.frame(fixef(b)) 
  bb$variable = rownames(bb)
  rownames(bb) <- NULL
  
  bb$mod <- as.character(b$formula)[1]
  bb$waic <- brms::waic(b)$waic
  
  mods <- bind_rows(mods, bb)
  return(mods)
}


b <- brm(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected +
           Disch_cv + (1|site_name), data = dd)
post <- brms::posterior_samples(b,
              pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'drainage_density_connected', 'Disch_cv'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_drainage_density_connected', 'b_Disch_cv')))
p1 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

mods <- add_model_to_table(b, mods)

b2 <- brm(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected +
           max_interstorm + (1|site_name), data = dd)

post <- brms::posterior_samples(b2,
              pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'drainage_density_connected', 'max_interstorm'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_drainage_density_connected', 'b_max_interstorm')))
p2 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

mods <- add_model_to_table(b2, mods)

b3 <- brm(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + width_to_area +
           max_interstorm + (1|site_name), data = dd)

post <- brms::posterior_samples(b3,
              pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'max_interstorm'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_width_to_area', 'b_max_interstorm')))
p3 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())


mods <- add_model_to_table(b3, mods)

b4 <- brm(log_PR ~ Stream_PAR_sum + MOD_ann_NPP + width_to_area +
           Disch_cv + (1|site_name), data = dd)

post <- brms::posterior_samples(b4,
              pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'Disch_cv'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_width_to_area', 'b_Disch_cv')))
p4 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())



mods <- add_model_to_table(b4, mods)

mods <- read_csv('data_working/brms_models.csv')

ggpubr::ggarrange(p1,p2,p3,p4, ncol = 4)

#NEP
b <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected +
           Disch_cv + (1|site_name), data = dd)

post <- brms::posterior_samples(b,
                                pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'drainage_density_connected', 'Disch_cv'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_drainage_density_connected', 'b_Disch_cv')))
p1 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())



# mods <- add_model_to_table(b, mods)

b2 <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected +
           max_interstorm + (1|site_name), data = dd)
post <- brms::posterior_samples(b2,
                                pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'drainage_density_connected', 'max_interstorm'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_drainage_density_connected', 'b_max_interstorm')))
p2 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# mods <- add_model_to_table(b2, mods)

b3 <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + width_to_area +
           max_interstorm + (1|site_name), data = dd)
post <- brms::posterior_samples(b3,
                                pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'max_interstorm'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_width_to_area', 'b_max_interstorm')))
p3 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0)+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# mods <- add_model_to_table(b3, mods)

b4 <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + width_to_area +
           Disch_cv + (1|site_name), data = dd)


post <- brms::posterior_samples(b4,
                                pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'Disch_cv'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples') %>%
  mutate(parameter = factor(parameter, levels = c('b_Stream_PAR_sum', 'b_MOD_ann_NPP', 'b_width_to_area', 'b_Disch_cv')))
p4 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.5) +
  geom_hline(yintercept = 0) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())


ggpubr::ggarrange(p1,p2,p3,p4, ncol = 4)


post <- brms::posterior_samples(b4,
              pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'Disch_cv'))
pp <- post %>% pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'samples')
p4 <- ggplot(pp, aes(x=parameter, y=samples)) + 
  geom_boxplot(fill='steelblue') +
  coord_flip() +
  ylim(-0.3, 0.3)

mods <- add_model_to_table(b4, mods)

write_csv(mods, 'data_working/brms_models.csv')
plot(b4, pars = c('Stream_PAR_sum', 'MOD_ann_NPP', 'width_to_area', 'Disch_cv'))
p[[1]] + xlim(-.5, .5)
  
  plot(conditional_effects(b), points = TRUE)
pp_check(b)
b_NEP <- brm(NEP ~ Stream_PAR_sum + MOD_ann_NPP + drainage_density_connected + 
      (1|site_name), data = dd)

l <- lmer(GPP ~ Stream_PAR_sum + MOD_ann_NPP + Disch_cv +
            drainage_density_connected + (1|site_name), data = dd)

l_nep <- lmer(NEP ~ Stream_PAR_sum + MOD_ann_NPP + 
            drainage_density_connected + (1|site_name), data = dd)

summary(l)


sum(is.na(dat$MOD_ann_NPP))
