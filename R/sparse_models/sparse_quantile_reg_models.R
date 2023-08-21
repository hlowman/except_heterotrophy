# This code is adapted from the constrained quantile regression model 
#   It refits the quantile regressions using a built in lasso method and 
#   including the whole dataset.

# A Carter
library(tidyverse)
library(quantreg)
library(corrplot)
library(viridis)

dat <- read_csv('data_working/across_sites_model_data.csv')

# filter dataframe to remove NA's in the relevant variables:
apply(dat, 2, function(x) sum(is.na(x)))

dd <- dat %>%
    filter(!is.na(site_name))%>%
    group_by(site_name) %>%
    summarize(across(where(is.numeric), median, na.rm = T))

# remove % cover geologic variables as these are highly zero inflated and likely
# not relevant predictors. Log transform highly skewed variables
geol_vars <- grep('^Pct(?!_)', colnames(dd), perl = TRUE)

dd <- dd[,-geol_vars] %>%
    select(-reach_proportion) %>%
    rename(NEP = ann_NEP_C)%>%
    mutate(PR = (-ann_GPP_C/ann_ER_C)) %>%
    filter(!is.na(drainage_density),
           !is.na(PrecipWs),
           !is.na(Stream_PAR_sum),
           !is.na(width_to_area)) %>%
    # select(-ann_GPP_C, -ann_ER_C, -PR, -NHD_TIDAL) %>%
    select(-ann_GPP_C, -ann_ER_C, -NHD_TIDAL) %>%
    mutate(across(c(Disch_mean, Width, NHD_SLOPE, ws_area_km2, ElevWs,
                    starts_with(c('Pct', 'NLCD_Pct', 'Dam')), OmWs, HydrlCondWs,
                    NWs, CaOWs, P2O5Ws, SWs, precip_runoff_ratio,
                    Inorg_N_fert_kgNhayr, Org_N_fert_kgNhayr, PopDen2010Ws,
                    Waste_point_srcs_perkm2, connected_flow_length, 
                    total_flow_length), ~log(.+1)))%>%
    mutate(across(c(-site_name, -year), ~scale(.)[,1]))

corr <- cor(dd[,c(3,4,7:ncol(dd))])
png('figures/annual_corrplot.png', width = 10, height = 10, units = 'in', res = 100)
  corrplot::corrplot(corr, method = 'color', type = 'lower',  diag = FALSE,
                     order = 'hclust', tl.col = 'black', tl.cex = 0.6)
dev.off()
# hist2 <- function(x,y){
#   hist(x, main = y)
# }
# par(mar = c(2,2,1,1),
#     mfrow = c(5,5))
# tmp <- select(dd, where(is.numeric)) 
# mapply(hist2, x = tmp, y = colnames(tmp))
# test sparse quantile regression:####
# summary(lm(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv, dd))

calc_mod_stats <- function(qmod, tau, data, y = 'NEP'){
  
  AIC = AIC(qmod)
  R1 = calc_R1(qmod, tau, data, y)

  return(list(AIC = AIC, R1 = R1))
}

calc_R1 <- function(qmod, tau, data, y = 'NEP'){
  rho <- function(u,tau=.5) u*(tau - (u < 0))
  qmod0 <- quantreg::rq(PR ~ 1,
                        tau = tau, data = data)
  
  V_tilde <- sum(rho(qmod0$resid, qmod0$tau))
  V_hat <- sum(rho(qmod$resid, qmod$tau))
  
  R1 = 1 - V_hat/V_tilde
  
  return(R1)
}



# generate formula that includes all covariates:
fmlqmod_PR <- as.formula(
  paste('PR ~', paste(colnames(dd[,c(3,4,7:ncol(dd))]), collapse = '+'), sep = ' '))
fmlqmod_NEP <- as.formula(
  paste('NEP ~', paste(colnames(dd[,c(3,4,7:ncol(dd))]), collapse = '+'), sep = ' '))

# regular quantile regression
qmod <- quantreg::rq(fmlqmod_PR,
                     tau =  0.95, 
                     data = dd)

calc_mod_stats(qmod, tau = 0.95, dd, y = 'PR')

qmod <- quantreg::rq(fmlqmod_NEP,
                     tau =  0.95, 
                     data = dd)

calc_mod_stats(qmod, tau = 0.95, dd, y = 'NEP')

# Sparse quantile regression
lqmod_NEP <- quantreg::rqss(fmlqmod_NEP,
                        tau =  0.95, method = 'lasso', lambda = 10,
                        data = dd)
lqmod_PR <- quantreg::rqss(fmlqmod_PR,
                        tau =  0.95, method = 'lasso', lambda = 10,
                        data = dd)

calc_mod_stats(lqmod_NEP, tau = 0.95, dd, 'NEP')
calc_mod_stats(lqmod_PR, tau = 0.95, dd, 'PR')
summary(lqmod_PR, se = 'boot')
summary(lqmod_NEP, se = 'boot')
qq_PR <- summary(lqmod_PR, se = 'boot')
qq_NEP <- summary(lqmod_NEP, se = 'boot')
dd$NEPpred_NEP <- predict(lqmod_NEP, dd[,c(-5, -6)])
dd$NEPpred_PR <- predict(lqmod_PR, dd[,c(-5, -6)])

# model with just the significant variables:
# qmod <- quantreg::rq(NEP ~ MOD_ann_NPP + Width + max_interstorm + PAR_kurt + ElevWs,
#                      tau = 0.95, data = dd)
qmod <- quantreg::rq(PR ~ Disch_amp + Width + max_interstorm + PAR_kurt + ElevWs + Stream_PAR_sum + TmeanWs,
                     tau = 0.95, data = dd)
R1_PR <- calc_mod_stats(qmod, tau = 0.95, dd, y = 'PR')$R1
qmod <- quantreg::rq(NEP ~ MOD_ann_NPP + Width + max_interstorm + PAR_kurt + ElevWs,
                     tau = 0.95, data = dd)
R1_NEP <- calc_mod_stats(qmod, tau = 0.95, dd, y = 'NEP')$R1

# # One attempt at showing model fits:
# png('figures/sparse_quantile_model_predictions.png',
#     width = 6, height = 4.5, units = 'in', res = 300)
#   dd %>%
#     select(Width, Elevation = ElevWs, `Terrestrial NPP` = MOD_ann_NPP, 
#            `Max interstorm interval` = max_interstorm) %>%
#     mutate(across(.fns = ~ (. - min(.))/(max(.)-min(.)))) %>% 
#     bind_cols(select(dd, NEP, NEPpred)) %>% 
#     pivot_longer(cols = any_of(c('Width', 'Elevation', 'Terrestrial NPP', 
#                           'Max interstorm interval')),
#                  values_to = 'Value', names_to = 'covariate') %>%
#   ggplot(aes(NEP, NEPpred, col = Value))+
#     geom_abline(slope = 1, intercept = 0, lty = 2, col = 'grey50')+
#     geom_point() +
#     scale_color_continuous(type = 'viridis',
#                            name= "Value of each predictor")+
#     ylab(expression('Predicted NEP (g '~ O[2]~ m^2~ d^-1* ')'))+
#     xlab(expression('NEP (g '~ O[2]~ m^2~ d^-1* ')'))+
#     scale_x_continuous(breaks=c(-4,-2,0,2,4), limits=c(-4,4))+
#     scale_y_continuous(breaks=c(-4,-2,0,2,4), limits=c(-4,4))+
#     facet_wrap(.~covariate)+
#     theme_classic() +
#     theme(panel.border = element_rect(fill = NA),
#           legend.position = "bottom")
# dev.off()

colnames(qq_NEP$coef)<- c('value', 'se', 't_val', 'p_val')
coefs_NEP <- data.frame(qq_NEP$coef) %>%
  arrange(-(abs(value))) %>%
  filter(p_val < 0.99) %>% # don't keep effects with p values of 1
  slice(-1) # remove the intercept from dataframe
colnames(qq_PR$coef)<- c('value', 'se', 't_val', 'p_val')
coefs_PR <- data.frame(qq_PR$coef) %>%
  arrange(-(abs(value))) %>%
  filter(p_val < 0.99) %>% # don't keep effects with p values of 1
  slice(-1) # remove the intercept from dataframe

coefs_NEP$var = row.names(coefs_NEP)
coefs_PR$var = row.names(coefs_PR)

coefs <- full_join(coefs_PR, coefs_NEP, by = 'var') #%>%
  # select(-starts_with(c('t_val', 'p_val')))
vars <- data.frame(var = c('MOD_ann_NPP', 'TmeanWs', 'PAR_kurt',
                           'Stream_PAR_sum', 'ElevWs', 'Disch_amp',
                           'max_interstorm', 'Width'),
                   lab = c("Terrestrial NPP", "Mean Air Temp", "PAR kurtosis", 
                           "Stream Surface PAR ", "Elevation", "Annual Discharge Amp",
                           "Max interstorm", "Width"),
                   col = c('red', 'grey', 'grey', 'black', 'grey', 'red', 
                           'black', 'black'))

coefs <- full_join(vars, coefs, by = 'var') %>%
  mutate(order = c(6,1,2,3,4,5,7,8)) %>%
  arrange(order)
min <- min(c(coefs$value.x - coefs$se.x, coefs$value.y - coefs$se.y),
                 na.rm = T) * 1.15
max <- max(c(coefs$value.x + coefs$se.x, coefs$value.y + coefs$se.y), 
           na.rm = T) * 1.15

png('figures/sparse_quantile_regression_comparison.png', width = 8, height = 4, 
    units = 'in', res = 300)
  par(mfrow = c(1,2),
      mar = c(4,0.9,1,0),
      oma = c(0,10,2,1))
  plot(coefs$value.y, rev(seq(1:nrow(coefs))), pch = 19, xlim = c(min,max),
       xlab = expression(paste( 'Coefficient estimate (', beta, ')')), 
       bty = 'n', yaxt = 'n', ylab = '', col = coefs$col)
  segments(x0 = coefs$value.y - coefs$se.y, y0 = rev(seq(nrow(coefs):1)),
           x1 = coefs$value.y + coefs$se.y, y1 = rev(seq(nrow(coefs):1)), 
           col = coefs$col)
  abline(v = 0)
  mtext('NEP Model', line = 1.5)
  mtext(paste0('R1 = ', round(R1_NEP, 2)), line = 0.5, cex = 0.8)
  axis(2, at = seq(nrow(coefs):1), 
       labels = rev(coefs$lab), las = 2)
  plot(coefs$value.x, rev(seq(1:nrow(coefs))), pch = 19, xlim = c(min,max),
       xlab = expression(paste( 'Coefficient estimate (', beta, ')')), 
       bty = 'n', yaxt = 'n', ylab = '', col = coefs$col)
  segments(x0 = coefs$value.x - coefs$se.x, y0 = rev(seq(nrow(coefs):1)),
           x1 = coefs$value.x + coefs$se.x, y1 = rev(seq(nrow(coefs):1)), 
           col = coefs$col)
  abline(v = 0)
  mtext('P:R Model', line = 1.5)
  mtext(paste0('R1 = ', round(R1_PR, 2)), line = 0.5, cex = 0.8)
dev.off()

s <- 0.1
png('figures/sparse_quantile_regression_comparison.png', width = 8, height = 4, 
    units = 'in', res = 300)
  par(mfrow = c(1,1),mar = c(4,0,1,0),
      oma = c(0,11,2,1))
  plot(coefs$value.y, rev(seq(1:nrow(coefs))) + s, pch = 19, xlim = c(min,max),
       xlab = expression(paste( 'Coefficient estimate (', beta, ')')), 
       bty = 'n', yaxt = 'n', ylab = '', col = coefs$col)
  segments(x0 = coefs$value.y - coefs$se.y, y0 = rev(seq(nrow(coefs):1))+s,
           x1 = coefs$value.y + coefs$se.y, y1 = rev(seq(nrow(coefs):1))+s, 
           col = coefs$col)
  abline(v = 0)
  axis(2, at = seq(nrow(coefs):1), 
       labels = rev(coefs$lab), las = 2)
  points(coefs$value.x, rev(seq(1:nrow(coefs)))-s, pch = 1, xlim = c(min,max),
       xlab = expression(paste( 'Coefficient estimate (', beta, ')')), 
       bty = 'n', yaxt = 'n', ylab = '', col = coefs$col)
  segments(x0 = coefs$value.x - coefs$se.x, y0 = rev(seq(nrow(coefs):1))-s,
           x1 = coefs$value.x + coefs$se.x, y1 = rev(seq(nrow(coefs):1))-s, 
           col = coefs$col)
  abline(v = 0)
  mtext('PR Coefficients', line = 1)
dev.off()


# make a table for the SI with all of the covariates included in the model and their estimates.
coef_tab <- coefs %>% 
  rename(NEP_est = value.y, NEP_se = se.y,
         NEP_pval = p_val.y, NEP_tval = t_val.y,
         PR_est = value.x, PR_se = se.x,
         PR_pval = p_val.x, PR_tval = t_val.x,
         variable = var) %>%
  select(-col, -order)
qrtab_NEP <- data.frame(qq_NEP$coef) %>%
  mutate(variable = row.names(qq_NEP$coef)) %>% tibble() 

tmp <- dat %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(across(where(is.numeric), median, na.rm = T)) %>%
  select(-reach_proportion) %>%
  rename(NEP = ann_NEP_C)%>%
  filter(!is.na(drainage_density),
         !is.na(PrecipWs),
         !is.na(Stream_PAR_sum),
         !is.na(width_to_area)) %>%
  select(-ann_GPP_C, -ann_ER_C, -NHD_TIDAL, -site_name, -year,
         -lat, -lon, -NEP, -PR)  %>%
  # summarize(across(everything(), min)) %>%
  summarize(across(everything(), .fns = list(~min(.), ~max(.)))) %>%
  pivot_longer(cols = everything(), values_to = 'val', 
               names_to = c('variable', 'stat'), 
               names_pattern = '^(.+)?_([12])$') %>%
  pivot_wider(id_cols = variable, names_from = stat, values_from = val) %>%
  rename(min = `1`, max = `2`)

qrtab <- qrtab %>%
  select(variable) %>%
  left_join(tmp, by = 'variable')
tab_ord <- read_csv('data_working/sparse_quantile_regression_results.csv')
qrtab <- left_join(tab_ord, qrtab)
write_csv(qrtab, 'data_working/sparse_quantile_regression_results.csv')


