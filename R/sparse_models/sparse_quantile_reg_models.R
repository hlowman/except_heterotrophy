# This code is adapted from the constrained quantile regression model 
#   It refits the quantile regressions using a built in lasso method and including the whole dataset.

# A Carter


library(tidyverse)
library(quantreg)
library(corrplot)

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
    filter(!is.na(drainage_density),
           !is.na(PrecipWs),
           !is.na(Stream_PAR_sum),
           !is.na(width_to_area)) %>%
    select(-ann_GPP_C, -ann_ER_C, -NHD_TIDAL) %>%
    mutate(across(c(Disch_mean, Width, NHD_SLOPE, ws_area_km2, ElevWs,
                    starts_with(c('Pct', 'NLCD_Pct', 'Dam')), OmWs, HydrlCondWs,
                    NWs, CaOWs, P2O5Ws, SWs, precip_runoff_ratio,
                    Inorg_N_fert_kgNhayr, Org_N_fert_kgNhayr, PopDen2010Ws,
                    Waste_point_srcs_perkm2, connected_flow_length, 
                    total_flow_length), ~log(.+1)))%>%
    mutate(across(c(-site_name, -year), ~scale(.)[,1]))

hist2 <- function(x,y){
  hist(x, main = y)
}
par(mar = c(2,2,1,1),
    mfrow = c(5,5))
tmp <- select(dd, where(is.numeric)) 
mapply(hist2, x = tmp, y = colnames(tmp))
# test sparse quantile regression:####
summary(lm(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv, dd))
fmlqmod <- as.formula(
  paste('NEP ~', paste(colnames(dd[,7:ncol(dd)]), collapse = '+'), sep = ' '))
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv,
qmod <- quantreg::rq(fmlqmod,
                     tau =  0.95, 
                     data = dd)

lqmod <- quantreg::rqss(fmlqmod,
                        tau =  0.95, method = 'lasso', lambda = 10,
                        data = dd)
summary(lqmod, se = 'boot')
qq <- summary(lqmod, se = 'boot')
colnames(qq$coef)<- c('value', 'se', 't_val', 'p_val')
coefs <- data.frame(qq$coef) %>%
  arrange((abs(value))) %>%
  filter(p_val < 0.99) %>%
  slice(-6) 

coefs$col <- c('red', 'black', 'black', 'grey', 'grey')

png('figures/sparse_quantile_regression_results.png', width = 5, height = 4, units = 'in', res = 300)
par(mar = c(4,7,1,1))
plot(coefs$value, seq(nrow(coefs):1), pch = 19, xlim = c(-0.2, 0.25),
     xlab = expression(paste( 'Coefficient estimate (', beta, ')')), 
     bty = 'n', yaxt = 'n', ylab = '', col = coefs$col)
segments(x0 = coefs$value - coefs$se, y0 = seq(nrow(coefs):1),
         x1 = coefs$value + coefs$se, y1 = seq(nrow(coefs):1), 
         col = coefs$col)
axis(2, at = seq(nrow(coefs):1), 
     labels = rev(c("Elevation", "PAR kurtosis", 
                "Max interstorm", "Width", 
                "Terrestrial NPP")), las = 2)
abline(v = 0)
dev.off()


# make a table for the SI with all of the covariates included in the model and their estimates.
qrtab <- data.frame(qq$coef) %>%
  mutate(variable = row.names(qq$coef)) %>% tibble() 

tmp <-dat %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(across(where(is.numeric), median, na.rm = T)) %>%
  select(-reach_proportion) %>%
  rename(NEP = ann_NEP_C)%>%
  filter(!is.na(drainage_density),
         !is.na(PrecipWs),
         !is.na(Stream_PAR_sum),
         !is.na(width_to_area)) %>%
  select(-ann_GPP_C, -ann_ER_C, -NHD_TIDAL, -site_name, -year, -lat, -lon, -NEP, -PR)  %>%
  summary()

  summarize_all(.funs = list(min, max)) %>%
  pivot_longer(cols = everything(), values_to = 'val', 
               names_to = c('Variable', 'stat'), 
               names_pattern = '(^[A-Za-z0-9*]+)_(fn[12]$)')
  
write_csv(qrtab, 'data_working/sparse_quantile_regression_results_raw.csv')

qrtab <- read_csv('data_working/sparse_quantile_regression_results.csv')
install.packages('kableExtra')
library(kableExtra)

table <- qrtab %>% select(-variable) %>%
  mutate(estimate = round(estimate, digits = 3),
         se = paste0(" $pm ",round(se, digits = 2)),
         t_val = round(t_val, digits = 2),
         p_val = round(p_val, digits = 2)) 

table %>%
  kbl(caption = 'Covariates and coefficient estimates for sparse quantile regression model',
      format = 'latex',
      col.names = c("Variable", "Estimate",  "$pm se", "t value", "p value"), 
      align = c('l', 'r', 'l', 'r', 'r')) %>%
  kable_classic(full_width = F, html_font = 'helvetica')
