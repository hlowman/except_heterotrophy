# Make a table formatted in LateX of model results for the manuscript or SI

# install.packages('kableExtra')
library(kableExtra)
library(tidyverse)
options(scipen = 0)


# read in sparse model results
qrtab <- read_csv('data_working/sparse_quantile_regression_results.csv')
qrcoefs <- read_csv('data_working/sparse_model_coefficient_estimates.csv')

# First, a table of just the sparse model results which could be in the SI:
table <- qrtab %>% select(-variable) %>%
  mutate(estimate = round(estimate, digits = 3),
         se = paste0(" $pm ",round(se, digits = 2)),
         t_val = round(t_val, digits = 2),
         p_val = round(p_val, digits = 2),
         min = as.character(signif(min, digits = 2)),
         max = paste0('- ', as.character(signif(max, digits = 2)))) %>%
  relocate(Units, .after = max)

table %>%
  kbl(caption = 'Covariates and coefficient estimates for sparse quantile regression model',
      format = 'latex',
      col.names = c("Variable", "Estimate",  "$pm se", "t value", "p value",
                    "", "range", ""), 
      align = c('l', 'r', 'l', 'r', 'r', 'r', 'l', 'l')) %>%
  kable_classic(full_width = F, html_font = 'helvetica')


# A table of combined results:
con_vars <- c('Stream_PAR_sum',
              unique(na.omit(c(con_qrtab$D_var, con_qrtab$C_var))),
              'TmeanWs', 'PAR_kurt', 'ElevWs', 'Disch_amp', 'Width')
con_vars <- data.frame(var = con_vars,
                       Variable = c('Total PAR at stream surface', 
                                    'Discharge coefficient of variation',
                                    'Richards Baker Flashiness Index',
                                    'Maximum interstorm interval', 
                                    'Drainage density of connected reaches',
                                    'MODIS annual NPP',
                                    'Stream width to watershed area ratio',
                                    'Total Precipitation',
                                    'Mean Air Temp',
                                    'PAR kurtosis',
                                    'Elevation', 
                                    'Annual Discharge Amplitude', 
                                    'Width'))

con_qrtab <- read_csv('data_working/constrained_quantile_regression_results_PR.csv')
tmp <- con_qrtab %>% filter(!is.na(L_mean)) %>%
  arrange(resp_var, AIC)

dd <- slice(tmp, 1) %>%
  mutate(Variable = con_vars$Variable[1],
         median_est = median(tmp$L_mean),
         median_se = median(tmp$L_se)) %>%
  select(Variable, best_est = L_mean, best_se = L_se,
         median_est, median_se) 

for(i in 2:4){
  
  tmp <- con_qrtab %>% filter(D_var == con_vars$var[i]) %>%
    arrange(AIC)
  
  dd2 <- slice(tmp, 1) %>%
    mutate(Variable = con_vars$Variable[i],
           median_est = median(tmp$D_mean),
           median_se = median(tmp$D_se)) %>%
    select(Variable, best_est = D_mean, best_se = D_se,
           median_est, median_se) 
  dd <- bind_rows(dd, dd2)

}

for(i in 5:8){
  
  tmp <- con_qrtab %>% filter(C_var == con_vars$var[i]) %>%
    arrange(AIC)
  
  dd2 <- slice(tmp, 1) %>%
    mutate(Variable = con_vars$Variable[i],
           median_est = median(tmp$C_mean),
           median_se = median(tmp$C_se)) %>%
    select(Variable, best_est = C_mean, best_se = C_se,
           median_est, median_se) 
  dd <- bind_rows(dd, dd2)

}
dd_pr <- dd %>%
  select(Variable, 
         value_con_pr = median_est,
         se_con_pr = median_se)

con_qrtab <- read_csv('data_working/constrained_quantile_regression_results_NEP.csv')
tmp <- con_qrtab %>% filter(!is.na(L_mean)) %>%
  arrange(resp_var, AIC)

dd <- slice(tmp, 1) %>%
  mutate(Variable = con_vars$Variable[1],
         median_est = median(tmp$L_mean),
         median_se = median(tmp$L_se)) %>%
  select(Variable, best_est = L_mean, best_se = L_se,
         median_est, median_se) 

for(i in 2:4){
  
  tmp <- con_qrtab %>% filter(D_var == con_vars$var[i]) %>%
    arrange(AIC)
  
  dd2 <- slice(tmp, 1) %>%
    mutate(Variable = con_vars$Variable[i],
           median_est = median(tmp$D_mean),
           median_se = median(tmp$D_se)) %>%
    select(Variable, best_est = D_mean, best_se = D_se,
           median_est, median_se) 
  dd <- bind_rows(dd, dd2)

}

for(i in 5:8){
  
  tmp <- con_qrtab %>% filter(C_var == con_vars$var[i]) %>%
    arrange(AIC)
  
  dd2 <- slice(tmp, 1) %>%
    mutate(Variable = con_vars$Variable[i],
           median_est = median(tmp$C_mean),
           median_se = median(tmp$C_se)) %>%
    select(Variable, best_est = C_mean, best_se = C_se,
           median_est, median_se) 
  dd <- bind_rows(dd, dd2)

}

dd_nep <- dd %>%
  select(Variable, 
         value_con_nep = median_est,
         se_con_nep = median_se)

qrcoefs <- left_join(qrcoefs, con_vars) %>%
  select(-starts_with(c('p_val', 't_val')), -order, -col, -lab, -var)

table <- full_join(dd_nep, dd_pr,  by = 'Variable') %>%
  full_join(qrcoefs, by = 'Variable') %>%
  mutate(across(starts_with('value'), ~round(., digits = 3)),
         across(starts_with('se'), ~case_when(is.na(.) ~ ' ',
                                              TRUE ~ paste(' $pm',round(., digits = 2)))))
 

table <- left_join(table, select(qrtab, Variable, min, max, Units), by = 'Variable')

write_csv(table, 'data_working/model_summary_table.csv')

table %>%
  kbl(caption = 'Covariates and coefficient estimates for constrained and sparse quantile regression model',
      format = 'latex',
       align = c('l', 'r', 'l', 'r', 'l', 'r', 'l', 'r', 'l', 'r', 'l', 'l')) %>%
  kable_classic(full_width = F, html_font = 'helvetica')

slice(table, -c(1:11)) %>%
  select(Variable) %>% c()
unique(table$Variable)

