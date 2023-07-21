# Make a table formatted in LateX of model results for the manuscript or SI

# install.packages('kableExtra')
library(kableExtra)
options(scipen = 0)


# read in sparse model results
qrtab <- read_csv('data_working/sparse_quantile_regression_results.csv')
con_qrtab <- read_csv('data_working/constrained_quantile_regression_results.csv')


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
              unique(na.omit(c(con_qrtab$D_var, con_qrtab$C_var))))
con_vars <- data.frame(var = con_vars,
                       Variable = c('Total PAR at stream surface', 
                                    'Discharge coefficient of variation',
                                    'Richards Baker Flashiness Index',
                                    'Maximum interstorm interval', 
                                    'Drainage density of connected reaches',
                                    'MODIS annual NPP',
                                    'Stream width to watershed area ratio',
                                    'Total Precipitation'))

tmp <- con_qrtab %>% filter(!is.na(L_mean)) %>%
  arrange(AIC)

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

table <- full_join(dd, arrange(qrtab, desc(abs(estimate))), by = 'Variable') %>%
  mutate(best_est = round(best_est, digits = 3),
         best_se = paste0(" $pm ", round(se, digits = 2)),
         median_est = round(median_est, digits = 3),
         median_se = paste0(" $pm ", round(median_se, digits = 2)),
         estimate = round(estimate, digits = 3),
         se = paste0(" $pm ",round(se, digits = 2)),
         p_val = round(p_val, digits = 2),
         min = as.character(signif(min, digits = 2)),
         max = paste0('- ', as.character(signif(max, digits = 2)))) %>%
  filter(Variable != 'Intercept') %>%
  rename(sparse_est = estimate, sparse_se = se) %>%
  select(-t_val, -variable) %>%
  relocate(best_est, best_se, median_est, median_se, .after = Variable) %>%
  relocate(Units, .after = max) 

table %>%
  slice(1:11) %>%
  kbl(caption = 'Covariates and coefficient estimates for constrained and sparse quantile regression model',
      format = 'latex',
      col.names = c("Variable", "best Estimate",  "$pm se",
                    "median Estimate",  "$pm se",
                    "sparse Estimate",  "$pm se", "p value",
                    "", "range", ""), 
      align = c('l', 'r', 'l', 'r', 'l', 'r', 'l', 'r', 'r', 'l', 'l')) %>%
  kable_classic(full_width = F, html_font = 'helvetica')

slice(table, -c(1:11)) %>%
  select(Variable) %>% c()
unique(table$Variable)
