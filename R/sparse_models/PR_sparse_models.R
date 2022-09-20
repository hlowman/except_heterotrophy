## Compiled by IAO 2022-09-20 with scripts from the Modelscape sparse modeling group

require('susieR')
require('tidyverse')
require('monomvn') 
require('ggrepel')
require('ggplot2') 


select <- dplyr::select

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/across_sites_model_data.csv')
key <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')
d <- dat %>% filter(!is.na(Stream_PAR_sum), !is.na(PrecipWs)) %>%
  select(-NLCD_LUCat, -starts_with("IGBP"))
summary(d)

#Response
PR <- -d$ann_GPP_C/d$ann_ER_C

#Predictors
preds <- select(d, -site_name, -year, -ann_GPP_C, -ann_ER_C) %>% 
  mutate(across(everything(), scale))
X <- as.matrix(preds)


###load our internal functions 
# source("/project/modelscape/analyses/sparseillustrateeval_bitbucket/analyses/gemma_function.R")
source("R/sparse_models/monomvn_function.R")
# source("/project/modelscape/analyses/sparseillustrateeval_bitbucket/analyses/random_forest_function.R")
source("R/sparse_models/performance_metrics_functions.R")




# First run the models with all sites -------------------------------------


out_list_monomvn_lasso <- monomvn_function(X, PR, method='lasso', revjump = T, training_perc = 1)
out_list_monomvn_ridge <- monomvn_function(X, PR, method='ridge', revjump = T, training_perc = 1)
out_list_susie <-susie(X,PR, L = 10,max_iter = 1000)


### Extract & store betas ###
betas <- data.frame(predictor = colnames(preds))
betas$susie_beta = coef(out_list_susie)[-1]
betas$lasso_beta = out_list_monomvn_lasso$betas
betas$ridge_beta = out_list_monomvn_ridge$betas

#Do the models arrive at similar 'answers' e.g., which predictors are most important?

#lasso vs susie
betas %>%
  ggplot(aes(x=susie_beta,y=lasso_beta, label=predictor)) +
  geom_point()+
  geom_text_repel() +
  geom_abline(slope=1, intercept=0)

#ridge vs susie
betas %>%
  ggplot(aes(x=susie_beta,y=ridge_beta, label=predictor)) +
  geom_point()+
  geom_text_repel() +
  geom_abline(slope=1, intercept=0)

#lasso vs ridge
betas %>%
  ggplot(aes(x=ridge_beta,y=lasso_beta, label=predictor)) +
  geom_point()+
  geom_text_repel() +
  geom_abline(slope=1, intercept=0)


# What do the predictions actually look like? -----------------------------


### Extract & store betas ###
predictions <- data.frame(PR_observed=PR)
predictions$susie_PR = out_list_susie$fitted
predictions$lasso_PR = out_list_monomvn_lasso$y.sim
predictions$ridge_PR = out_list_monomvn_ridge$y.sim


#Do the models arrive at similar 'answers' e.g., actual vs. predicted values?

predictions %>%
  pivot_longer(-1) %>%
  ggplot(aes(x=value,y=PR_observed)) +
  geom_point()+
  geom_abline(slope=1, intercept=0)+
  facet_wrap(~name)+
  xlim(0, 2)+
  ylim(0, 2)

