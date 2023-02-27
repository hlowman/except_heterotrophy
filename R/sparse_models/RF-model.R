## Compiled by IAO 2022-02-17, RF modeling
## Notes: a classification workflow tutorial using tidymodels: https://www.kirenz.com/post/2021-02-17-r-classification-tidymodels/
## Look at this tutorial for regression problems (haven't read this thoroughly but I know glmnet is one of the engines interestingly enough)
# https://www.tidymodels.org/learn/models/parsnip-ranger-glmnet/


if (!require('pacman')) install.packages('pacman'); library('pacman')

pacman::p_load("tidyverse",
               "patchwork",
               "tidymodels",
               "themis",
               "ranger")


# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/across_sites_model_data.csv')
key <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')

# filter dataframe to remove NA's in the relevant variables:
apply(dat, 2, function(x) sum(is.na(x)))

dd <- dat %>%
  filter(!is.na(site_name))%>%
  group_by(site_name) %>%
  summarize(across(where(is.numeric), median, na.rm = T)) %>%
  mutate(NEP_cat = case_when(ann_NEP_C<0 ~ "heterotrophic",
                             TRUE ~ "autotrophic"))

dd_trim <- dd %>%
  dplyr::select(site_name, NEP_cat, ann_NEP_C, lat:width_to_area) %>%
  drop_na()

cor.test(dd_trim$ElevWs, dd_trim$MOD_ann_NPP)
plot(dd_trim$PAR_kurt, dd_trim$MOD_ann_NPP)

preds <- c("lat", "lon", "PAR_sum", "Stream_PAR_sum", "Wtemp_mean", "Wtemp_cv",
           "Wtemp_skew", "Wtemp_kurt", "Wtemp_amp", "Wtemp_phase", "Wtemp_ar1",
           "Disch_mean", "Disch_cv", "Disch_skew", "Disch_kurt", "Disch_amp",
           "Disch_phase", "Disch_ar1", "PAR_mean", "PAR_cv", "PAR_skew", "PAR_kurt",
           "PAR_amp", "PAR_phase", "PAR_ar1", "LAI_mean", "LAI_cv", "LAI_skew",
           "LAI_kurt", "LAI_amp", "LAI_phase", "LAI_ar1", "Width", "MOD_ann_NPP",
           "IGBP_LC_Type1Class2018", "reach_proportion", "NHD_STREAMORDE",
           "NHD_SLOPE", "NHD_TIDAL", "ws_area_km2", "ElevWs", "PrecipWs",
           "TminWs", "TmaxWs", "TmeanWs", "RunoffWs", "HUDen2010Ws", "PopDen2010Ws",
           "Inorg_N_WetDep_kgNhayr_2008", "RdDensWs", "PctCarbResidWs",
           "PctNonCarbResidWs", "PctAlkIntruVolWs", "PctSilicicWs", "PctExtruVolWs",
           "PctColluvSedWs", "PctGlacTilClayWs", "PctGlacTilLoamWs", "PctGlacTilCrsWs",
           "PctGlacLakeCrsWs", "PctGlacLakeFineWs", "PctHydricWs", "PctEolCrsWs",
           "PctEolFineWs", "PctSalLakeWs", "PctAlluvCoastWs", "PctCoastCrsWs", "PctWaterWs",
           "WtDepWs", "OmWs", "PermWs", "RckDepWs", "BFIWs", "HydrlCondWs",
           "NWs", "Al2O3Ws", "CaOWs", "Fe2O3Ws", "K2OWs", "MgOWs", "Na2OWs",
           "P2O5Ws", "SWs", "SiO2Ws", "precip_runoff_ratio", "NLCD_PctUrban",
           "NLCD_PctAgriculture", "NLCD_PctForest", "NLCD_PctBarren", "NLCD_PctWater",
           "NLCD_PctWetland", "NLCD_PctGrassland", "Inorg_N_fert_kgNhayr",
           "Org_N_fert_kgNhayr", "Dam_densityperkm2", "Dam_total_vol_m3km2",
           "Dam_normal_vol_m3km2", "Waste_point_srcs_perkm2", "Pct_impcov",
           "connected_flow_length", "total_flow_length", "drainage_density",
           "drainage_density_connected", "med_interstorm", "max_interstorm",
           "width_to_area")

#What's the ratio of auto to heterotrophic sites?
dd_trim %>% count(NEP_cat)
#Yikes... 

split_d <- initial_split(dd_trim, strata = NEP_cat, prop=0.75)
train_d <- training(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
test_d<- testing(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
val_d <- validation_split(train_d, 
                          strata = NEP_cat, 
                          prop = 0.8)

train_d <- train_d %>% select(ann_NEP_C, site_name, all_of(preds))
test_d <- test_d %>% select(ann_NEP_C, site_name, all_of(preds))






# Try random forest REGRESSION --------------------------------------------


cores <- parallel::detectCores()
cores

rf_mod <-
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_engine("ranger", num.threads = cores) %>%
  set_mode("regression") #classification if categorical
# 
rf_recipe <-
  recipe(ann_NEP_C ~ ., data = train_d) %>%#Unlike MLR, doesn't require dummy or normalized predictor variables
  update_role(site_name, new_role = "ID") #Specify that this is an identifier
# 
rf_workflow <-
  workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(rf_recipe)
# 
# #TRAIN AND TUNE THE MODEL
rf_mod #we  have 2 hyperparameters for tuning
rf_mod %>%
  parameters()
# 
# #Use a space-filling design to tune, with 25 candidate models
# set.seed(345)
rf_res <-
  rf_workflow %>%
  tune_grid(val_d,
            grid = 25,
            control = control_grid(save_pred = TRUE))
#In regression, the root mean squared error and coefficient of determination are computed

#Top 5 model configuration based on rsq
rf_res %>%
  show_best(metric = c("rsq"))
autoplot(rf_res)+geom_smooth(method="lm")
##Looks like rsq/rmse improve are best with min_n =3 and mtry =106

rf_best <-
  rf_res %>%
  select_best(metric = "rsq")
rf_best #selects best model based on rsq

# #Calculate the data needed to plot the ROC curve. Possible after tuning with control_grid(save_pred=TRUE)
rf_res %>%
  collect_predictions()

## To filter the predictions for only our best random forest model, we can use the parameters argument and pass it our tibble with the best hyperparameter values from tuning, which we called rf_best:

# rf_auc <-
#   rf_res %>%
#   collect_predictions(parameters = rf_best) %>%
#   roc_curve(NEP_cat, .pred_autotrophic:.pred_heterotrophic) %>%
#   mutate(model = "Random Forest")

rf_res %>%
  collect_predictions(parameters = rf_best) %>%
  ggplot(aes(ann_NEP_C, .pred)) +
  geom_abline(lty = 2, color = "gray80", size = 1) +
  geom_point() +
  labs(
    x = "Observed ann_NEP_C",
    y = "Predicted ann_NEP_C"
  )+
  ggpubr::theme_pubr()

 
# the last model
last_rf_mod <-
  rand_forest(mtry = 103, min_n = 6, trees = 1000) %>%
  set_engine("ranger", num.threads = cores) %>%
  set_mode("regression")
# the last workflow
last_rf_workflow <-
  rf_workflow %>%
  update_model(last_rf_mod)
# # the last fit
set.seed(345)
last_rf_fit <-
  last_rf_workflow %>%
  last_fit(split_d)
last_rf_fit %>%
  collect_metrics()

# Judging model effectiveness for regression problems
# https://www.tmwr.org/performance.html


# Try random forest CLASSIFICATION --------------------------------------------
split_d <- initial_split(dd_trim, strata = NEP_cat, prop=0.75)
train_d <- training(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
test_d<- testing(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
val_d <- validation_split(train_d, 
                          strata = NEP_cat, 
                          prop = 0.8)

train_d <- train_d %>% select(NEP_cat, site_name, all_of(preds))
test_d <- test_d %>% select(NEP_cat, site_name, all_of(preds))


#Try random forest again, this time as classification
cores <- parallel::detectCores()
cores

rf_mod <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger", num.threads = cores) %>% 
  set_mode("classification")

rf_recipe <- 
  recipe(NEP_cat ~ ., data = train_d) %>%#Unlike MLR, doesn't require dummy or normalized predictor variables
  update_role(site_name, new_role = "ID") #Specify that this is an identifier

rf_workflow <- 
  workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe)

#TRAIN AND TUNE THE MODEL
rf_mod #we  have 2 hyperparameters for tuning
rf_mod %>%    
  parameters() 

#Use a space-filling design to tune, with 25 candidate models
set.seed(345)
rf_res <- 
  rf_workflow %>% 
  tune_grid(val_d,
            grid = 25,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc))

#Top 5 models
rf_res %>%
  show_best(metric = "roc_auc")
autoplot(rf_res)+geom_smooth(method="lm")
## Uhhh?


#Looks like roc_auc decreaes as minimum node size increases, similar pattern for # randomly selected predictors
rf_best <-
  rf_res %>% 
  select_best(metric = "roc_auc")
rf_best #selects best model based on roc_auc

#Calculate the data needed to plot the ROC curve. Possible after tuning with control_grid(save_pred=TRUE)
rf_res %>% 
  collect_predictions()

#To filter the predictions for only our best random forest model, we can use the parameters argument and pass it our tibble with the best hyperparameter values from tuning, which we called rf_best:
rf_auc <- 
  rf_res %>% 
  collect_predictions(parameters = rf_best) %>% 
  roc_curve(NEP_cat, .pred_autotrophic:.pred_heterotrophic) %>% 
  mutate(model = "Random Forest") 
#^^ Doesn't work and I'm not sure why yet

#So start by rebuilding parsnip model object from scratch and add a new argument (impurity) to get VI scores
# the last model
last_rf_mod <- 
  rand_forest(mtry = 47, min_n = 7, trees = 1000) %>% 
  set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
  set_mode("classification")

# the last workflow
last_rf_workflow <- 
  rf_workflow %>% 
  update_model(last_rf_mod)

# the last fit
set.seed(345)
last_rf_fit <- 
  last_rf_workflow %>% 
  last_fit(split_d)

last_rf_fit %>% 
  collect_metrics()

#Get VI scores
library(vip)
vip_plot<-last_rf_fit %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip(num_features = 10)
vip_plot


#Plot ROC curve
last_rf_fit %>% 
  collect_predictions() %>% 
  roc_curve(NEP_cat, .pred_autotrophic:.pred_heterotrophic) %>% 
  autoplot()
# I guess roc_curve just doesn't work for binary classifications?

#This might give us some options for evaluating performance:
#https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html