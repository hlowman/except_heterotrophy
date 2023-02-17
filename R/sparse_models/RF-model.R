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
  dplyr::select(site_name, ann_NEP_C, PAR_sum:width_to_area) %>%
  drop_na()

#What's the ratio of auto to heterotrophic sites?
dd_trim %>% count(NEP_cat)
#Yikes... 

split_d <- initial_split(dd_trim, strata = ann_NEP_C, prop=0.50)
train_d <- training(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
test_d<- testing(split_d)%>%  mutate_if(is.numeric, round, digits=2) 
## I doubled checked and at least 25% of each Trend group is set aside for validation
val_d <- validation_split(train_d, 
                          strata = ann_NEP_C, 
                          prop = 0.8)


#Try random forest again
cores <- parallel::detectCores()
cores
rf_mod <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger", num.threads = cores) %>% 
  set_mode("regression") #classification if categorical

rf_recipe <- 
  recipe(ann_NEP_C ~ ., data = train_d) %>%#Unlike MLR, doesn't require dummy or normalized predictor variables
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
## ^^ Note this doesn't work (yet)

#So start by rebuilding parsnip model object from scratch and add a new argument (impurity) to get VI scores
# the last model
last_rf_mod <- 
  rand_forest(mtry = 4, min_n = 3, trees = 1000) %>% 
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



# Plots -------------------------------------------------------------------

#Get VI scores

# vip_plot<-last_rf_fit %>% 
#   pluck(".workflow", 1) %>%   
#   extract_fit_parsnip() %>% 
#   vip::vip(num_features = 10)
# 
# vip_plot
# 
# 
# 
# #Plot ROC curve
# last_rf_fit %>% 
#   collect_predictions() %>% 
#   roc_curve(NEP_cat, .pred_autotrophic) %>% 
#   autoplot()
# fit_rf<-as.data.frame(last_rf_fit %>%
#                         pluck(".predictions"))
# require(multiROC)
# true_label <- data.frame(dummies::dummy(test_d$Trend_new))
# colnames(true_label) <- c("Negative","NoTrend","Positive")
# colnames(true_label) <- paste(colnames(true_label), "_true", sep="")
# rf_pred <- fit_rf[,1:3]
# colnames(rf_pred) <- c("Negative","NoTrend","Positive")
# colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF", sep="")
# final_df <- cbind(true_label, rf_pred)
# roc_res <- multi_roc(final_df, force_diag=T)
# plot_roc_df <- plot_roc_data(roc_res)
# aucs <- plot_roc_df %>%
#   select(AUC, Method, Group) %>%
#   filter(!Group %in% c('Micro','Macro'))%>%
#   distinct()
# ROC<-plot_roc_df %>%
#   filter(!Group %in% c('Micro','Macro'))%>%
#   ggplot(., aes(x=1-Specificity,y=Sensitivity,color = Group)) +
#   geom_step() +
#   geom_text(data=aucs[aucs$Group=='NoTrend',], aes(x=0.2,y=1, label=paste0('AUC = ',round(AUC,2))), show.legend = FALSE, size=3) +
#   geom_text(data=aucs[aucs$Group=='Negative',], aes(x=0.2,y=.95, label=paste0('AUC = ',round(AUC,2))), show.legend = FALSE, size=3) +
#   geom_text(data=aucs[aucs$Group=='Positive',], aes(x=0.2,y=.9, label=paste0('AUC = ',round(AUC,2))), show.legend = FALSE, size=3) +
#   scale_color_manual(values = trendColors_a) +
#   geom_abline(slope=1,intercept=0, linetype="dashed")+
#   theme_few()
# ROC
# # Confusion matrix
# confMatRF<-confusionMatrix(fit_rf$'.pred_class', test_d$Trend_new)
# confMatRF
