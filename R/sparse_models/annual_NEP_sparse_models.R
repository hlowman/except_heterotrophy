# Sparse Model sandbox for annual Metabolism as a function of watershed characteristics
# A Carter
# March 2022

library(tidyverse)
library(susieR)
library(glmnet)

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/across_sites_model_data.csv')
key <- read_csv('data_356rivers/streamcat_variablelist_quickreference.csv')
d <- dat %>% filter(!is.na(Stream_PAR_sum), !is.na(PrecipWs)) %>%
  select(-NLCD_LUCat, -starts_with("IGBP"))
summary(d)

d %>%  
  mutate(NEP = ann_GPP_C + ann_ER_C,
         PR = -ann_GPP_C/ann_ER_C) %>%
  ggplot(aes(ElevWs, log(ann_GPP_C), col = lat)) +
  geom_point(size = 1) +
  # scale_color_gradientn(colors = viridis::magma(5)) +
  # xlab("GPP (g O2/m2/y)") + ylab("NEP (g O2/m2/y)") +
  theme(legend.position = 'none') +
  theme_minimal() 


GPP <- log(d$ann_GPP_C)
PR <- -d$ann_GPP_C/d$ann_ER_C

preds <- select(d, -site_name, -year, -ann_GPP_C, -ann_ER_C) %>% 
  mutate(across(everything(), scale))
X <- as.matrix(preds)
# fit test sparse models ####
# linear regression
# 
# linreg <- lm(NEP~X)
# summary(linreg)
# betas <- data.frame(predictor = colnames(preds),
#            lin_reg_beta = coef(linreg)[-1])

# susie
res <- susie(X, PR, L=20)
betas <- data.frame(predictor = colnames(preds))
betas$susie_L20_beta = coef(res)[-1]

# susie GPP
res <- susie(X, GPP, L=20)
betas$susie_L20_betag = coef(res)[-1]

# lasso
lasso <- glmnet(X, PR, alpha = 1, nlambda = 5)
betas$lasso_a1nl5_beta1 = coef(lasso)[-1,1]
betas$lasso_a1nl5_beta2 = coef(lasso)[-1,2]
betas$lasso_a1nl5_beta3 = coef(lasso)[-1,3]

# lasso GPP
lasso <- glmnet(X, GPP, alpha = 1, nlambda = 5)
betas$lasso_a1nl5_beta1g = coef(lasso)[-1,1]
betas$lasso_a1nl5_beta2g = coef(lasso)[-1,2]
betas$lasso_a1nl5_beta3g = coef(lasso)[-1,3]


par(mfrow = c(2,2), mar = c(2,4,1.5,1), oma = c(3,1,1,1))

plot(betas$susie_L20_beta, betas$lasso_a1nl5_beta2, pch = 19, xlim = c(-.2,.2),
     ylim = c(-.2,.2),     ylab = 'Lasso beta',  main = 'P/R coefficients')
abline(0,1)
identify(betas$susie_L20_beta, betas$lasso_a1nl5_beta2, labels = betas$predictor,
         atpen = T)

# points(betas$susie_L20_beta, betas$lasso_a1nl5_beta2, pch = 19, col = 2)
# points(betas$susie_L20_beta, betas$lasso_a1nl5_beta1, pch = 19, col = 3)

plot(betas$susie_L20_betag, betas$lasso_a1nl5_beta2g, pch = 19, xlim = c(-.4,.65), ylim = c(-.4,.65),
     ylab = 'Lasso beta',  main = 'GPP coefficients')
abline(0,1)
identify(betas$susie_L20_betag, betas$lasso_a1nl5_beta2g, labels = betas$predictor,
         atpen = T)

# ridge
ridge <- glmnet(X, PR, alpha = 1, nlambda = 10)
betas$ridge_a1nl10_beta1 = coef(ridge)[-1,1]
betas$ridge_a1nl10_beta2 = coef(ridge)[-1,2]
betas$ridge_a1nl10_beta3 = coef(ridge)[-1,3]

plot(betas$susie_L20_beta, betas$ridge_a1nl10_beta2, pch = 19, ylim = c(-0.2,.2), xlim = c(-0.2,.2),
     ylab = 'Ridge beta')
abline(0, 1)
identify(betas$susie_L20_beta, betas$ridge_a1nl10_beta2, atpen = T,
         labels = betas$predictor)
# points(betas$susie_L20_beta, betas$ridge_a1nl10_beta2, pch = 19, col = 2)
# points(betas$susie_L20_beta, betas$ridge_a1nl10_beta1, col = 3, pch = 19)
mtext('Susie beta', 1, 1, outer = T, )

# ridge GPP
ridge <- glmnet(X, GPP, alpha = 1, nlambda = 10)
betas$ridge_a1nl10_beta1g = coef(ridge)[-1,1]
betas$ridge_a1nl10_beta2g = coef(ridge)[-1,2]
betas$ridge_a1nl10_beta3g = coef(ridge)[-1,3]

plot(betas$susie_L20_betag, betas$ridge_a1nl10_beta2g, pch = 19, xlim = c(-0.4,.65), ylim = c(-0.4,.65),
     ylab = 'Ridge beta')
abline(0, 1)
identify(betas$susie_L20_betag, betas$ridge_a1nl10_beta2g, atpen = T,
         labels = betas$predictor)
# points(betas$susie_L20_betag, betas$ridge_a1nl10_beta2g, pch = 19, col = 2)
# points(betas$susie_L20_betag, betas$ridge_a1nl10_beta1g, col = 3, pch = 19)



# Split dataframe into autotrophic and heterotrophic: ####
aut <- filter(d, ann_GPP_C + ann_ER_C >0)
het <- filter(d, ann_GPP_C + ann_ER_C <0)

GPP_aut <- log(aut$ann_GPP_C)
GPP_het <- log(het$ann_GPP_C)

preds_aut <- select(aut, -site_name, -year, -ann_GPP_C, -ann_ER_C, -PctCoastCrsWs) %>% 
  mutate(across(everything(), scale))
X_aut <- as.matrix(preds_aut)
summary(X_aut)

preds_het <- select(het, -site_name, -year, -ann_GPP_C, -ann_ER_C) %>% 
  mutate(across(everything(), scale))
X_het <- as.matrix(preds_het)
summary(X_het)

betas_aut <- data.frame(predictor = colnames(preds_aut))
betas_het <- data.frame(predictor = colnames(preds_het))
# susie
res_aut <- susie(X_aut, GPP_aut, L=20)
betas_aut$susie_L20_beta = coef(res_aut)[-1]

res_het <- susie(X_het, GPP_het, L=20)
betas_het$susie_L20_beta = coef(res_het)[-1]

# lasso
lasso_aut <- glmnet(X_aut, GPP_aut, alpha = 1, nlambda = 5)
betas_aut$lasso_a1nl5_beta = coef(lasso_aut)[-1,2]

lasso_het <- glmnet(X_het, GPP_het, alpha = 1, nlambda = 5)
betas_het$lasso_a1nl5_beta = coef(lasso_het)[-1,2]

# ridge
ridge_aut <- glmnet(X_aut, GPP_aut, alpha = 1, nlambda = 10)
betas_aut$ridge_a1nl10_beta = coef(ridge_aut)[-1,3]

ridge_het <- glmnet(X_het, GPP_het, alpha = 1, nlambda = 10)
betas_het$ridge_a1nl10_beta = coef(ridge_het)[-1,3]


par(mfrow = c(2,2), mar = c(2,4,1.5,1), oma = c(3,1,3,1))

plot(betas_aut$susie_L20_beta, betas_aut$lasso_a1nl5_beta, pch = 19, xlim = c(-.6,.6),
     ylab = 'Lasso beta',  main = 'autotrophic rivers', ylim = c(-.6,.6))
abline(0,1)
identify(betas_aut$susie_L20_beta, betas_aut$lasso_a1nl5_beta, labels = betas_aut$predictor,
         atpen = T)

plot(betas_het$susie_L20_beta, betas_het$lasso_a1nl5_beta, pch = 19, xlim = c(-.6,.6),
     ylab = 'Lasso beta',  main = 'hetotrophic rivers', ylim = c(-.6,.6))
abline(0,1)
identify(betas_het$susie_L20_beta, betas_het$lasso_a1nl5_beta, labels = betas_het$predictor,
         atpen = T)

plot(betas_aut$susie_L20_beta, betas_aut$ridge_a1nl10_beta, pch = 19, xlim = c(-.6,.6),
     ylab = 'Ridge beta', ylim = c(-.6,.6))
abline(0,1)
identify(betas_aut$susie_L20_beta, betas_aut$ridge_a1nl10_beta, labels = betas_aut$predictor,
         atpen = T)

plot(betas_het$susie_L20_beta, betas_het$ridge_a1nl10_beta, pch = 19, xlim = c(-.6,.6),
     ylab = 'Ridge beta', ylim = c(-.6,.6))
abline(0,1)
identify(betas_het$susie_L20_beta, betas_het$ridge_a1nl10_beta, labels = betas_het$predictor,
         atpen = T)

mtext('Susie beta', 1, 1, outer = T )
mtext('Coefficients for log(GPP)', 3, 1, outer = T, cex = 1.6)

