# Sparse Model sandbox for annual Metabolism as a function of watershed characteristics
# A Carter
# March 2022

library(tidyverse)
library(susieR)
library(glmnet)

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/across_sites_model_data.csv')
d <- dat %>% filter(!is.na(Stream_PAR_sum), !is.na(PrecipWs)) %>%
  select(-NLCD_LUCat, -starts_with("IGBP"))
summary(d)

GPP <- d$ann_GPP_C
ER <- d$ann_ER_C
NEP <- scale(GPP + ER)
GPP <- scale(GPP)
preds <- select(d, -site_name, -year, -ann_GPP_C, -ann_ER_C) %>% 
  mutate(across(everything(), scale))
X <- as.matrix(preds)
# fit test sparse models ####
# linear regression

linreg <- lm(NEP~X)
summary(linreg)
betas <- data.frame(predictor = colnames(preds),
           lin_reg_beta = coef(linreg)[-1])

# susie
res <- susie(X, NEP, L=20)
betas$susie_L20_beta = abs(coef(res)[-1])

par(mfrow = c(3,2), mar = c(2,4,1.5,1), oma = c(3,1,1,1))
plot(betas$susie_L20_beta, betas$lin_reg_beta, ylab = 'linear regression beta',
     xlab = '', xlim = c(0,0.6), pch = 19, main = 'NEP coefficients')
abline(0,1)
identify(betas$susie_L20_beta, betas$lin_reg_beta, labels = betas$predictor,
         atpen = T)

# GPP
linreg <- lm(GPP~X)
betas$lin_reg_betag = coef(linreg)[-1]

# susie
res <- susie(X, GPP, L=20)
betas$susie_L20_betag = abs(coef(res)[-1])

plot(betas$susie_L20_betag, betas$lin_reg_betag, ylab = 'linear regression beta',
     xlab = '', xlim = c(0,0.7), pch = 19, main = 'GPP coefficients')
abline(0,1)
identify(betas$susie_L20_betag, betas$lin_reg_betag, labels = betas$predictor,
         atpen = T)
# lasso
lasso <- glmnet(X, NEP, alpha = 1, nlambda = 5)
betas$lasso_a1nl5_beta1 = abs(coef(lasso)[-1,1])
betas$lasso_a1nl5_beta2 = abs(coef(lasso)[-1,2])
betas$lasso_a1nl5_beta3 = abs(coef(lasso)[-1,3])

plot(betas$susie_L20_beta, betas$lasso_a1nl5_beta3, pch = 19, xlim = c(0,.6),
     ylab = 'Lasso beta')
abline(0,1)
identify(betas$susie_L20_beta, betas$lasso_a1nl5_beta3, 
         labels = betas$predictor, atpen = T)

points(betas$susie_L20_beta, betas$lasso_a1nl5_beta2, pch = 19, col = 2)
points(betas$susie_L20_beta, betas$lasso_a1nl5_beta1, pch = 19, col = 3)

# lasso GPP
lasso <- glmnet(X, GPP, alpha = 1, nlambda = 5)
betas$lasso_a1nl5_beta1g = abs(coef(lasso)[-1,1])
betas$lasso_a1nl5_beta2g = abs(coef(lasso)[-1,2])
betas$lasso_a1nl5_beta3g = abs(coef(lasso)[-1,3])

plot(betas$susie_L20_betag, betas$lasso_a1nl5_beta3g, pch = 19, xlim = c(0,.7),
     ylab = 'Lasso beta')
abline(0,1)
identify(betas$susie_L20_betag, betas$lasso_a1nl5_beta3g, 
         labels = betas$predictor, atpen = T)

points(betas$susie_L20_betag, betas$lasso_a1nl5_beta2g, pch = 19, col = 2)
points(betas$susie_L20_betag, betas$lasso_a1nl5_beta1g, pch = 19, col = 3)

# ridge
ridge <- glmnet(X, NEP, alpha = 1, nlambda = 10)
betas$ridge_a1nl10_beta1 = abs(coef(ridge)[-1,1])
betas$ridge_a1nl10_beta2 = abs(coef(ridge)[-1,2])
betas$ridge_a1nl10_beta3 = abs(coef(ridge)[-1,3])

plot(betas$susie_L20_beta, betas$ridge_a1nl10_beta3, pch = 19, xlim = c(0,.6),
     ylab = 'Ridge beta')
abline(0, 1)
identify(betas$susie_L20_beta, betas$ridge_a1nl10_beta3, atpen = T,
         labels = betas$predictor)
points(betas$susie_L20_beta, betas$ridge_a1nl10_beta2, pch = 19, col = 2)
points(betas$susie_L20_beta, betas$ridge_a1nl10_beta1, col = 3, pch = 19)
mtext('Susie beta', 1, 1, outer = T, )

# ridge GPP
ridge <- glmnet(X, GPP, alpha = 1, nlambda = 10)
betas$ridge_a1nl10_beta1g = abs(coef(ridge)[-1,1])
betas$ridge_a1nl10_beta2g = abs(coef(ridge)[-1,2])
betas$ridge_a1nl10_beta3g = abs(coef(ridge)[-1,3])

plot(betas$susie_L20_betag, betas$ridge_a1nl10_beta3g, pch = 19, xlim = c(0,.7),
     ylab = 'Ridge beta')
abline(0, 1)
identify(betas$susie_L20_betag, betas$ridge_a1nl10_beta3g, atpen = T,
         labels = betas$predictor)
points(betas$susie_L20_betag, betas$ridge_a1nl10_beta2g, pch = 19, col = 2)
points(betas$susie_L20_betag, betas$ridge_a1nl10_beta1g, col = 3, pch = 19)
