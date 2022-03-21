# Sparse Model sandbox for annual Metabolism as a function of watershed characteristics
# A Carter
# March 2022

library(tidyverse)
library(susieR)
library(glmnet)

# load data
# this subset contains only the siteyears that have decent data and coverage (according to Bernhardt 2022)
dat <- read_csv('data_working/annual_summary_data.csv')


# simulate correlated predictor variables from multivariate normal distribution
n <- 1000
p <- 50
mu <- rnorm(p, 0, 1)                 # means
Sigma <- rethinking::rlkjcorr(1, p)  # covariance matrix
X <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)  # draw predictors from multivariate normal distribution  

# regression coefficients; most close to zero but some large values
beta <- rgamma(p, 0.1, 0.1) * sample(c(-1, 1), p, replace = T)
plot(ecdf(beta))

# simulate response variable
y <- X %*% beta + rnorm(n)

# fit models and plot actual vs. estimated coefficients -----------------------

# linear regression
linreg <- lm(y~X)
plot(abs(beta), abs(coef(linreg)[-1])); abline(0, 1)
cor(abs(beta), abs(coef(linreg)[-1]))

# susie
res <- susie(X, y, L=20)
plot(abs(beta), abs(coef(res)[-1])); abline(0, 1)
cor(abs(beta), abs(coef(res)[-1]))
# lasso
lasso <- glmnet(X, y, alpha = 1, nlambda = 5)
plot(abs(beta), abs(coef(lasso)[-1,1])); abline(0, 1)
cor(abs(beta), abs(coef(lasso)[-1,1]))
plot(abs(beta), abs(coef(lasso)[-1,2])); abline(0, 1)
cor(abs(beta), abs(coef(lasso)[-1,2]))
plot(abs(beta), abs(coef(lasso)[-1,3])); abline(0, 1)
cor(abs(beta), abs(coef(lasso)[-1,3]))

# ridge
ridge <- glmnet(X, y, alpha = 1, nlambda = 10)
plot(abs(beta), abs(coef(ridge)[-1,1])); abline(0, 1)
cor(abs(beta), abs(coef(ridge)[-1,1]))
plot(abs(beta), abs(coef(ridge)[-1,2])); abline(0, 1)
cor(abs(beta), abs(coef(ridge)[-1,1]))
plot(abs(beta), abs(coef(ridge)[-1,3])); abline(0, 1)
cor(abs(beta), abs(coef(ridge)[-1,1]))


