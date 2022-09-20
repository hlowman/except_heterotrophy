require(Metrics)


# True positive rate
tpr <- function(actual, estimated) {
  p <- as.numeric(actual != 0)
  pp <- as.numeric(estimated != 0)
  tp <- sum(p * pp)
  return(tp/sum(p))
}

# False positive rate
fpr <- function(actual, estimated) {
  n <- as.numeric(actual == 0)
  pp <- as.numeric(estimated != 0)
  fp <- sum(n * pp)
  return(fp/sum(n))
}

# Precision
precision <- function(actual, estimated) {
  p <- as.numeric(actual != 0)
  pp <- as.numeric(estimated != 0)
  tp <- sum(p * pp)
  return(tp/sum(pp))
}

# F1 score
f1 <- function(actual, estimated) {
  precision <- precision(actual, estimated)
  recall <- tpr(actual, estimated)
  return(2 / (1/precision + 1/recall))
}


# Calculate true and false positive rates
error_rates <- function(actual, estimated) {
  
  # actual: vector of true parameter values
  # estimated: vector or matrix (N parameters x N samples) of parameter estimates
  
  if(is.vector(estimated)) {
    tpr <- tpr(actual, estimated)
    fpr <- fpr(actual, estimated)
  } else {
    tpr <- apply(estimated, 2, function(i) tpr(actual, i))
    fpr <- apply(estimated, 2, function(i) fpr(actual, i))
  }
  return(list(true_positive_rate = tpr, false_positive_rate = fpr))
}


# Calculate precision, recall, and F1 score
F1 <- function(actual, estimated) {
  
  # actual: vector of true parameter values
  # estimated: vector or matrix (N parameters x N samples) of parameter estimates
  
  if(is.vector(estimated)) {
    precision <- precision(actual, estimated)
    recall <- tpr(actual, estimated)
    f1 <- f1(actual, estimated)
  } else {
    precision <- apply(estimated, 2, function(i) precision(actual, i))
    recall <- apply(estimated, 2, function(i) tpr(actual, i))
    f1 <- apply(estimated, 2, function(i) f1(actual, i))
  }
  return(list(F1 = f1, precision = precision, recall = recall))
}


# Calculate RMSE
rmse2 <- function(actual, estimated) {
  
  # actual: vector of true values
  # estimated: vector or matrix (N parameters x N samples) of predictions/estimates
  
  if(is.vector(estimated)) {
    return(Metrics::rmse(actual, estimated))
  } else {
    return(apply(estimated, 2, function(i) Metrics::rmse(actual, i)))
  }
}


# Calculate coverage of parameter estimates 
# (proportion of parameters whose uncertainty intervals include the actual value)
coverage <- function(actual, estimated, interval = 0.95) {
  
  # actual: vector of true parameter values
  # estimated: matrix (N parameters x N samples) of parameter estimates
  
  ci <- apply(estimated, 1, function(i) quantile(i, c((1-interval)/2, 1-((1-interval)/2))))
  in_interval <- sapply(1:length(actual), function(i) actual[i] >= ci[1,i] & actual[i] <= ci[2,i])
  return(sum(in_interval)/length(actual))
}


# Calculate where (which percentile) actual parameter values fall within estimated uncertainty intervals
percentile <- function(actual, estimated) {
  
  # actual: vector of true parameter values
  # estimated: matrix (N parameters x N samples) of parameter estimates
  
  return(sapply(1:length(actual), function(i) mean(estimated[i,] < actual[i])))
}


# Calculate R squared
rsq <- function(actual, estimated) {
  
  # actual: vector of actual values
  # estimated: vector or matrix of predictions/estimates
  
  if(is.vector(estimated)) {
    return(cor(actual, estimated)^2)
  } else {
    return(apply(estimated, 2, function(i) cor(actual, i)^2))
  }
}


# Calculate performance metrics (true and false positive rates, RMSE, coverage) for parameter estimates
parameter_metrics <- function(actual, estimated, ...) {
  
  # actual: vector of true parameter values
  # estimated: vector or matrix (N parameters x N samples) of parameter estimates
  
  er <- error_rates(actual, estimated)
  rmserr <- rmse2(actual, estimated)
  
  tpr <- mean(er$true_positive_rate)
  tpr_sd <- if(is.matrix(estimated)) sd(er$true_positive_rate) else NA
  fpr <- mean(er$false_positive_rate)
  fpr_sd <- if(is.matrix(estimated)) sd(er$false_positive_rate) else NA
  rmse <- mean(rmserr)
  rmse_sd <- if(is.matrix(estimated)) sd(rmserr) else NA
  coverage <- if(is.matrix(estimated)) coverage(actual, estimated, ...) else NA
  
  return(c(tpr = tpr, 
           tpr_sd = tpr_sd, 
           fpr = fpr, 
           fpr_sd = fpr_sd,
           rmse = rmse,
           rmse_sd = rmse_sd,
           coverage = coverage))
}


# Calculate performance metrics (RMSE, R2) for predictions
prediction_metrics <- function(y, ypred) {
  
  # y: actual values of y
  # ypred: vector or matrix (N observations x N samples) of predictions
  
  r_sq <- rsq(y, ypred)
  rmserr <- rmse2(y, ypred)
  
  rmse <- mean(rmserr)
  rmse_sd <- if(is.matrix(ypred)) sd(rmserr) else NA
  r2 <- mean(r_sq)
  r2_sd <- if(is.matrix(ypred)) sd(r_sq) else NA
  
  
  return(c(rmse = rmse,
           rmse_sd = rmse_sd,
           r2 = r2,
           r2_sd = r2_sd))
}