# Modified functions from BernhardtMetabolism package for getting magnificent 7 metrics
calculate_mag7 <- function (compiled_data, var, standardize) {
    if (all(is.na(compiled_data[, var])) == FALSE) {
      calc_mag7 <- mag7_fun(compiled_data, var = var, standardize = standardize)
    }
    if (all(is.na(compiled_data[, var])) == TRUE) {
      calc_mag7 <- setNames(data.frame(t(rep(NA, 7))), c("lam1", 
                                                         "tau2", "tau3", "tau4", "amplitude", "phase", "ar1"))
    }
    return(calc_mag7)
  }
  
mag7_fun <- function (timeseries, var, standardize) {
    x <- timeseries[, var]
    if(is.list(x)) x <- unlist(x)
    l_mom <- lmomco::lmom.ub(x)
    lam1 <- round(l_mom$L1, digits = 2)
    tau2 <- round(l_mom$LCV, digits = 2)
    tau3 <- round(l_mom$TAU3, digits = 2)
    tau4 <- round(l_mom$TAU4, digits = 2)
    seasonality <- get_seasonality(timeseries, var, standardize)
    ar1 <- ar_fun(timeseries, var)
    mag7 <- data.frame(lam1, tau2, tau3, tau4, seasonality, 
                       ar1)
    colnames(mag7) <- c("mean", "cv", "skew", "kurt", "amp", 
                     "phase", "ar1")
    return(mag7)
}

get_seasonality <- function (timeseries, var, standardize) {
  decimal_year <- as.numeric(timeseries[, "year"]) + (timeseries[,"DOY"]/365.25)
  ifelse(standardize == "yes", standard <- scale(timeseries[, 
                                                            var], center = TRUE, scale = TRUE), standard <- timeseries[, 
                                                                                                                       var])
  x_mat <- cbind(1, sin(2 * pi * decimal_year), cos(2 * pi * 
                                                      decimal_year))
  seasonfit <- .lm.fit(x_mat, standard)
  b1 <- as.vector(coef(seasonfit)[2])
  b2 <- as.vector(coef(seasonfit)[3])
  amplitude <- round(sqrt((b2^2) + (b1^2)), digits = 2)
  MaxDay <- function(b1, b2) {
    version1 <- 365.25 * ((pi/2) - atan(b2/b1))/(2 * pi)
    version2 <- 365.25 * ((pi/2) - pi - atan(b2/b1))/(2 * 
                                                        pi)
    MaxDay <- if (b1 > 0) 
      version1
    else 365.25 + version2
    MaxDay <- if (b1 == 0 & b2 > 0) 
      365.25
    else MaxDay
    MaxDay <- if (b1 == 0 & b2 < 0) 
      365.25/2
    else MaxDay
    MaxDay <- if (b1 == 0 & b2 == 0) 
      NA
    else MaxDay
    return(MaxDay)
  }
  phase <- MaxDay(b1, b2)
  get_seasonalityv <- cbind(amplitude, phase)
  return(get_seasonalityv)
}

 
ar_fun <- function (timeseries, var) {
   timeseries$Month <- as.numeric(format(as.Date(timeseries[,"DOY"] - 1, 
                                                 origin = paste(timeseries[, "year"],
                                                                "-01-01", sep = "")), 
                                         format = "%m"))
   monmeans <- aggregate(timeseries[, var], list(timeseries[, 
                                                            "Month"]), FUN = mean, na.rm = TRUE)
   mon_merge <- merge(timeseries, monmeans, by.x = "Month", 
                      by.y = "Group.1")
   mon_merge$deseason <- mon_merge[, var] - mon_merge[, "x"]
   ordered <- mon_merge[order(mon_merge$year, mon_merge$DOY), ]
   stand <- scale(ordered[, "deseason"], center = TRUE, scale = TRUE)
   ar_mod <- ar(stand, aic = FALSE, order.max = 1, method = "yule-walker")
   ar_cor <- round(ar_mod$ar, digits = 3)
   return(ar_cor)
 }

