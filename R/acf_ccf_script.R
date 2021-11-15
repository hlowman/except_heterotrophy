## Correlation Analysis Script
## October 20, 2021
## Betsy Summers, Heili Lowman

# Load packages
library(tidyverse)
library(here)
library(tseries)

# Load data
dat <- read_csv("selected_autotrophic_rivers_daily.csv")

# Select Canyon Creek data
cc11 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2011)

cc12 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012)

# Make a quick plot of GPP and ER
(figure1 <- ggplot(cc11, aes(x = Date)) +
  geom_line(aes(y = GPP), color = "green") +
  geom_line(aes(y = ER), color = "black") +
  labs(y = "ER / GPP"))

# We will be basing our analysis off of Diamond et al. 2021.
# But here are some additional resources:
# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-correlation-within-and-among-time-series.html
# http://r-statistics.co/Time-Series-Analysis-With-R.html

# Tests for stationarity
# Augmented Dickey-Fuller test
adf.test(cc11$GPP) # apparently doesn't play nice with missing data
# so, to run some of the same analyses, we will need to impute as
# Diamond did

# Kwiatkowski-Phillips-Schmidt-Shin test
kpss.test(cc11$GPP) # skips missing data
# p < 0.01 which would suggest the dataset is non-stationary
# kpss tests check for stationarity, but can do so in the presence
# of a deterministic trend: https://www.machinelearningplus.com/time-series/kpss-test-for-stationarity/

# We'll calculate the lagged differences, as Diamond did.
# Only doing it for one year for now (2011).
diff_c11_er <- diff(-cc11$ER,1)
diff_c11_gpp <- diff(cc11$GPP,1)

# And create vectors for the raw data.
cc11_er <- -cc11$ER
cc11_gpp <- cc11$GPP

# Then, we'll calculate the cross-correlation of both the raw data,
# and of the lagged differences.
# Note, the default type = "correlation"
# na.action = na.pass will skip all NAs in the dataset
ccf_11 <- ccf(cc11_gpp, cc11_er, na.action = na.pass)
ccf_11_diff <- ccf(diff_c11_gpp, diff_c11_er, na.action = na.pass)

# Transform data into dataframes
df_ccf_11 <- data.frame(year = "2011", Lag = ccf_11$lag, CCF = ccf_11$acf)
df_ccf_11_diff <- data.frame(year = "2011", Lag = ccf_11_diff$lag, CCF = ccf_11_diff$acf)

# Plot the data
ggplot(df_ccf_11) +
  geom_col(aes(x = Lag, y = CCF), width = 0.2) +
  labs(title = "Daily cross-correlation - raw data") +
  theme_bw()

ggplot(df_ccf_11_diff) +
  geom_col(aes(x = Lag, y = CCF), width = 0.2) +
  labs(title = "Daily cross-correlation - lagged differences") +
  theme_bw()

# We can also examine the autocorrelation of the data
acf(cc11$GPP)
acf(diff_c11_gpp)
# Like the adf.test function above, the acf function does not play
# nice with datasets with missing data.

# End of script.
