## Correlation Analysis Script
## October 20, 2021
## Betsy Summers, Heili Lowman

# Load packages
library(tidyverse)
library(here)

# Load data
dat <- read_csv("selected_autotrophic_rivers_daily.csv")

# Select Canyon Creek data
cc11 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2011)

cc12 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012)

# Make a quick plot of GPP and ER
figure1 <- ggplot(cc11, aes(x = Date)) +
  geom_line(aes(y = GPP), color = "green") +
  geom_line(aes(y = ER), color = "black") +
  labs(y = "ER / GPP")

# We will be basing our analysis off of Diamond et al. 2021.
# But here are some additional resources:
# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-correlation-within-and-among-time-series.html
# http://r-statistics.co/Time-Series-Analysis-With-R.html

# Canyon Creek 2011 data
# First, we'll calculate the lagged differences, as Diamond did.
diff_c11_er <- diff(-cc11$ER,1)
diff_c11_gpp <- diff(cc11$GPP,1)

# Then, we'll calculate the cross-correlation of both the raw data,
# and of the lagged differences.
# Note, the default type = "correlation"
ccf_11 <- ccf(cc11$GPP, -cc11$ER)
ccf_11_diff <- ccf(diff_c11_gpp, diff_c11_er)

# End of script.

