## Correlation Analysis Script
## October 20, 2021
## Betsy Summers, Heili Lowman

# Load packages
library(tidyverse)
library(here)

# Load data
dat <- read_csv("selected_autotrophic_rivers_daily.csv")

# Select Canyon Creek data
cc <- dat %>%
  filter(sitecode == "nwis_10133650")

# Make a quick plot of GPP and ER
figure1 <- ggplot(cc, aes(x = Date)) +
  geom_line(aes(y = GPP), color = "green") +
  geom_line(aes(y = ER), color = "black") +
  labs(y = "ER / GPP")

# We will be basing our analysis off of Diamond et al. 2021.
# But here are some additional resources:
# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-correlation-within-and-among-time-series.html
# http://r-statistics.co/Time-Series-Analysis-With-R.html

