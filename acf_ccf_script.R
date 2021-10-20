## Correlation Analysis Script
## October 20
## Betsy Summers, Heili Lowman

# Load packages
library(tidyverse)
library(here)

# Load data
dat <- read_csv("selected_autotrophic_rivers_daily.csv")

# Select Canyon Creek data
cc <- dat %>%
  filter(sitecode == "nwis_10133650")
