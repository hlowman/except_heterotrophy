# Autotrophy Modelscape Project
# August 21, 2023
# Heili Lowman

# The following script will be used to calculate summaries for sites
# that are autotrophic on an annual scale.

# Load packages.
library(here)
library(tidyverse)

# Load dataset with site-years for which sites are net autotrophic.
dat <- read_csv("data_working/across_sites_model_data.csv")

# End of script.