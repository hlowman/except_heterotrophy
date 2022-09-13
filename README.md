# Autotrophy in Rivers

This Project is to study rivers that are exceptions to the long-standing dogma in stream ecology that streams tend to be heterotrophic, respiring organic carbon that they receive from upstream ecosystems. We are using data from a [USGS Powell Center Synthesis](https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389) that includes metabolism estimates for 356 streams throughout the United States. 


**Contents**
  
1. [Data Sets](#data-sets-description)
    - [Data 356 Rivers](#data-356-rivers)
    - [Working Datasets](#working-datasets)

<!-- Data Sets description -->
## Data Sets Description

1. Powell Center Dataset
2. StreamPulse data release
3. National Hydrography Dataset HR
4. StreamCat


<!-- Data 356 Rivers -->
### Data 356 Rivers

**1. high_quality_daily_metabolism_with_SP_covariates.rds**  -  Daily metabolism estimates only for sites that meet the quality filtering requirements: 1) days with poor fits for GPP, ER, or K600 are removed, 1) site years with a high correlation between K600 and ER are removed.  
**2. watershed_summary_data.csv**  -  NHD and StreamCat data summarized for all Powell Center Synthesis and StreamPulse sites.  
**3. streamcat_variablelist_quickreference.csv**  -  Reference list for streamcat variable names (found in watershed_summary_data.csv).  
**4. site_data.tsv**  -  Metadata for Powell Center and StreamPulse sites.  

<!-- Working data -->
### Working datasets
**1. autotrophic_siteyears_daily.csv**  -  Daily data from Powell Center and StreamPulse sites from site-years that are autotrophic at the annual timescale.  
**2. autotrophic_siteyears_annual.csv**  -  Annual summaries of the above dataset.  
**3. across_sites_model_data.csv**  -  summary data from all site years with a minimum of 60% annual coverage with high quality days (see above). Summary metabolism values are calculated from gap-filled data using the workflow from [Bernhardt et al 2022]. Watershed data from NHD and streamcat are also included.  


