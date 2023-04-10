# Autotrophy in Rivers

This Project is to study rivers that are exceptions to the long-standing dogma in stream ecology that streams tend to be heterotrophic, respiring organic carbon that they receive from upstream ecosystems.  that includes metabolism estimates for 356 streams throughout the United States. 


**Contents**
  
1. [Data Sets](#data-sets-description)  
    - [Data 356 Rivers](#data-356-rivers)  
    - [Working Datasets](#working-datasets)  
2. [Workflow](#workflow)  
    - [Data Collection and Munging](#data-collection-and-munging)  
    - [Description of Dataset](#description-of-dataset)  
    - [Researcher Constrained Model](#researcher-constrained-model)  
    - [Sparse Model](#sparse-model)  

<!-- Data Sets description -->
## Data Sets Description

1.  [USGS Powell Center Synthesis](https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389) - Oxygen data from 356 USGS sites with metabolism estimates based on StreamMetabolizer [(Appling et al 2018)](https://www.nature.com/articles/sdata2018292). 
2. [StreamPulse data release](https://figshare.com/articles/software/Code_and_RDS_data_for_Bernhardt_et_al_2022_PNAS_/19074140?backTo=/collections/Data_and_code_for_Bernhardt_et_al_2022_PNAS_/5812160) and associated publication [Bernhardt et al 2022](https://www.pnas.org/doi/abs/10.1073/pnas.2121976119). We follow the workflow in this paper for filtering the Powell Center Metabolism data and for gap-filling site years for annual synthesis. Data for watershed terrestrial NPP are from this data release.
3. National Hydrography Dataset [(NHDplusHR)](https://www.usgs.gov/national-hydrography/nhdplus-high-resolution) - Each site is paired to a location on a river network in this database in order to extract covariates and delineate watersheds. 
4. [StreamCat](https://www.epa.gov/national-aquatic-resource-surveys/streamcat-dataset) - Land use and geologic characteristics summarized over the contributing watershed areas. 
5. [Global River Dissolved Oxygen Dataset](https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5) - includes covariates for landuse and watershed characteristics from Hydroatlas and NLCD.


<!-- Data 356 Rivers -->
### Data 356 Rivers

**1. high_quality_daily_metabolism_with_SP_covariates.rds**  -  Daily metabolism estimates only for sites that meet the quality filtering requirements: 1) days with poor fits for GPP, ER, or K600 are removed, 2) site years with a high correlation between K600 and ER are removed. 3) site years with less than 60% coverage of high quality days were removed. Includes raw metabolism values, filtered values, and gap-filled values according to the Bernhardt 2022 workflow.  
**2. watershed_summary_data.csv**  -  NHD and StreamCat data summarized for all Powell Center Synthesis and StreamPulse sites.  
**3. streamcat_variablelist_quickreference.csv**  -  Reference list for streamcat variable names (found in watershed_summary_data.csv).  
**4. site_data.tsv**  -  Metadata for Powell Center and StreamPulse sites.  

<!-- Working data -->
### Working datasets
**1. autotrophic_siteyears_daily.csv**  -  Daily data from Powell Center and StreamPulse sites from site-years that are autotrophic at the annual timescale.  
**2. autotrophic_siteyears_annual.csv**  -  Annual summaries of the above dataset.  
**3. across_sites_model_data.csv**  -  summary data from all site years with a minimum of 60% annual coverage with high quality days (see above). Summary metabolism values are calculated from gap-filled data using the workflow from [Bernhardt et al 2022](https://www.pnas.org/doi/abs/10.1073/pnas.2121976119). Watershed data from NHD and streamcat are also included.  


<!-- Workflow -->
## Workflow

<!-- data collection and munging -->
### Data Collection and Munging

**1. summarize_annual_data.R **  
    - Summarize annual data from StreamPulse data release for site years with at least 60% coverage of high quality estimates.   
    - Summaries are based on gap-filled data calculated using the Bernhardt 2022 workflow.  
    
**2. filter_powell_estimates_merge_streampulse.R**  
    - Remove metabolism estimates where there was poor model convergence ($\hat{r} > 1.05$)  
    - Remove estimates where: $GPP < 0$ or $ER > 0$  
    - Total estimates removed: $ER = 13\%, GPP = 14\%$  
    - Pairs Powell center estimates with annual summary covariates (calculated in step 1) from Bernhardt et al 2022  
    
**3. bulk_download_nhd_streamcat.R**  
    - Pairs all sites to NHD comids and VPUs  
    - NHD watersheds are discritized by reach, we calculate the distance along a reach where the site falls to correct landcover estimates from streamcat  
    - Downloads NHD covariates  
    - Downloads StreamCat covariates   
    
**4. baseflow_separation_calc_Q_metrics.R **
    - uses FlowScreen to separate baseflow based on the Eckhardt Method
    - contains an analysis of different baseflow index thresholds for calculating interstorm intervals
    - calculated max and median interstorm interval for each siteyear
    - calculated Richards-Baker Flashiness index
    
**5. prep_watershed_data_for_model.R**  
    - Condenses streamcat and NHD variables into categories (eg. %Urban-high, %Urban-mid, and %Urban-low all combine to become %Urban). See script for specific combinations.  
    - Assigns temporally changing variables to site years based on most recent data availability (eg. 2011 - 2015 years get %Urban2011, years 2016 and on get %Urban2016).  
    - Generates the across_sites_model_data.csv file which will be used for model building.  
    
    
<!-- Description of Dataset -->
### Description of Dataset

**1. describe_plot_annual_timeseries.R**    
    - Descriptive summary of the complete site years  
    
<!-- Researcher Constrained Model -->
### Researcher Constrained Model

**1. explore_connectivity_metrics.R**  
    - test out different derived metrics of river-landscape Carbon connectivity.
    
**2. test_quantile_reg_models.R**
    - runs regressions on the 95% quantile using the quantreg package in R. 
    - models are run using different combinations of light, disturbance, and connectivity metrics
    - plot of model coefficients sorted by Model AIC
    
<!-- Sparse Model -->
### Sparse Model

