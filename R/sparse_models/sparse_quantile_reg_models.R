# This code is adapted from the constrained quantile regression model 
#   It refits the quantile regressions using a built in lasso method and including the whole dataset.

# A Carter


# update the two quantile model files to reflect that
# summarize the issue you are having with exponentially distributed data and message bella and Christa

library(tidyverse)
library(quantreg)
library(corrplot)

dat <- read_csv('data_working/across_sites_model_data.csv')

# filter dataframe to remove NA's in the relevant variables:
apply(dat, 2, function(x) sum(is.na(x)))

dd <- dat %>%
    filter(!is.na(site_name))%>%
    group_by(site_name) %>%
    summarize(across(where(is.numeric), median, na.rm = T))

# log transform highly skewed variables
dd <- dd %>%
    rename(NEP = ann_NEP_C)%>%
    filter(!is.na(drainage_density),
           !is.na(PrecipWs),
           !is.na(Stream_PAR_sum),
           !is.na(width_to_area)) %>%
    select(-ann_GPP_C, -ann_ER_C) %>%
    mutate(across(c(Disch_mean, Width, NHD_SLOPE, ws_area_km2, ElevWs, PopDen2010Ws,
                    starts_with(c('Pct', 'NLCD_Pct', 'Dam')), OmWs, HydrlCondWs, NWs, 
                                CaOWs, P2O5Ws, SWs, precip_runoff_ratio,
                                Inorg_N_fert_kgNhayr, Org_N_fert_kgNhayr,
                                Waste_point_srcs_perkm2, connected_flow_length,
                                total_flow_length), ~log(.+1)))%>%
    mutate(across(c(-site_name, -year), ~scale(.)[,1]))

hist2 <- function(x,y){
  hist(x, main = y)
}
par(mar = c(2,2,1,1),
    mfrow = c(5,5))
tmp <- select(dd, where(is.numeric)) 
mapply(hist2, x = tmp, y = colnames(tmp))
# test sparse quantile regression:####
summary(lm(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv, dd))
# qmod <- quantreg::rq(NEP ~ Stream_PAR_sum + PrecipWs + Disch_cv,
qmod <- quantreg::rq(NEP ~ lat+lon+PAR_sum+Stream_PAR_sum+Wtemp_mean+Wtemp_cv+Wtemp_skew+Wtemp_kurt+Wtemp_amp+Wtemp_phase+Wtemp_ar1+Disch_mean+Disch_cv+Disch_skew+Disch_kurt+Disch_amp+Disch_phase+Disch_ar1+PAR_mean+PAR_cv+PAR_skew+PAR_kurt+PAR_amp+PAR_phase+PAR_ar1+LAI_mean+LAI_cv+LAI_skew+LAI_kurt+LAI_amp+LAI_phase+LAI_ar1+Width+MOD_ann_NPP+IGBP_LC_Type1Class2018+reach_proportion+NHD_STREAMORDE+NHD_SLOPE+NHD_TIDAL+ws_area_km2+ElevWs+PrecipWs+TminWs+TmaxWs+TmeanWs+RunoffWs+HUDen2010Ws+PopDen2010Ws+Inorg_N_WetDep_kgNhayr_2008+RdDensWs+PctCarbResidWs+PctNonCarbResidWs+PctAlkIntruVolWs+PctSilicicWs+PctExtruVolWs+PctColluvSedWs+PctGlacTilClayWs+PctGlacTilLoamWs+PctGlacTilCrsWs+PctGlacLakeCrsWs+PctGlacLakeFineWs+PctHydricWs+PctEolCrsWs+PctEolFineWs+PctSalLakeWs+PctAlluvCoastWs+PctCoastCrsWs+PctWaterWs+WtDepWs+OmWs+PermWs+RckDepWs+BFIWs+HydrlCondWs+NWs+Al2O3Ws+CaOWs+Fe2O3Ws+K2OWs+MgOWs+Na2OWs+P2O5Ws+SWs+SiO2Ws+precip_runoff_ratio+NLCD_PctUrban+NLCD_PctAgriculture+NLCD_PctForest+NLCD_PctBarren+NLCD_PctWater+NLCD_PctWetland+NLCD_PctGrassland+Inorg_N_fert_kgNhayr+Org_N_fert_kgNhayr+Dam_densityperkm2+Dam_total_vol_m3km2+Dam_normal_vol_m3km2+Waste_point_srcs_perkm2+Pct_impcov+connected_flow_length+total_flow_length+drainage_density+drainage_density_connected+med_interstorm+max_interstorm+width_to_area,
                     tau =  0.95, 
                     data = dd)

fmlqmod <- paste('NEP ~', paste(colnames(dd[,5:ncol(dd)]), collapse  = '+'), sep = ' ')
lqmod <- quantreg::rqss(NEP ~ lat+lon+PAR_sum+Stream_PAR_sum+Wtemp_mean+Wtemp_cv+Wtemp_skew+Wtemp_kurt+Wtemp_amp+Wtemp_phase+Wtemp_ar1+Disch_mean+Disch_cv+Disch_skew+Disch_kurt+Disch_amp+Disch_phase+Disch_ar1+PAR_mean+PAR_cv+PAR_skew+PAR_kurt+PAR_amp+PAR_phase+PAR_ar1+LAI_mean+LAI_cv+LAI_skew+LAI_kurt+LAI_amp+LAI_phase+LAI_ar1+Width+MOD_ann_NPP+IGBP_LC_Type1Class2018+reach_proportion+NHD_STREAMORDE+NHD_SLOPE+NHD_TIDAL+ws_area_km2+ElevWs+PrecipWs+TminWs+TmaxWs+TmeanWs+RunoffWs+HUDen2010Ws+PopDen2010Ws+Inorg_N_WetDep_kgNhayr_2008+RdDensWs+PctCarbResidWs+PctNonCarbResidWs+PctAlkIntruVolWs+PctSilicicWs+PctExtruVolWs+PctColluvSedWs+PctGlacTilClayWs+PctGlacTilLoamWs+PctGlacTilCrsWs+PctGlacLakeCrsWs+PctGlacLakeFineWs+PctHydricWs+PctEolCrsWs+PctEolFineWs+PctSalLakeWs+PctAlluvCoastWs+PctCoastCrsWs+PctWaterWs+WtDepWs+OmWs+PermWs+RckDepWs+BFIWs+HydrlCondWs+NWs+Al2O3Ws+CaOWs+Fe2O3Ws+K2OWs+MgOWs+Na2OWs+P2O5Ws+SWs+SiO2Ws+precip_runoff_ratio+NLCD_PctUrban+NLCD_PctAgriculture+NLCD_PctForest+NLCD_PctBarren+NLCD_PctWater+NLCD_PctWetland+NLCD_PctGrassland+Inorg_N_fert_kgNhayr+Org_N_fert_kgNhayr+Dam_densityperkm2+Dam_total_vol_m3km2+Dam_normal_vol_m3km2+Waste_point_srcs_perkm2+Pct_impcov+connected_flow_length+total_flow_length+drainage_density+drainage_density_connected+med_interstorm+max_interstorm+width_to_area,
                        tau =  0.95, method = 'lasso', lambda = 2,
                        data = dd)
qq <- summary(lqmod, se = 'boot')
colnames(qq$coef)<- c('value', 'se', 't_val', 'p_val')
data.frame(qq$coef) %>%
  arrange(desc(value))
summary(lqmod, se = 'boot')

plot(qmod)
