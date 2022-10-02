# Test for relationships to light and flow across years within a site 
# This will help determine the model structure we decide to use

# A Carter

library(tidyverse)

dat <- read_csv('data_working/across_sites_model_data.csv')

glimpse(dat)

sites <- dat %>% 
    group_by(site_name)%>%
    summarize(n_years = n()) %>%
    filter(n_years > 2)

for(i in 1:nrow(sites)){
  if(i == 20) next
    d <- filter(dat, site_name == sites$site_name[i])  %>%
        mutate(Y = -ann_GPP_C/ann_ER_C) %>%
        # mutate(Y = ann_GPP_C+ann_ER_C) %>%
      select(Y, Disch_cv, PrecipWs, Stream_PAR_sum) %>%
      scale() %>%
      as.data.frame()
    
    m <- lm(Y ~ Disch_cv, data = d)
    sites$Q_CV_slope[i] = m$coefficients[2]
    f <- summary(m)$fstatistic
    try(sites$Q_CV_pval[i] <- pf(f[1],f[2],f[3],lower.tail=F))

    
    try(m <- lm(Y ~ PrecipWs, data = d))
    try(sites$precip_slope[i] <- m$coefficients[2])
    try(f <- summary(m)$fstatistic)
    try(sites$precip_pval[i] <- pf(f[1],f[2],f[3],lower.tail=F))
    
    try(m <- lm(Y ~ Stream_PAR_sum, data = d))
    try(sites$par_slope[i] <- m$coefficients[2])
    try(f <- summary(m)$fstatistic)
    try(sites$par_pval[i] <- pf(f[1],f[2],f[3],lower.tail=F))
}

plot(density(sites$Q_CV_slope, na.rm = T))
plot(density(sites$precip_slope, na.rm = T))
plot(density(sites$par_slope, na.rm = T))


ggplot(sites, aes(Q_CV_pval, Q_CV_slope))+
  geom_point()
ggplot(sites, aes(precip_pval, precip_slope))+
  geom_point()
ggplot(sites, aes(par_pval, par_slope))+
  geom_point()

 # it doesn't seem like there is sufficient evidence to include within site variation in this dataset

