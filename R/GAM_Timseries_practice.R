library(ggplot2)
library(dplyr)
library(epitools)

riverdaily<-read.csv("data_working/autotrophic_siteyears_daily.csv")

riversubset<-subset(riverdaily, site_name == "nwis_14211010")

riversubset$date<-as.Date(riversubset$date)

riversubset2<-riversubset[riversubset$date > "2010-01-01" & riversubset$date < "2011-01-01", ]


ERGPP<-ggplot(data=riversubset2) +   geom_line(aes(x = date, y = GPP), color="green")+  geom_line(aes(x = date, y = ER), color = "brown")+ facet_wrap(~ site_name, ncol = 5)

###dates###

riversubset2$Month <- format(as.Date(riversubset2$date), "%m")
riversubset2$Day <- format(as.Date(riversubset2$date), "%d")

riversubset2$jday <- format(riversubset2$date, "%j")


##GAM####

riversubset2$jday<-as.numeric(riversubset2$jday)

riverGAM<- gam(GPP~s(jday,  k = 50, m=2), data = riversubset2, method = "REML")
plot(riverGAM, shade = TRUE)


fd_GS <- derivatives(met_modGS, term = "site")

gam.check(riverGAM)

###derivative of GAM##
derGAM<- derivatives(riverGAM)

###plot data and GAM##

GPPGAM<-ggplot(data=riversubset2, aes(x = jday, y = GPP)) +   geom_line(color="green")+geom_smooth(method = 'gam', formula = y~s(x,  k = 50, m=2))+ggtitle("GPP CLACKAMAS RIVER")

GAMderivative<-ggplot(data=derGAM, aes(x = data, y = derivative)) +   geom_line(color="black")+ggtitle("GPP DERIVATIVE CLACKAMAS RIVER")

