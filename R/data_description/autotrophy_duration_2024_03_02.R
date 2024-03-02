##==============================================================================
## Script for autotrophy duration 
## Code author: J.R. Blaszczak
## Last Edited: March 1, 2024
##
## Changed NEP to P:R as the metric by which we are quanitfying autotrophy
##==============================================================================

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse"), require, character.only=T)

## import most up-to-date dataset
dat <- readRDS("../../data_356rivers/high_quality_daily_metabolism_with_SP_covariates.rds")
df <- dat[,c("site_name","date","GPP","ER")]
df <- na.omit(df)

## split to list by site
df_list <- split(df, df$site_name)

###########################################
## Extract events and duration of events
#############################################
duration_calc <- function(df_site){

  d <- df_site[,c("site_name","date","GPP","ER")]
  d$PR <- d$GPP/abs(d$ER)

  # First calc time difference and split to segments to avoid NA days
  d$diff_time <- NA
  d$diff_time[1] <- 0
  
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$date[i], time2 = d$date[(i-1)], units="days")
  }
  
  d$diff_time <- as.character(as.numeric(d$diff_time))
  d$seq <- NA
  d$seq[1] <- 1
  
  for(i in 2:nrow(d)){
    if(d$diff_time[i] %in% c("1")){
      d$seq[i] = d$seq[(i-1)]
    } else{
      d$seq[i] = d$seq[(i-1)]+1
    }
  }
  
  lseq <- split(d, as.factor(d$seq))
  
  
  events_calc <- function(z, t) {
    zz <- z %>% 
      #add id for different periods/events
      mutate(PR_above = PR > t, id = data.table::rleid(PR_above)) %>% 
      # keep only periods with autotrophy
      filter(PR_above) %>%
      # for each period/event, get its duration & magnitude
      group_by(id) %>%
      summarise(event_duration = difftime(last(date), first(date), units = "days"),
                start_date = first(date),
                end_date = last(date),
                PR_mean = mean(PR))
    
    zz[nrow(zz)+1,] <- NA
    
    return(zz)
  }
  
  event_above1 <- ldply(lapply(lseq, function(x) events_calc(x, 1)), data.frame);event_above1$PR_thresh <- 1
  
    ## subset
  events_df <- event_above1[,c("event_duration","start_date","end_date","PR_mean","PR_thresh")]
  events_df$site_name <- d$site_name[1]
  events_df <- na.omit(events_df)
  
  return(events_df)
  
}


auto_events <- lapply(df_list, function(x) duration_calc(x))
auto_df <- ldply(auto_events, data.frame)
## check event duration
auto_df[which(auto_df$event_duration < 0),]

## Add 1 to event duration
auto_df$event_duration <- auto_df$event_duration+1
auto_df$event_dur <- as.numeric(auto_df$event_duration)
head(auto_df);tail(auto_df)

## Visualize
ggplot(auto_df, aes(event_dur, fill=PR_thresh))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~PR_thresh,ncol=1)+
  theme_bw()

## save
saveRDS(auto_df, "../../data_356rivers/autotrophic_event_duration_PR.rds")

#############################################
## Calculations for results
#############################################
auto_df <- readRDS("../../data_356rivers/autotrophic_event_duration_PR.rds")

## 1 ## What % of rivers experienced at least one autotrophic event
length(levels(as.factor(auto_df$site_name))) # 212 sites
length(levels(as.factor(df$site_name))) # 223 sites
212/223 #95%

## 2 ## what percentage of days are autotrophic
sum(auto_df$event_dur) #46,042
sum(auto_df$event_dur)/nrow(df) #19%

## 3 ## Number of events per year
# Compare mean +/- events per year between 1-3 days and 1-3 months
#first compare if any events cross years
auto_df$year_start <- year(auto_df$start_date)
auto_df$year_end <- year(auto_df$end_date)
auto_df$year_diff <- auto_df$year_end - auto_df$year_start
nrow(auto_df[which(auto_df$year_diff > 0),]) ## only 18 sites; will attribute to start year
#classify if an event is 1-3 days or 1-3 months
duration_days<-c(1, 3, 7, 14, 30, 90)
auto_df$duration_cat <- factor(findInterval(auto_df$event_dur,duration_days))
auto_df$duration_length <- revalue(auto_df$duration_cat, c("1" = "1 day to 3 days",
                                              "2" = "3 days to 1 week",
                                              "3" = "1 week to 2 weeks",
                                              "4" = "2 weeks to 1 month",
                                              "5" = "1 month to 3 months"))
#group by site and year and calculate mean duration_length category per year per site
length_site_year <- auto_df %>%
  group_by(site_name, year_start,duration_length) %>%
  count()
#expand this to include 0 for each length and site and year
#only include 1-3 days and 1-3 months
length(levels(as.factor(auto_df$site_name))) #212
levels(as.factor(length_site_year$year_start)) #2008-2016
# Create a site-year index
length_site_year$site_year <- paste(length_site_year$site_name,
                                    length_site_year$year_start,sep="_")
# For every site-year index, create a data frame with both 1-3 days and 1-3 months
events <- NULL
events <- as.data.frame(rep(levels(as.factor(length_site_year$site_year)),2))
colnames(events) <- "site_year"
events$duration_length <- c(rep(levels(as.factor(auto_df$duration_length))[1], 1616/2),
                                rep(levels(as.factor(auto_df$duration_length))[5], 1616/2))
head(events); tail(events)
#merge
combined <- merge(events, length_site_year,
                  by = c("site_year","duration_length"),
                  all = TRUE)
#change all NA n values to 0
combined[is.na(combined$n),]$n <- 0
#summarize
mean_sd_length_year <- combined %>%
  group_by(duration_length) %>%
  summarize(mean = mean(n), sd = sd(n))
# 1 to 3 days = 9.8 +/- 8.4 events/year
# 1 to 3 months = 0.04 +/- 0.21 events/year

## 4 ## Month of onset and termination
auto_df$onset_month <- month(auto_df$start_date)
auto_df$end_month <- month(auto_df$end_date)

ggplot(auto_df, aes(as.factor(onset_month)))+
    geom_bar(alpha=0.4, color="black", position="identity")+
    geom_bar(aes(as.factor(end_month)), alpha=0.4, color="blue", position="identity")+
    facet_wrap(~as.factor(duration_length), ncol=1, scales = "free_y")+
    theme_bw()

PC$start_hour <- hour(as.POSIXct(PC$start_date, format="%Y-%m-%d %H:%M:%S"))
PC$end_hour <- hour(as.POSIXct(PC$end_date, format="%Y-%m-%d %H:%M:%S"))

ggplot(PC_sub2[which(PC_sub2$quant_val %in% quant_val_sub_list),], aes(start_hour))+
  geom_bar(fill="#010D26", alpha=0.4, color="black")+
  geom_bar(aes(end_hour), fill="#4CBFBB", alpha=0.5, color="black")+
  labs(x="Hour of Day", y="Number of Events")+
  facet_wrap(~quant_val, nrow=1)+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        strip.background = element_rect(fill="white"))
time_of_day  

## All for SI
ggplot(PC_sub2, aes(start_hour))+
  geom_bar(fill="#010D26", alpha=0.4, color="black")+
  geom_bar(aes(end_hour), fill="#4CBFBB", alpha=0.5, color="black")+
  labs(x="Hour of Day", y="Number of Events")+
  facet_wrap(~quant_val, nrow=4)+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        strip.background = element_rect(fill="white"))






## 5 ## Magnitude of P:R during events








#############################################
## Which sites have long periods of NEP > 0
#############################################

# Read in dataset created above
auto_df <- readRDS("data_working/autotrophic_event_durations.rds")

auto_df[which(auto_df$event_dur > 30),]

## Group them
quantiles<-c(1, 3, 7, 14, 30, 90)
auto_df$quant <- factor(findInterval(auto_df$event_dur,quantiles))
auto_df$quant_val <- revalue(auto_df$quant, c("1" = "1 day to 3 days",
                                              "2" = "3 days to 1 week",
                                    "3" = "1 week to 2 weeks",
                                    "4" = "2 weeks to 1 month",
                                    "5" = "1 month to 3 months"))

## Plot
levels(factor(auto_df$NEP_thresh))
auto_df$NEP_thresh_name <- factor(auto_df$NEP_thresh, 
                                  levels = c("0" = "NEP > 0",
                                             "0.5" = "NEP > 0.5",
                                             "1" = "NEP > 1",
                                             "5" = "NEP > 5"))

(fig1 <- ggplot(auto_df, aes(quant_val, fill=as.factor(NEP_thresh)))+
  geom_bar(alpha=0.4, color="black", position="identity")+
  theme_bw()+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=14, angle=35, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.position = "top")+
  labs(x="Event duration", y="Number of events"))

# ggsave(("figures/auto_events_duration.png"),
#        width = 25,
#        height = 15,
#        units = "cm"
# )

#############################
## What month is the onset?
#############################

auto_df$month <- month(auto_df$start_date)

(fig2 <- ggplot(auto_df, aes(as.factor(month)))+
  geom_bar(alpha=0.4, color="black", position="identity")+
  facet_wrap(~as.factor(quant_val), ncol=1, scales = "free_y")+
  theme_bw())

# ggsave(("figures/auto_events_onset.png"),
#        width = 25,
#        height = 15,
#        units = "cm"
# )

############################
## Mean duration per site
###########################

auto_mean <- auto_df %>%
  group_by(SiteID, NEP_thresh) %>%
  summarize_at(.vars = "event_dur", .funs = mean)

auto_1 <- auto_mean[which(auto_mean$NEP_thresh == "1"),]

ggplot(auto_1, aes(event_dur))+
  geom_histogram()

## load more packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "plotrix", "data.table","ggmap","maps","mapdata",
         "ggsn","wesanderson"), require, character.only=T)

## merge with site_info
# data available here: https://www.sciencebase.gov/catalog/item/59bff64be4b091459a5e098b
# But file is small enough and has been added to "data_356rivers" folder
site_info <- read.table("data_356rivers/site_data.tsv",sep = "\t", header=T)
auto_1$site_name <- auto_1$SiteID
auto_1 <- merge(auto_1, site_info, by="site_name")

(fig3 <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = auto_1, aes(x = lon, y = lat, 
                                 fill=event_dur, size=event_dur),
             shape=21)+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")+
  scale_fill_gradient("Mean Autotrophic Event (days)",
                      low = "blue", high = "red",
                      breaks=c(1, 7, 14),
                      labels=c("1 day", "1 week", "2 weeks"))+
  scale_size_continuous("Mean Event Duration",
                        breaks = c(1,7,14),
                        labels=c("1 day", "1 week", "2 weeks")))

# ggsave(("figures/auto_events_USmap.png"),
#        width = 25,
#        height = 15,
#        units = "cm"
# )

# End of script.
