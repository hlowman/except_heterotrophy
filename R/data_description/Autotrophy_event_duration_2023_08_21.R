# Autotrophy Event Duration
# August 21, 2023
# Joanna Blaszczak

######################################
## Autotrophy duration distribution
######################################

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse", "data.table"), require, character.only=T)

## Import lotic_standardized_full from Bernhardt metabolism 2022 data release
## Download from: https://figshare.com/articles/software/Code_and_RDS_data_for_Bernhardt_et_al_2022_PNAS_/19074140?backTo=/collections/Data_and_code_for_Bernhardt_et_al_2022_PNAS_/5812160
## Setwd to output_data
lotic_standardized_full <- readRDS("lotic_standardized_full.rds")

## Subset data frame for test
df <- lotic_standardized_full[1:3]

## Extract events and duration of events
duration_calc <- function(df){

  d <- df[,c("Site_ID","Date","GPP","ER")]
  d$NEP <- d$GPP - abs(d$ER)
  
  d <- na.omit(d)
  
  # First calc time difference and split to segments to avoid NA days
  d$diff_time <- NA
  d$diff_time[1] <- 0
  
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$Date[i], time2 = d$Date[(i-1)], units="days")
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
      mutate(NEP_above = NEP > t, id = rleid(NEP_above)) %>% 
      # keep only periods with autotrophy
      filter(NEP_above) %>%
      # for each period/event, get its duration
      group_by(id) %>%
      summarise(event_duration = difftime(last(Date), first(Date), units = "days"),
                start_date = first(Date),
                end_date = last(Date))
    
    zz[nrow(zz)+1,] <- NA
    
    return(zz)
  }
  
  event_above0 <- ldply(lapply(lseq, function(x) events_calc(x, 0)), data.frame);event_above0$NEP_thresh <- 0
  event_above0.5 <- ldply(lapply(lseq, function(x) events_calc(x, 0.5)), data.frame);event_above0.5$NEP_thresh <- 0.5
  event_above1 <- ldply(lapply(lseq, function(x) events_calc(x, 1)), data.frame);event_above1$NEP_thresh <- 1
  event_above5 <- ldply(lapply(lseq, function(x) events_calc(x, 5)), data.frame);event_above5$NEP_thresh <- 5
  
  events <- rbind(event_above0,
                  event_above0.5,
                  event_above1,
                  event_above5)
  
  ## subset
  events_df <- events[,c("event_duration","start_date","end_date","NEP_thresh")]
  events_df$SiteID <- d$Site_ID[1]
  events_df <- na.omit(events_df)
  
  return(events_df)
  
}

auto_events <- lapply(df, function(x) duration_calc(x))
auto_df <- ldply(auto_events, data.frame)
head(auto_df);tail(auto_df)

## Add 1 to event duration
auto_df$event_duration <- auto_df$event_duration+1
auto_df$event_dur <- as.numeric(auto_df$event_duration)

## Visualize
ggplot(auto_df, aes(event_dur, fill=NEP_thresh))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~NEP_thresh,ncol=1)+
  theme_bw()

saveRDS(auto_df, "data_working/autotrophic_event_durations.rds")

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

fig1 <- ggplot(auto_df, aes(quant_val, fill=as.factor(NEP_thresh)))+
  geom_bar(alpha=0.4, color="black", position="identity")+
  theme_bw()+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=14, angle=35, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.position = "top")+
  labs(x="Event duration", y="Number of events")
fig1

# ggsave(("figures/auto_events_duration.png"),
#        width = 25,
#        height = 15,
#        units = "cm"
# )

#############################
## What month is the onset?
#############################

auto_df$month <- month(auto_df$start_date)

fig2 <- ggplot(auto_df, aes(as.factor(month)))+
  geom_bar(alpha=0.4, color="black", position="identity")+
  facet_wrap(~as.factor(quant_val), ncol=1, scales = "free_y")+
  theme_bw()
fig2

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
