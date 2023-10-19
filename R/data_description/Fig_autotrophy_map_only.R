# Annual autotrophy in streams
# May 26, 2022
# J. Blaszczak, C. Barbosa, M. Desiervo, H. Lowman

## load more packages
lapply(c("lubridate","cowplot",
         "tidyverse", "reshape2",
         "ggmap","maps","mapdata", "mapproj",
         "ggsn","wesanderson", "sf"), require, character.only=T)

## Import Data
GPP_info <- read_csv("data_working/across_sites_model_data.csv") %>%
  mutate(logPtoP = log(PR))

sites <- read_tsv("data_356rivers/site_data.tsv")

# First, check to see if there are any sites in AK, HI, or PR.
site_name <- unique(GPP_info$site_name)

my_sites <- as.data.frame(site_name)

check_sites <- left_join(my_sites, sites) # None.

# Calculate mean annual GPP and mean annual P:R.
Annual_info <- GPP_info %>%
  group_by(site_name, lat, lon) %>%
  summarize(mean_GPP_annual = mean(ann_GPP_C),
            mean_ER_annual = mean(ann_ER_C),
            mean_NEP_annual = mean(ann_NEP_C)) %>%
  ungroup() %>%
  mutate(mean_PR_annual = mean_GPP_annual/(-mean_ER_annual))

# Also identify the autotrophic sites (n = 15).
aut_sites <- GPP_info %>% 
  select(site_name, ER = ann_ER_C, GPP = ann_GPP_C, PR, NEP = ann_NEP_C) %>%
  group_by(site_name) %>%
  # Using Alice C's workflow to calculate median PR/NEP across all site-years
  summarize(across(everything(), median, na.rm = T)) %>%
  ungroup() %>%
  filter(NEP > 0)

aut_site_names <- unique(aut_sites$site_name)

Annual_info$class <- case_when(Annual_info$site_name %in% aut_site_names ~ "Aut",
                               TRUE ~ "Het")

## Generate map of NEP
(NEP_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                             maptype='toner'))+
    geom_point(data = Annual_info, aes(x = lon, y = lat, 
                                  fill=mean_NEP_annual, size=mean_NEP_annual), shape=21)+
    theme(legend.position = "right")+
    labs(x="Longitude", y="Latitude")+
    scale_fill_gradient("Mean annual NEP (units C)",
                        low = "blue", high = "red")+
    scale_size_continuous("Mean annual NEP (units C)"))

## Generate map of P:R
(PR_map_fig <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = Annual_info, aes(x = lon, y = lat, 
                                fill=mean_PR_annual, size=mean_PR_annual), shape=21)+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")+
  scale_fill_gradient("P:R",
                      low = "brown", high = "green")+
  scale_size_continuous("P:R"))

## Generate map sized by mean annual GPP and colored by P:R.
# make data sf object
AI_sf <- st_as_sf(Annual_info,
                    coords = c("lon", "lat"), # always put lon (x) first
                    remove = F, # leave the lat/lon columns in too
                    crs = 4269) # projection: NAD83

# make base US map using code from https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf
states <- ggplot2::map_data("state")
state_sf <- st_as_sf(states,
                     coords = c("long", "lat"),
                     remove = F,
                     crs = 4269)

# Set color legend break points.
my_breaks <- c(0.05, 0.15, 0.4, 1)

# Map with points colored on a scale of P:R
(sitemap <- ggplot(state_sf) + # base plot
    geom_polygon(aes(x = long, y = lat, group = group), 
                 fill = "white", color = "black") + # map of states
    geom_point(data = AI_sf, aes(x = lon, y = lat,
                                 color = mean_PR_annual, 
                                 size = mean_GPP_annual),
               alpha = 0.95) + # map of sites
    labs(x = "Longitude", y = "Latitude") +
    scale_color_gradientn("Mean Annual P:R", colors = c("#F0C9C0", "#EDB48E", 
                                             "#E8C32E", "#E6E600", "#00A600"),
                                             trans = "log",
                          breaks = my_breaks, labels = my_breaks) +
    # based on terrain.colors(n = 10) # custom colors
    scale_size_continuous("Mean Annual\nCumulative GPP") +
    theme_classic()) # remove grid

(sitemap_binary <- ggplot(state_sf) + # base plot
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") + # map of states
  geom_point(data = AI_sf, aes(x = lon, y = lat,
                               #fill = mean_PR_annual, 
                               size = mean_GPP_annual,
                               color = class,
                               stroke = class),
             shape = 21, alpha = 0.95) + # map of sites
  labs(x = "Longitude", y = "Latitude") +
  # all points will be on same P:R color scale
  # scale_fill_gradientn("Mean Annual P:R", colors = c("#F0C9C0", "#EDB48E", 
  #                                         "#E8C32E", "#E6E600", "#00A600"),
  #                                         trans = "log",
  #                       breaks = my_breaks, labels = my_breaks) +
                        # based on terrain.colors(n = 10) # custom colors
  scale_size_continuous("Mean Annual\nCumulative GPP") +
  # autotrophic sites will be in green, het will be pink
  scale_color_manual("Category", values = c("#00A600", "#EDB48E")) +
  # autotrophic sites will have a thicker outline
  scale_discrete_manual(aesthetics = "stroke", values = c(2, 1)) +
  theme_classic()) # remove grid

# export exploratory figures
# ggsave(("figures/Annual_GPP_PR_USmap_082123.png"),
#        width = 22,
#        height = 11,
#        units = "cm"
# )

# End of script.
