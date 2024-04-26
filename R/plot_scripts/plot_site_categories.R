
library(tidyverse)
library(ggpattern)



lit <- read_csv('data_working/literature_streams_data_for_PCA.csv') %>%
  mutate(across(c('Type1', 'Type2'),
                ~case_when(. == 'Desert' ~ 'd. Desert',
                           . == 'Urban' ~ 'e. Urban', 
                           . == 'Agricultural' ~ 'f. Agricultural',
                           . == 'Snowmelt' ~ 'c. Snowmelt', 
                           . == 'Spring' ~ 'a. Spring', 
                           . == 'Dam tail waters' ~ 'b. Dam tail \nwaters',
                           TRUE ~ NA_character_))) 

dat <- read_csv('data_working/Autotrophic_site_classifications.csv')  %>%
  mutate(across(c('primary_class', 'secondary_class'),
                ~case_when(. == 'desert' ~ 'd. Desert',
                           . == 'urban' ~ 'e. Urban', 
                           . == 'agri' ~ 'f. Agricultural',
                           . == 'snowmelt' ~ 'c. Snowmelt', 
                           . == 'dam' ~ 'b. Dam tail \nwaters',
                           TRUE ~ NA_character_))) %>%
  rename(Type1 = primary_class, Type2 = secondary_class) %>%
  filter(!is.na(Type1))


lit_types <-data.frame(type = c(lit$Type1, lit$Type2)) %>%
  group_by(type) %>%
  summarize(lit_count = n()) %>%
  mutate(Literature = lit_count/nrow(lit)*100)
dat_types <-data.frame(type = c(dat$Type1, dat$Type2)) %>%
  group_by(type) %>%
  summarize(dat_count = n(),
            'This Study' = dat_count/nrow(dat)*100)
  

types <- left_join(lit_types, dat_types, by = 'type') %>%
  filter(!is.na(type)) %>%
  select(-ends_with('count')) %>%
  mutate(type = as.factor(type),
         type = fct_reorder(type, desc(type))) %>%
  pivot_longer(c('Literature', 'This Study'), 
               names_to = "category", values_to = 'percent') %>%
  mutate(category = factor(category, levels = c('This Study', 'Literature' )),
         patt = case_when(category == 'This Study' ~ 'none', 
                          category == 'Literature' ~ 'stripe'))
  
# Barplot

png('figures/site_categorization_barplot.png', 
    width = 2, height = 5, units = 'in', res = 300)
  ggplot(types, aes(x=type, y=percent, col = category, fill = category)) + 
    # geom_bar(stat = "identity",position="dodge") +
    geom_col_pattern(
      aes(pattern = patt, col = category, fill = category 
          ), pattern_fill = 'darkgreen', pattern_colour = 'darkgreen',
      position = 'dodge', width = 0.7) +
    scale_pattern_identity() +
    coord_flip() +
    scale_color_discrete("", type = c('forestgreen', 'darkgreen'))+
    scale_fill_discrete("", type = c('#2DB60080', 'white'))+
    theme_classic() +
    xlab('') + ylab('Percent \nof Sites')+
    theme(legend.position='top',
          legend.title = element_text(size = 8), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.6, "lines")) +
    guides(fill=guide_legend(nrow=2),
           col = guide_legend(nrow = 2)) 
  
dev.off()
