library(tidyverse)

batch_qcs <- read_tsv("Data/CD_QCs_analysis.csv", col_types = cols(Name = col_character())) 

CDdata_gathered <- batch_qcs %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  unique() %>% 
    filter(!is.na(Area)) %>% 
  select("Name", 'Formula','RT [min]', 'Area (Max.)', 'Sample') %>% 
  group_by(Name, `RT [min]`) %>% 
  mutate(count = n()) %>% 
  filter(count >= 10)


  unique('Sample')

  group_by(Name, `RT [min]`)%>%
