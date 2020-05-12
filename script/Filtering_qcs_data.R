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

  group_by(Name, `RT [min]`)

    
#or
  CDdata_gathered <- batch_qcs %>% 
    gather(Sample, Area, contains("Area: ")) %>%
    filter(!is.na(Area)) %>%
    group_by(Name, Formula, `RT [min]`, `Molecular Weight`,`Area (Max.)`) %>% 
    summarise(
      count = n()
    ) %>% 
    mutate(Class = str_extract(`Name`, "[a-zA-Z]+"))
  
  CD_filtered <- CDdata_gathered %>%  
    filter(count >= 10)
  
    
   ggplot(CD_filtered, aes(`RT [min]`, `Molecular Weight`)) +
     geom_point()
   