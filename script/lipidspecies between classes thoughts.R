# test for confirming lipid species in two batches
library(tidyverse)
batch_01 <- read_tsv("Data/CD_results_Batch01._final_library.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT01 = RT)%>% 
  mutate(batch = "batch_01")

batch_02 <- read_tsv("Data/CD_results_Batch02._final_library.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT02 = RT)%>% 
  mutate(batch = "batch_02")

b01 <- batch_01 %>% 
  select(Name, RT_01 = `RT [min]`)
b02 <- batch_02 %>% 
  select(Name, RT_02 = `RT [min]`)

full_batch <- full_join(b01, b02) %>% 
  group_by(Name)
  
outer_batch <- anti_join(b01, b02)

bind_batch <- rbind(select(batch_01, Name, `RT [min]`, batch), select(batch_02, Name, `RT [min]`, batch)) %>%
  select(Name, batch, RT = `RT [min]`) %>% 
  group_by(Name) %>% 
  summarise(
    n = n(),
    meanRT = mean(RT),
    sdRT = sd(RT)
  )
  

comm <- na.omit(bind_batch)
