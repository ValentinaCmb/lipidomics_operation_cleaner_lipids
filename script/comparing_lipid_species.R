library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(stringr)


batch_01 <- read_tsv("Data/CD_results_b01_unfiltered.csv", col_types = cols(Name = col_character())) 

batch_02 <- read_tsv("Data/CD_results_b02_unfiltered.csv", col_types = cols(Name = col_character())) 

batch_03 <- read_tsv("Data/CD_results_b03_unfiltered.csv", col_types = cols(Name = col_character())) 

batch_04 <- read_tsv("Data/CD_results_b04_unfiltered.csv", col_types = cols(Name = col_character()))


B01 <- batch_01 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_01 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)

B02 <- batch_02 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_02 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)


B12 <- inner_join(B01, B02, by= "Name")
B12_not_common <- anti_join(B01, B02, by= "Name")

# 

B03 <- batch_03 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_03 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)

B04 <- batch_04 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_04 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)


B34 <- inner_join(B03, B04, by= "Name")
B34_not_common <- anti_join(B03, B04, by= "Name")

B12vs34_common <- inner_join(B12, B34, by= "Name")
B12vs34_not_common <- anti_join(B12, B34, by= "Name")
