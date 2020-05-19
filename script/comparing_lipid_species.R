library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(stringr)


batch_01 <- read_tsv("Data/CD_results_b01_unfiltered.csv", col_types = cols(Name = col_character())) 
batch_01_tags <- read_tsv("Data/CD_results_B01_extraTAGs.csv", col_types = cols(Name = col_character()))
# tags indiccate data set wit extra HH3 ion added to ID tags

batch_02 <- read_tsv("Data/CD_results_b02_unfiltered.csv", col_types = cols(Name = col_character())) 

batch_03 <- read_tsv("Data/CD_results_b03_unfiltered.csv", col_types = cols(Name = col_character())) 

batch_04 <- read_tsv("Data/CD_results_b04_unfiltered.csv", col_types = cols(Name = col_character()))

# compariosn bt B01 and B01+ extra tags

B01 <- batch_01 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_01 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)

B01_tags <- batch_01_tags %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, RT_02 = `RT [min]`) %>% 
  distinct(Name, .keep_all=TRUE)


B11tags <- inner_join(B01, B01_tags, by= "Name") 
B1tags1_not_common <- anti_join(B01_tags, B01, by= "Name")

# comparison btw lipid species in B1 and 2

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

# comparison btw lipid species in B3 and 4

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

#comparison btw lipid species btw subbatches 12 and 34

B12vs34_common <- inner_join(B12, B34, by= "Name")
B12vs34_not_common <- anti_join(B12, B34, by= "Name")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    