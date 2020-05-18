library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(stringr)

extract_TC <- function(molecule){
  lipidclass <- str_extract(molecule, "[a-zA-Z]+")
  # df <- as.data.frame(str_extract_all(molecule, "\\d+[:]\\d+")) %>%
  #   separate(1, c("C", "DB"), ":", convert = T) %>%
  #   summarise(
  #     Carbons = sum(C),
  #     DB = sum(DB)
  #   )
  C <- sum(as.numeric(unlist(str_extract_all(molecule, "\\d+(?=\\:)"))))
  DB <- sum(as.numeric((unlist(str_extract_all(molecule, "(?<=\\:)\\d+")))))
  modification <- str_extract(str_extract(molecule, "(?<=\\().+?(?=\\))"), "[a-z]")
  # species <- ifelse(modification == "t"|modification == "d", paste(modification, df$Carbons, ":", df$DB, sep = ""), 
  #                   ifelse(modification == "e"| modification == "p", paste(df$Carbons, modification, ":", df$DB, sep = ""), 
  #                   paste(df$Carbons, df$DB, sep = ":")))
  
  species <- case_when(modification == "t"|modification == "d" ~ paste(modification, C, ":", DB, sep = ""), 
                       modification == "e"| modification == "p" ~  paste(C, ":", DB, modification, sep = ""), 
                       TRUE ~ paste(C, DB, sep = ":"))
  return(species)
}


batch_01 <- read_tsv("Data/CD_results_b01_unfiltered.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT01 = RT)%>% 
  mutate(batch = "batch_01")

batch_02 <- read_tsv("Data/CD_results_b02_unfiltered.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT02 = RT)%>% 
  mutate(batch = "batch_02")

B01 <- batch_01 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, Formula, `Molecular Weight`, `RT [min]`, `Area (Max.)`, Sample, Area) %>% 
  unique() %>% 
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`) %>% 
  filter(!is.na(Area),
         `RT [min]` < 25,
         !is.na(Name)) %>% 
  mutate(
    count = n()
  ) %>%  
  mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(population = tolower(population)) %>% 
  filter(!is.na(population)) %>% 
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`, Sample, Area, count) %>% 
  mutate(class = str_extract(`Name`, "[a-zA-Z]+"),
         species = extract_TC(Name)) %>% 
  # mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`species`, 1, 1, sep = ":"),
         DB = word(`species`, -1, -1, sep = ":")) %>% 
    group_by(Sample) %>%
  mutate(Area2 = Area/sum(Area) * 1000000) %>% 
  ungroup() %>% 
  group_by(Formula, `Molecular Weight`, `RT [min]`) %>%
  mutate(reference = median(Area2),
         quotient = Area2/reference) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(quotient.median = median(quotient),
         pqn = Area2/quotient.median)%>% 
  mutate(batch = "batch_01") %>% 
  filter((class =="CL" & !`RT [min]` < 15) | (class =="TG" & !`RT [min]` < 10)) %>% 
  ungroup() %>% 
  # group_by(class, Sample) %>%  
  # mutate(sumpqn = sum(pqn)) %>% 
  # ungroup() %>% 
  select(Name, RT_01 = `RT [min]`, population, sample, class, pqn, batch)

B02 <- batch_02 %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  select(Name, Formula, `Molecular Weight`, `RT [min]`, `Area (Max.)`, Sample, Area) %>% 
  unique() %>% 
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`) %>% 
  filter(!is.na(Area),
         `RT [min]` < 25,
         !is.na(Name)) %>% 
  mutate(
    count = n()
  ) %>%  
  mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(population = tolower(population)) %>% 
  filter(!is.na(population)) %>% 
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`, Sample, Area, count) %>% 
  mutate(class = str_extract(`Name`, "[a-zA-Z]+"),
         species = extract_TC(Name)) %>% 
  # mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`species`, 1, 1, sep = ":"),
         DB = word(`species`, -1, -1, sep = ":")) %>% 
  group_by(Sample) %>%
  mutate(Area2 = Area/sum(Area) * 1000000) %>% 
  ungroup() %>% 
  group_by(Formula, `Molecular Weight`, `RT [min]`) %>%
  mutate(reference = median(Area2),
         quotient = Area2/reference) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(quotient.median = median(quotient),
         pqn = Area2/quotient.median)%>% 
  mutate(batch = "batch_01") %>% 
  filter((class =="CL" & !`RT [min]` < 15) | (class =="TG" & !`RT [min]` < 10)) %>% 
  ungroup() %>% 
  # group_by(class, Sample) %>%  
  # mutate(sumpqn = sum(pqn)) %>% 
  # ungroup() %>% 
  select(Name, RT_01 = `RT [min]`, population, sample, class, pqn, batch)

B12 <- inner_join(B01, B02, by= "Name")
