library(tidyverse)

batch_qcs <- read_csv("Data/CD_QCs_analysis.csv", col_types = cols(Name = col_character())) 
  
tidier_batch <- batch_01 %>% 
  rbind(batch_02, batch_03, batch_04) %>% 
  mutate(RT = round(RT, 1)) %>% 
  select(-"X11", -Window) %>% 
  arrange(CompoundName, RT) %>% 
  rename(Name = CompoundName, Formula = ChemicalFormula) %>% 
  unique() %>% 
  mutate(Class = str_extract(`Name`, "[a-zA-Z]+"))  
# mutate(TC = word(`Name`, 1, 1, sep = "_"),
#        DB = word(`Name`, -1, -1, sep = "_"))


cleaned_batch <- tidier_batch %>% 
  mutate(Reject = case_when(Adduct == "Na-Gain" ~ "rej",
                            Class == "CL" & RT < 10 ~ "rej",
                            TRUE ~ ""))