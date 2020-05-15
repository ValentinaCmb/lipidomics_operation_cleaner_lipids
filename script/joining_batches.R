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


batch_01 <- read_tsv("Data/CD_results_Batch01._final_library.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT01 = RT)%>% 
  mutate(batch = "batch_01")

batch_02 <- read_tsv("Data/CD_results_Batch02._final_library.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT02 = RT)%>% 
  mutate(batch = "batch_02")

batch_03 <- read_tsv("Data/CD_results_Batch03._final_library.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT03 = RT)%>% 
  mutate(batch = "batch_03")
batch_04 <- read_tsv("Data/CD_results_Batch04._final_library.csv", col_types = cols(Name = col_character())) %>%
  # rename(RT04 = RT)%>%  
  mutate(batch = "batch_04")


CDdata_gathered_03 <- batch_03 %>% 
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
  filter(count >= 5) %>%
  mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(population = tolower(population)) %>% 
  filter(!is.na(population)) %>% 
  filter(!population == "syd_160") %>% # for batch03 only
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`, Sample, Area, count) %>% 
  mutate(class = str_extract(`Name`, "[a-zA-Z]+"),
         species = extract_TC(Name)) %>% 
  # mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`species`, 1, 1, sep = ":"),
         DB = word(`species`, -1, -1, sep = ":")) %>% 
  #!!  mutate(TC = str_extract(TC, "\\d+(?=\\:)"), #TC and DB instead of species
  #       DB = str_extract(DB, "(?<=\\:)\\d+"))
  # C = sum(as.numeric(unlist(str_extract_all(Name, "\\d+(?=\\:)")))),
  # DB = sum(as.numeric((unlist(str_extract_all(Name, "(?<=\\:)\\d+")))))
  # mutate(TC = word(`Name`, 1, 1, sep = "_"),
  #        DB = word(`Name`, -1, -1, sep = "_"))
  # filter(!grepl("e|p|Q|O", `Lipid species`)) %>% # this excludes the plasmalogen but not Ceramides and SPHM (kept by deleting d and t from grepl)
  # mutate(DB = as.numeric(str_extract(DB, "\\d+"))) %>%   #this separate the letters of plasmalogen from the DB IF run
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
  mutate(batch = "batch_03")



CDdata_gathered_04 <- batch_04 %>% 
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
  filter(count >= 5) %>%
  mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(population = tolower(population)) %>% 
  filter(!is.na(population)) %>% 
  #filter(!population == "syd_160") %>% # for batch03 only
  group_by(Name, Formula, `RT [min]`, `Area (Max.)`, Sample, Area, count) %>% 
  mutate(class = str_extract(`Name`, "[a-zA-Z]+"),
         species = extract_TC(Name)) %>% 
  # mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`species`, 1, 1, sep = ":"),
         DB = word(`species`, -1, -1, sep = ":")) %>% 
  #!!  mutate(TC = str_extract(TC, "\\d+(?=\\:)"), #TC and DB instead of species
  #       DB = str_extract(DB, "(?<=\\:)\\d+"))
  # C = sum(as.numeric(unlist(str_extract_all(Name, "\\d+(?=\\:)")))),
  # DB = sum(as.numeric((unlist(str_extract_all(Name, "(?<=\\:)\\d+")))))
  # mutate(TC = word(`Name`, 1, 1, sep = "_"),
  #        DB = word(`Name`, -1, -1, sep = "_"))
  # filter(!grepl("e|p|Q|O", `Lipid species`)) %>% # this excludes the plasmalogen but not Ceramides and SPHM (kept by deleting d and t from grepl)
  # mutate(DB = as.numeric(str_extract(DB, "\\d+"))) %>%   #this separate the letters of plasmalogen from the DB IF run
  group_by(Sample) %>%
  mutate(Area2 = Area/sum(Area) * 1000000) %>% 
  ungroup() %>% 
  group_by(Formula, `Molecular Weight`, `RT [min]`) %>%
  mutate(reference = median(Area2),
         quotient = Area2/reference) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(quotient.median = median(quotient),
         pqn = Area2/quotient.median) %>% 
  mutate(batch = "batch_04")



batch_final03 <-  CDdata_gathered_03  %>% 
  group_by(class, Sample) %>%  
  mutate(sumpqn = sum(pqn))

batch_final04 <-  CDdata_gathered_04  %>% 
  group_by(class, Sample) %>%  
  mutate(sumpqn = sum(pqn))


# FULL JOIN here

batch_3_4 <- full_join(batch_final03, batch_final04)
#   full_join(batch_03) %>% 
#   full_join(batch_04) %>% 
#   select(-"X11", -Window)
# mutate(population = tolower(population))

#write_csv(batch01_0420, path = "Data/batch_____.csv") # file saved :)

# plot batch 02 all lipids minus plasmalogens

batch_0304_plot <-  ggplot(batch_3_4, aes(`class`, log(sumpqn)))+
  geom_boxplot(aes(fill = `population`), width =1) +
  facet_wrap(~batch) +
  scale_fill_brewer(palette = "Paired") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch 03 and 04 Lipid Classes", y = "Log normalised area", x = "Lipid Classes") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11)) 

batch_0304_plot

# plots using the for loop 





 

# batch_bind <- batch_01 %>% 
#   rbind(batch_02, batch_03, batch_04) %>% 
#   mutate(RT = round(RT, 1)) %>% 
#   unique()
# 
# # prepare  the data to make a library for CD that contains only Name and Formula
# batchall <- batch_bind %>% 
#   select(-"X11", -Window, -batch) %>% 
#   unique() %>% 
#   arrange(CompoundName, RT) %>% 
#   rename(Name = CompoundName)
# 
# CD_library <- batchall %>% 
#   rename(Name = CompoundName) %>% 
#   rename( Formula = ChemicalFormula) %>% 
#   select(Name, Formula) %>% 
#   unique()

write_csv(CD_library, path = "data/CD_library_noRT.csv")




# over_rep <- batchall %>% 
#   group_by(CompoundName, ChemicalFormula, Ionization, Polarity, Adduct, ChargeState, ExtractedMass) %>% 
#   filter(RT < 25,
#          ChargeState == 1) %>% 
#   summarise(
#     n = n(),
#     meanRT = mean(RT),
#     sdRT = sd(RT),
#     minRT = min(RT),
#     maxRT = max(RT)
#   )
# 
# View(filter(over_rep, n > 4))
# 
# 
# ggplot(over_rep, aes(meanRT, ExtractedMass, colour = Polarity)) +
#   geom_point()
# 
# View(filter(over_rep, ExtractedMass > 1300))


