library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(stringr) 

#function to transform name in lipid species, TC and DB

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


#read batch

batch <- read_tsv("Data/CD_results_Batch02._final_library.csv", col_types = cols(Name = col_character())) 
#key <- read_csv("Data/masslist_for_CD.csv", col_types = cols(Name = col_character()))   %>% 
#   select(Formula, species) 
  # unique()

# run
CDdata_gathered <- batch %>% 
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
         pqn = Area2/quotient.median)

var <-  unique(CDdata_gathered$class)
#[1]Batch 1 "Cer"   "CerPE" "CL"    "DG"    "LPC"   "PA"    "PC"    "PE"    "PG"    "PI"    "PS"    "TG" 
#[1]Batch 2 "PC"    "TG"    "PE"    "CL"    "PS"    "CerPE" "PG"    "PA"    "LPC"   "PI"    "Cer"   "DG" 
# sum the areas belonging to the same class of lipids (i.e. the areas of the different species)
# batch_0420 is the final data for statistical analysis of ALL lipids based on normalised areas

batch_0520 <-  CDdata_gathered  %>% 
  group_by(class, Sample) %>%  
  mutate(sumpqn = sum(pqn))
 # mutate(population = tolower(population))

#write_csv(batch01_0420, path = "Data/batch_____.csv") # file saved :)

# plot batch 02 all lipids minus plasmalogens

batch_0520_plot <-  ggplot(batch_0520, aes(`class`, log(sumpqn)))+
  geom_boxplot(aes(fill = `population`), width =1) +
  scale_fill_brewer(palette = "Paired") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch 02 Lipid Classes", y = "Log normalised area", x = "Lipid Classes") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11)) 

batch_0520_plot

# plots using the for loop 

library(ggplot2)
#lapply(var, function(x))

df <- batch_0520
var <-  unique(df$class)


for (i in var) {
      plot =  ggplot(data= subset(df, class == i)) +  
          geom_boxplot(aes(x= species, y= log(pqn), fill = population )) +
          ggtitle (i) +
          theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.text = element_text(size = 10, vjust = -0.25), 
          axis.text.x = element_text(angle = 40))
      
      print(plot)
}

# ----------------------------------------------
# old long version for plots here below

df <- batch01_0520

species_plot<- function(x){
  df <- df %>% 
    filter(class == x)
  # mutate(C = str_extract(species, "\\d+(?=\\:)"),
  #        DB = str_extract(species, "(?<=\\:)\\d+"))
  #plotcolour <- plotcolour
  ggplot(df, aes(species, log(pqn))) +
    geom_boxplot(aes(fill = population))+
    labs(title = x)+
    theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.text = element_text(size = 10, vjust = -0.25), 
          axis.text.x = element_text(angle = 40))
}

species_plot( "TG")
species_plot( "CL")
species_plot( "DG")
species_plot( "PG")
species_plot( "PE")
species_plot( "PI")
species_plot( "PC")
species_plot( "PS")
species_plot( "LPC")
species_plot( "PA")
species_plot( "Cer")
species_plot( "CerPE")
species_plot( "CerP")
species_plot( "SM")
species_plot( "SPH")
species_plot( "FA")

# -------------------------------------------

