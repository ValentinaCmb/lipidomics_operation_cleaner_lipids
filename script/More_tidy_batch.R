library(tidyverse)

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

batch_01 <- read_csv("Data/LS_batch01.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT01 = RT)%>% 
  mutate(batch = "batch_01")

batch_02 <- read_csv("Data/LS_batch02.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT02 = RT)%>% 
  mutate(batch = "batch_02")

batch_03 <- read_csv("Data/LS_batch03.csv", col_types = cols(Name = col_character())) %>% 
  # rename(RT03 = RT)%>% 
  mutate(batch = "batch_03")
batch_04 <- read_csv("Data/LS_batch04.csv", col_types = cols(Name = col_character())) %>%
  # rename(RT04 = RT)%>%  
  mutate(batch = "batch_04")

# batch <- full_join(batch_01, batch_02) %>% 
#   full_join(batch_03) %>% 
#   full_join(batch_04) %>% 
#   select(-"X11", -Window) 

tidier_batch <- batch_01 %>% 
  rbind(batch_02, batch_03, batch_04) %>% 
  mutate(RT = round(RT, 1)) %>% 
  select(-"X11", -Window) %>% 
  arrange(CompoundName, RT) %>% 
  rename(Name = CompoundName, Formula = ChemicalFormula) %>% 
  unique() %>% 
  group_by(Name, RT, batch, Adduct) %>% 
  mutate(Class = str_extract(`Name`, "[a-zA-Z]+"),
         species = extract_TC(Name)
         # C = sum(as.numeric(unlist(str_extract_all(Name, "\\d+(?=\\:)")))),
         # DB = sum(as.numeric((unlist(str_extract_all(Name, "(?<=\\:)\\d+")))))
         )
        
  # mutate(TC = word(`Name`, 1, 1, sep = "_"),
  #        DB = word(`Name`, -1, -1, sep = "_"))

# plot based on batches
ggplot(tidier_batch, aes(RT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16))

#remove unwanted lipids
cleaned_batch <- tidier_batch %>% 
  mutate(Reject = case_when(RT > 25 ~ "rej",
                            ChargeState != 1 ~ "rej",
                            Adduct == "Na-Gain" ~ "rej",
                            Class == "CL" & RT < 15 ~ "rej",
                            Class == "TG" & RT < 10 ~"rej",
                         TRUE ~ ""))

# plot based on cleaned batches
ggplot(cleaned_batch %>% filter(Reject != "rej"), aes(RT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16))


#The adduct mass is reported - need to make this MW
H <- 1.007825
N <- 14.003074
electron <- 0.000548#0.000548579909

Adduct_spread <- cleaned_batch %>% 
  filter(Reject != "rej") %>% 
  select(Name, Class, Formula, batch, RT, ExtractedMass, Adduct) %>% 
  spread(Adduct, ExtractedMass) %>% 
  mutate(MW = case_when(!is.na(`H-Gain`) ~ `H-Gain` - H + electron,
                        !is.na(`H-Loss`) ~ `H-Loss` + H - electron,
                        !is.na(`NH4-Gain`) ~ `NH4-Gain` -(N + 4 * H) + electron,
                        !is.na(`HCO2-Gain`) ~ `HCO2-Gain` - 45),
         MW = round(MW, 4))

RT_test <- Adduct_spread %>% 
  group_by(Name, Class, Formula, MW) %>% 
  summarise(
    n = n(),
    meanRT = mean(RT),
    sdRT = sd(RT),
    minRT = min(RT),
    maxRT = max(RT)
  ) %>% 
  mutate(species = 
  extract_TC(Name))  

ggplot(RT_test, aes(meanRT, MW, colour = Class)) +
  geom_point()

RT_test%>% filter(sdRT>2)

Class_plot<- function(df, lipidClass, plotcolour = C){
  df <- df %>% 
    filter(Class == lipidClass) %>% 
    mutate(C = str_extract(species, "\\d+(?=\\:)"),
           DB = str_extract(species, "(?<=\\:)\\d+"))
  plotcolour <- plotcolour
  ggplot(df, aes(meanRT, MW, colour = DB)) +
    geom_point() +
    labs(title = lipidClass)
}

Class_plot(RT_test, "TG")
Class_plot(RT_test, "CL")
Class_plot(RT_test, "DG")
Class_plot(RT_test, "PG")
Class_plot(RT_test, "PE")
Class_plot(RT_test, "PI")
Class_plot(RT_test, "PC")
Class_plot(RT_test, "PS")
Class_plot(RT_test, "LPC")
Class_plot(RT_test, "PA")
Class_plot(RT_test, "Cer")
Class_plot(RT_test, "CerPE")
Class_plot(RT_test, "CerP")
Class_plot(RT_test, "SM")
Class_plot(RT_test, "SPH")
Class_plot(RT_test, "FA")


new_list <- RT_test %>% 
  ungroup() %>% 
  mutate(Formula = case_when(Class == "TG" ~ paste(Formula, "N", "H3", sep = " "),
                             Class == "DG" ~ paste(Formula, "N", "H3", sep = " "),
                             TRUE~Formula)) %>% 
  arrange(Class, species)

over_rep <- cleaned_batch %>% 
  group_by(Name, Class, Formula, Ionization, Polarity, Adduct, ChargeState, ExtractedMass, batch) %>% 
  filter(Reject != "rej") %>% 
  summarise(
    n = n(),
    meanRT = mean(RT),
    sdRT = sd(RT),
    minRT = min(RT),
    maxRT = max(RT)
  )

#View(filter(over_rep, n > 4))

# ????
over_rep02 <- over_rep %>% 
  ungroup() %>% 
  #group_by(Class) %>% 
  #mutate(batch = word(`batch`, 2,2, sep ="_")) %>% 
  spread(Adduct, ChargeState)

 
# plot based on Classes

lipid_class <- unique(over_rep02$Class)

class  <- over_rep02 %>% 
  filter(Class == "CL")

ggplot(class, aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16))

# trial to build the function
plots_classes <-  over_rep02 %>% 
  filter(Class == "SPH" ) %>%
 ggplot(aes(meanRT, ExtractedMass)) +
   geom_point(aes(colour = batch, shape = Polarity)) +
   scale_shape_manual(values = c(6,16)) +
  labs(title = "SPH") 
  
plots_classes

# eliminate CL under 15 min
plots_classes <-  over_rep02 %>% 
  filter(Class == "CL") %>%
  filter(!meanRT < 15 ) %>% 
  ggplot(aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16)) +
  labs(title = "CL") 

plots_classes

# eliminate PC with NA adducts
plots_classes <-  over_rep02 %>% 
  filter(Class == "PC") %>%
  filter(!Adduct =="Na-Gain" ) %>% 
  ggplot(aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16)) +
  labs(title = "PC") 

plots_classes

# now plotting all lipids with some of the modifications above
#filter(grepl( "CL", Class)) %>% 

all_lipids01 <- over_rep02 %>% 
  filter(!(Class =="CL" & meanRT < 15)) %>% 
  filter(!(Class=="PC" & Adduct =="Na-Gain"))  %>%
  ggplot(aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16)) +
  labs(title = "All lipids tider version")


all_lipids
#--------------------------------------------------------------------
# this function below is wrong why?
# '''
# No need to add an argument to the function with the object where you are going to store 
# the ggplot output. 
# 
# Before you had function(df, y){
# 
# ... 
# y <- ggplot()
# }
# 
# that is not need it, you can always pass the results of ggplot into an object withouth defining it previously. 
# That object would live withing the function. You can do this when you want to 
# keep working on it and will create an output later.
# 
# or when you want to store the result of the function into another object with return
# 
# '''
# this function below doesn't work
plot_class <-  function(df, class) {
    df <-  df %>% 
     filter(Class %in% c(class) )
    y <- ggplot(df, aes(meanRT, ExtractedMass)) +
    geom_point(aes(colour = batch, shape = Polarity)) +
    scale_shape_manual(values = c(6,16))
    
    return(y)
}

plot_class(over_rep02, "CL")


#same function below without object This works
plot_class <-  function(df) {
      df <-  df %>% 
         filter(Class == "CL" )
       y <- ggplot(df, aes(meanRT, ExtractedMass)) +
         geom_point(aes(colour = batch, shape = Polarity)) +
         scale_shape_manual(values = c(6,16)) + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) +
         labs(title = "CL")
       
          # y.classic <- y + theme_classic()
          # print(y)
          # print(y.classic)
     }


plot_class(over_rep02)

# -----------------------------------------------------
# this function is only taking a couple of classses, not all of them 
plot_class <-  function(data, var =c("x")) {
  
  data2 <- data %>% 
    filter(Class %in% var)
  
  ggplot(data2, aes(meanRT, ExtractedMass)) +
    geom_point(aes(colour)) +
    geom_point(aes(colour = batch, shape = Polarity)) +
    scale_shape_manual(values = c(6,16))
}


plot_class(over_rep02, var =c("PI","PG"))


extract_TC(molecule)

## a more dplyr like function  - but can't get it to work
extract_species <- function(dt, molecule){
  all_cols <- names(dt)[1:ncol(dt)]
  df <- dt %>% 
    ungroup() %>% 
    select(molecule) %>% 
    unique() %>% 
    group_by(1)
  df <- df %>% 
    mutate(class = str_extract(1, "[a-zA-Z]+"),
      # df <- as.data.frame(str_extract_all(molecule, "\\d+[:]\\d+")) %>%
      #   separate(1, c("C", "DB"), ":", convert = T) %>%
      #   summarise(
      #     Carbons = sum(C),
      #     DB = sum(DB)
      #   )
      C = sum(as.numeric(unlist(str_extract_all(1, "\\d+(?=\\:)")))),
      DB = sum(as.numeric((unlist(str_extract_all(1, "(?<=\\:)\\d+"))))),
      modification <- str_extract(str_extract(1, "(?<=\\().+?(?=\\))"), "[a-z]"),
      # species <- ifelse(modification == "t"|modification == "d", paste(modification, df$Carbons, ":", df$DB, sep = ""), 
      #                   ifelse(modification == "e"| modification == "p", paste(df$Carbons, modification, ":", df$DB, sep = ""), 
      #                   paste(df$Carbons, df$DB, sep = ":")))
      
      species = case_when(modification == "t"|modification == "d" ~ paste(modification, C, ":", DB, sep = ""), 
                           modification == "e"| modification == "p" ~  paste(C, modification, ":", DB, sep = ""), 
                           TRUE ~ paste(C, DB, sep = ":")))
  return(df)
}

extract_species(over_rep, Name)
molecule <- filter(over_rep, grepl("TG", Name)) 
molecule <- molecule$Name
                      
 C <- sum(as.numeric(unlist(str_extract("TG(18:1_18:2_16:0", "\\d+(?=\\:)"))))
 
