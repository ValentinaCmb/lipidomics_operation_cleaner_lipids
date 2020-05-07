library(tidyverse)
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

batch_bind <- batch_01 %>% 
  rbind(batch_02, batch_03, batch_04) %>% 
  mutate(RT = round(RT, 1)) %>% 
  unique()

batchall <- batch_bind %>% 
  select(-"X11", -Window) %>% 
  unique() %>% 
  arrange(CompoundName, RT) 

tidier_batch <- batchall %>% 
  rename(Name = CompoundName) %>% 
  rename( Formula = ChemicalFormula) %>% 
  #select(Name, Formula) %>% 
  unique() %>% 
  mutate(Class = str_extract(`Name`, "[a-zA-Z]+"))  
  # mutate(TC = word(`Name`, 1, 1, sep = "_"),
  #        DB = word(`Name`, -1, -1, sep = "_"))
  


over_rep <- tidier_batch %>% 
  group_by(Class, Formula, Ionization, Polarity, Adduct, ChargeState, ExtractedMass, batch) %>% 
  filter(RT < 25,
         ChargeState == 1) %>% 
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
  mutate(batch = word(`batch`, 2,2, sep ="_"))


# plot based on batches
ggplot(over_rep02, aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16))

# plot based on Classes


class  <- over_rep02 %>% 
  filter(Class == "CL")

ggplot(class, aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16))


# this function below is wrong why?
plot_class <-  function(x) {
  over_rep02 %>% 
    filter(Class == "x" ) %>% 
  ggplot(class, aes(meanRT, ExtractedMass)) +
    geom_point(aes(colour = batch, shape = Polarity)) +
    scale_shape_manual(values = c(6,16))
}

plot_class(PE)


