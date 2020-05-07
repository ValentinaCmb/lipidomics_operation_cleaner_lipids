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

# prepare  the data to make a library for CD that contains only Name and Formula
batchall <- batch_bind %>% 
  select(-"X11", -Window, -batch) %>% 
  unique() %>% 
  arrange(CompoundName, RT) %>% 
  rename(Name = CompoundName)

CD_library <- batchall %>% 
  rename(Name = CompoundName) %>% 
  rename( Formula = ChemicalFormula) %>% 
  select(Name, Formula) %>% 
  unique()

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


