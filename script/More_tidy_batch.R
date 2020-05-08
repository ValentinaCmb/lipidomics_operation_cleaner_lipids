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
    filter((!Class =="CL" & meanRT > 15)|(Class =="PC" & Adduct !="Na-Gain"))  %>%
    ggplot(aes(meanRT, ExtractedMass)) +
  geom_point(aes(colour = batch, shape = Polarity)) +
  scale_shape_manual(values = c(6,16)) +
  labs(title = "All lipids tid version") 

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
    geom_point(aes(colour = batch, shape = Polarity)) +
    scale_shape_manual(values = c(6,16))
}


plot_class(over_rep02, var =c("PI","PG"))



