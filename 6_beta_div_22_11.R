library(tidyverse)
library(Hmisc)

all <-read_csv("./new/tidyPE462_abundances_calc_v5.csv")%>%
  filter(timepoint.days.%nin% c("mock_DNA","Neg_PCR")) %>%
  filter(station %nin% c("mock_DNA","Neg_PCR")) %>%
  mutate(station = str_replace_all(station,c("C13" = "Coastal station"))) %>%
  mutate(station = str_replace_all(station,c("C05"= "open water station"))) %>%
  mutate_if(is.character,str_replace_all, "Glass_fiber", "SW") %>%
  mutate_if(is.character, as.factor)

coast <- all %>% filter(station %in% c("Coastal station")) 


# coast <- all %>% filter(station %in% c("Coastal station"))
# 
# pe0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
# peuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
# pet0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
# petuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
# ps0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
# psuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
# nylon0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
# nylonuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
#   mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
# wood0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
#   mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="wood")))
# coast2 <- coast %>% filter(polymer %nin% c("Glass_fiber")) %>%
#   bind_rows (pe0, peuv0, pet0, petuv0, ps0, psuv0, nylon0, nylonuv0, wood0)
