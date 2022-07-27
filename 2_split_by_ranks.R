######libraries#####
library(tidyverse)
library(Hmisc)
library("RColorBrewer")

############tidy data import
#import of final file created in script 1_PE462_data_prep.R
t3 <- read_csv("new/tidyPE462_abundances_calc_v5.csv")

#####split by ranks###########
Kingdom <- t3  %>%  select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,polymer_photo,
                           material, detail,Kingdom, Kingdom_rep_rel_abund,Kingdom_st_dev_abund_samples)%>% 
  distinct() 
#write_csv(Kingdom, "./new/Kingdom_PE462_v5.csv")

Phyla <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail, Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
                       distinct() 
#write_csv(Phyla, "./new/phyla_PE462_v5.csv")

Class <-t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail,
                      Class, Class_rep_rel_abund, st_dev_Class_abund )%>% 
distinct()  
#write_csv(Class, "./new/class_PE462_v5.csv")

Order <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail,
                       Order, Order_rep_rel_abund, st_dev_Order_abund )%>% 
  distinct() 
#write_csv(Order, "./new/order_PE462_v5.csv")

Family <- t3 %>% select(station,timepoint.days.,treatment,polymer, polymer_photo, material, detail,
                        Family, Family_rep_rel_abund, st_dev_Family_abund )%>% 
  distinct()
#write_csv(Family, "./new/family_PE462_v5.csv")

Genus <- t3%>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
#write_csv(Genus, "./new/genus_PE462_v5.csv")

Species <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                         polymer_photo, material, detail,  Species, Species_rep_rel_abund, st_dev_Species_abund)%>% 
  distinct()
#write_csv(Species, "./new/species_PE462_v5.csv") 
