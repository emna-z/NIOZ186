t3 <- read_csv("./new/tidyPE462_abundances_calc_w.csv")
t3$polymer <- if_else(t3$material=="wood", str_c(t3$material), t3$polymer)
Kingdom <- t3  %>%  select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,polymer_photo,
                           material, detail,Kingdom, Kingdom_rep_rel_abund,Kingdom_st_dev_abund_samples)%>% 
  distinct() 
#write_csv(Phyla, "./new/Kingdom_PE462.csv")

Phyla <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail, Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
  distinct() 
head(Phyla)
Phyla <- Phyla %>% filter(timepoint.days.%in% c("10","45")) %>% filter(polymer %nin% c("Neg_PCR","mock_DNA" ))
#write_csv(Phyla, "./new/phyla_PE462.csv")

Class <-t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail,
                      Class, Class_rep_rel_abund, st_dev_Class_abund )%>% 
distinct()  
#write_csv(Class, "./new/class_PE462.csv")
Order <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail,
                       Order, Order_rep_rel_abund, st_dev_Order_abund )%>% 
  distinct() 
#write_csv(Order, "./new/order_PE462.csv")
Family <- t3 %>% select(station,timepoint.days.,treatment,polymer, polymer_photo, material, detail,
                        Family, Family_rep_rel_abund, st_dev_Family_abund )%>% 
  distinct()
#write_csv(Family, "./new/family_PE462.csv")
Genus <- t3%>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
#write_csv(Genus, "./new/genus_PE462.csv")

Species <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                         polymer_photo, material, detail,  Species, Species_rep_rel_abund, st_dev_Species_abund)%>% 
  distinct()
#write_csv(Species, "./new/species_PE462.csv") 
