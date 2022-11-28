t <- read_csv ("./new/tidyPE462_abundances_calc_v5.csv", show_col_types = FALSE)
pdb <- read_lines("PlasticDB_Prokaryotic_genera.txt")
hcb <- read_lines("../Hydrocarbon_degraders_sorted_22_08.txt")
t <- t %>%
  filter(timepoint.days.%nin% c("mock_DNA","Neg_PCR")) %>% 
  filter(station%nin% c("mock_DNA","Neg_PCR"))%>% 
  mutate(station = str_replace_all(station,c("C13" = "Coastal station"))) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station"))) 
  
