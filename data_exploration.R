library(tidyverse)
############PE_analysis#############
pe <- read_csv("./PE462_PE.csv")
unique (pe$Family) ==
######subset per station#########
coast <- pe %>% filter(station %in% c("C05"))
os <- pe %>% filter(station %in% c("C13"))
#####subset stations per treatment######
coast_uv <- coast %>% filter(treatment %in% c("UV"))
coast_nuv <- coast %>% filter(treatment %in% c("no_UV"))
coast_n <- coast %>% filter(treatment %in% c("no"))
##same for 2nd station
os_uv <- os %>% filter(treatment %in% c("UV"))
os_nuv <- os %>% filter(treatment %in% c("no_UV"))
os_n <- os %>% filter(treatment %in% c("no"))
#length(unique(coast_uv$Genus_rep_rel_abund))
#length(unique(coast_uv$Description))
#length(unique(coast_uv$Genus)) 

############top5 taxa in subsets#############
coast_uv_g5 <- coast_uv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
              mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
              group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
coast_nuv_g5 <- coast_nuv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
coast_n_g5 <- coast_n %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)

unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus) )

####back to stations datasets to plot only top5 taxa for each condition
coast_top5 <- coast %>% select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PE_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus))))
#coast_top5_n <- coast_top5_n %>% 


g<- ggplot(data = coast_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
                        geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PE in coastal station") 
g
superheat(coast_top5)
