library(tidyverse)
library(Hmisc)

#PE_analysis==============
pe <- read_csv("./new/PE462_PE.csv")
pe <- pe %>% filter(Genus %nin% c("unassigned", "Unassigned"))
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

    ############top5 taxa in coast#############
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


pe_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PE in coastal station") 
pe_coast_top5

    ############top5 taxa in offshore#############
os_uv_g5 <- os_uv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_nuv_g5 <- os_nuv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_n_g5 <- os_n %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)

unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus) )

####back to stations datasets to plot only top5 taxa for each condition
os_top5 <- os %>% select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PE_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus))))

pe_os_top5<- ggplot(data = os_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PE in offshore station") 



########PET_analysis#############
pet <- read_csv("./new/PE462_PET.csv")
pet <- pet %>% filter(Genus %nin% c("unassigned", "Unassigned"))
    #subset petr station#########
coast <- pet %>% filter(station %in% c("C05"))
os <- pet %>% filter(station %in% c("C13"))
    #subset stations petr treatment######
coast_uv <- coast %>% filter(treatment %in% c("UV"))
coast_nuv <- coast %>% filter(treatment %in% c("no_UV"))
coast_n <- coast %>% filter(treatment %in% c("no"))
##same for 2nd station
os_uv <- os %>% filter(treatment %in% c("UV"))
os_nuv <- os %>% filter(treatment %in% c("no_UV"))
os_n <- os %>% filter(treatment %in% c("no"))

    #top5 taxa in coast#############
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
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PET_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus))))
#coast_top5_n <- coast_top5_n %>% 


pet_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PET in coastal station") 
pet_coast_top5

    #top5 taxa in offshore#############
os_uv_g5 <- os_uv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_nuv_g5 <- os_nuv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_n_g5 <- os_n %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)

unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus) )

####back to stations datasets to plot only top5 taxa for each condition
os_top5 <- os %>% select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PET_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus))))

pet_os_top5<- ggplot(data = os_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PET in offshore station") 

############PS_analysis#############
ps <- read_csv("./new/PE462_PS.csv")
ps <- ps %>% filter(Genus %nin% c("unassigned", "Unassigned"))
#subset per station#########
coast <- ps %>% filter(station %in% c("C05"))
os <- ps %>% filter(station %in% c("C13"))
#subset stations petr treatment######
coast_uv <- coast %>% filter(treatment %in% c("UV"))
coast_nuv <- coast %>% filter(treatment %in% c("no_UV"))
coast_n <- coast %>% filter(treatment %in% c("no"))
##same for 2nd station
os_uv <- os %>% filter(treatment %in% c("UV"))
os_nuv <- os %>% filter(treatment %in% c("no_UV"))
os_n <- os %>% filter(treatment %in% c("no"))

#top5 taxa in coast#############
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
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PS_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus))))
#coast_top5_n <- coast_top5_n %>% 


ps_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PS in coastal station") 
ps_coast_top5

#top5 taxa in offshore#############
os_uv_g5 <- os_uv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_nuv_g5 <- os_nuv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_n_g5 <- os_n %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)

unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus) )

####back to stations datasets to plot only top5 taxa for each condition
os_top5 <- os %>% select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "PS_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus))))

ps_os_top5<- ggplot(data = os_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for PS in offshore station") 

############Nylon_analysis#############
Nylon <- read_csv("./new/PE462_Nylon.csv")
Nylon <- Nylon %>% filter(Genus %nin% c("unassigned", "Unassigned"))
  #subset per station#########
  coast <- Nylon %>% filter(station %in% c("C05"))
  os <- Nylon %>% filter(station %in% c("C13"))
  #subset stations petr treatment######
  coast_uv <- coast %>% filter(treatment %in% c("UV"))
  coast_nuv <- coast %>% filter(treatment %in% c("no_UV"))
  coast_n <- coast %>% filter(treatment %in% c("no"))
  ##same for 2nd station
  os_uv <- os %>% filter(treatment %in% c("UV"))
  os_nuv <- os %>% filter(treatment %in% c("no_UV"))
  os_n <- os %>% filter(treatment %in% c("no"))

  #top5 taxa in coast#############
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
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "Nylon_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus))))
#coast_top5_n <- coast_top5_n %>% 


Nylon_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for Nylon in coastal station") 
Nylon_coast_top5

  #top5 taxa in offshore#############
os_uv_g5 <- os_uv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_nuv_g5 <- os_nuv %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)
os_n_g5 <- os_n %>% select(timepoint.days.,detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate(across(c(timepoint.days.),factor))%>% distinct() %>% 
  group_by(timepoint.days.) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5)

unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus) )

####back to stations datasets to plot only top5 taxa for each condition
os_top5 <- os %>% select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
  mutate( polymer_photo = replace(polymer_photo, polymer_photo == "Nylon_wood_no", "ctrl_wood")) %>% mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus))))

Nylon_os_top5<- ggplot(data = os_top5, mapping = aes(y = Genus, x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradient(name = "relative abundance",low = "#FFFFFF", high = "#483D8B")+
  facet_grid(~ polymer_photo, switch = "x", scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = "Top genera abundance for Nylon in offshore station") 


# NMDS --------------------------------------------------------------------
gen <- read_csv("./new/genus_PE462.csv")
length(unique(gen$Genus))
length(unique(gen$detail))
gen %>% filter(detail %nin% c("mock_DNA","Glass_fiber_no_C05_0", "Neg_PCR", "Glass_fiber_no_C13_0" )) %>% 
  select(detail, Genus, Genus_rep_rel_abund) %>% mutate(across(c(Genus),factor)) %>% 
  pivot_wider(names_from = detail, values_from = Genus_rep_rel_abund, values_fill = 0)
