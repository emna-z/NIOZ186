library(tidyverse)
library(Hmisc)
library("RColorBrewer")
library(mapsf)
library(ggtext)
library(glue)
library(ggpubr)

no_white <- mf_get_pal(n = c( 1,30), pal = c("Grays","Mako" ), rev =c(F,F))
#PE_analysis==============
pe <- read_csv("./new/PE462_PE_v5.csv")
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
  mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus))))
#coast_top5_n <- coast_top5_n %>% 


pe_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus))))

pe_os_top5<- ggplot(data = os_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
  ggtitle(label = "Top genera abundance for PE in offshore station") 
#############################################stations together
pe_top5 <- pe %>% mutate(across(c(timepoint.days.),factor)) %>% 
  filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus, coast_nuv_g5$Genus,coast_uv_g5$Genus,coast_n_g5$Genus)))) %>% 
  mutate(pol_photo_station = str_replace_all(pol_photo_station,c("C05"= "open water", "C13" = "Coastal", "PE_"=""))) %>% 
  distinct() %>% 
  mutate(across(c(pol_photo_station),factor)) %>% 
  mutate(pol_photo_station = (fct_relevel(pol_photo_station,"no_UV_Coastal","UV_Coastal","wood_Coastal", "no_UV_open water", "UV_open water", "wood_open water")))



pe_top5<- ggplot(data = pe_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.75))+
  facet_grid(~ pol_photo_station, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() + theme(legend.position = "bottom" )+
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
  ggtitle(label = "Top genera abundance for PE polymer") 

plot(pe_top5)


########PET_analysis#############
pet <- read_csv("./new/PE462_PET_v5.csv")
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

coast_top5 <- coast %>% 
              select(timepoint.days.,detail,polymer_photo, Genus, Genus_rep_rel_abund, st_dev_Genus_abund) %>% 
              mutate(across(c(timepoint.days.),factor)) %>% 
              distinct()%>%
              filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus)))) %>%
              mutate(across(c(Genus), factor))

#coast_top5_n <- coast_top5_n %>% 
pet_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance", colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus)))) %>%
  mutate(across(c(Genus), factor))

pet_os_top5<- ggplot(data = os_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
  ggtitle(label = "Top genera abundance for PET in offshore station") 

############PS_analysis#############
ps <- read_csv("./new/PE462_PS_v5.csv")
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus)))) %>%
  mutate(across(c(Genus), factor))
#coast_top5_n <- coast_top5_n %>% 


ps_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus)))) %>%
  mutate(across(c(Genus), factor))

ps_os_top5<- ggplot(data = os_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
  ggtitle(label = "Top genera abundance for PS in offshore station") 

############Nylon_analysis#############
Nylon <- read_csv("./new/PE462_Nylon_v5.csv")
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(coast_nuv_g5$Genus,coast_uv_g5$Genus, coast_n_g5$Genus)))) %>%
  mutate(across(c(Genus), factor))
#coast_top5_n <- coast_top5_n %>% 


Nylon_coast_top5<- ggplot(data = coast_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
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
   mutate(across(c(timepoint.days.),factor)) %>% 
  distinct()%>% filter(Genus %in% (unique(c(os_nuv_g5$Genus,os_uv_g5$Genus, os_n_g5$Genus)))) %>%
  mutate(across(c(Genus), factor))

Nylon_os_top5<- ggplot(data = os_top5, mapping = aes(y = fct_relevel(Genus,rev), x = timepoint.days., fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance",colours = no_white, limits=c(0,0.7))+
  facet_grid(~ polymer_photo, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"), strip.text.x = element_text(face = "bold")) +
  ggtitle(label = "Top genera abundance for Nylon in offshore station") 



plot(ggarrange(pet_coast_top5, pet_os_top5, ncol = 2, common.legend = T))
plot(ggarrange(ps_coast_top5, ps_os_top5, Nylon_coast_top5, Nylon_os_top5, ncol = 2, nrow = 2, common.legend = T))

