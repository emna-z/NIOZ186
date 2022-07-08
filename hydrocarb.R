library(tidyverse)
library(Hmisc)
hydrocarb <- read_lines("../Hydrocarbon degraders.txt")
Genus <- read_csv("./new/genus_PE462.csv")
hd <- Genus %>% filter (Genus %in% hydrocarb) %>% filter(station %in% c("C05","C13")) %>% filter(timepoint.days. %nin% c("0","15")) %>% 
  filter(Genus_rep_rel_abund > 0.01)
gen_keep <- unique(hd$Genus)
hd_coast <- Genus %>% filter (Genus %in% gen_keep) %>% filter(station %in% c("C13")) %>% filter(timepoint.days. %nin% c("0","15"))
hd_os <- Genus %>% filter (Genus %in% gen_keep) %>% filter(station %in% c("C05")) %>% filter(timepoint.days. %nin% c("0","15"))  
ggplot(hd, aes(x=polymer_photo, y= Genus_rep_rel_abund, fill=Genus))+
  geom_point()+theme_classic()+   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="right")+ 
  facet_grid (station~timepoint.days.)

heat_coast <- ggplot(hd_coast, aes(x = reorder(timepoint.days., as.numeric(timepoint.days.)), y = rev(Genus),
                             fill = Genus_rep_rel_abund)) + geom_tile() +  
              scale_fill_gradientn(colours = c(rep("#ecefb7", 2), "#7fc8b9", "#3a81b5"),
                                   na.value = "#FFFFE0",
                                   values = c(0, 0.05, 0.0501, 1),
                                   limits = c(0, 1)) +
              theme(
                axis.text.x=element_text(hjust=1,vjust=0),
                panel.background = element_blank()) +
  ylab("hydrocarbon degraders Genera") +
  xlab("")+facet_grid (~ polymer_photo)
heat_coast


bubble_plot <- ggplot(hd_coast,aes(x=reorder(timepoint.days., as.numeric(timepoint.days.)),y=reorder(Genus, as.character(Genus), decreasing = T))) +  
geom_point(aes(size=Genus_rep_rel_abund, colour = factor(treatment))) +  
  scale_fill_manual(values = c("green", "red", "yellow"))+
  theme_classic()+theme(
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  ylab("hydrocarbon degraders Genera") +
  xlab("") +  facet_grid (~ polymer_photo)
bubble_plot 
  # geom_bar(stat="identity", position="stack")+ scale_color_brewer()+
  #theme_classic()+   theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="right")+ facet_grid (timepoint.days.~station)
