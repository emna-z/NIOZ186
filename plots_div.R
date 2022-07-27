library(scales)
library(RColorBrewer)
library(ggthemes)
library(Polychrome)



k <- read_csv("./new/Kingdom_PE462_v5.csv")
k <- k %>% filter(station %in% c("C05","C13")) %>%
     filter(timepoint.days.%nin% c("0")) %>% 
     mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "Coastal station"))) %>% 
     mutate(across(c(polymer,polymer_photo,timepoint.days.,station),factor)) %>% 
     distinct() 
k

king <- ggplot(k, aes(x=polymer_photo, y= Kingdom_rep_rel_abund, fill=Kingdom))+
     geom_bar(stat="identity", position="stack")+ scale_fill_manual(values=c("#FF0000" ,"#00A08A" ,"#F2AD00" ))+
     scale_y_continuous(labels=percent)+
     theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     xlab("")+ylab("Relative Abundance")+
     facet_grid (fct_relevel(timepoint.days., c('5',"10", "30", "45"))~station)+ 
     theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=10, face="bold", angle = 0))
king
#ggexport(king,filename = "./plots/kingdom_rel_abund_2.pdf")

 
 
phyla <- read_csv("./new/phyla_PE462_v5.csv")
p <- phyla %>% filter(station %in% c("C05","C13")) %>%
  filter(timepoint.days.%nin% c("0")) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "Coastal station"))) %>% 
  mutate(across(c(polymer,polymer_photo,timepoint.days.,station),factor)) %>% 
  mutate(Phylum, Phylum = if_else(Phyla_rep_rel_abund < 0.05, str_c("others <5%"), Phylum)) %>% 
  mutate(Phylum = fct_relevel(Phylum,"others <5%", after = Inf))

ph <- ggplot(p, aes(x=polymer_photo, y= Phyla_rep_rel_abund, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+ scale_fill_ptol()+
  scale_y_continuous(labels=percent)+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint.days., c('5',"10", "30", "45"))~station)+ 
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=10, face="bold", angle = 0))
ph

o <- read_csv("./new/order_PE462_v5.csv")
o <- o %>% filter(station %in% c("C05","C13")) %>%
  filter(timepoint.days.%nin% c("0")) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "Coastal station"))) %>% 
  mutate(across(c(polymer,polymer_photo,timepoint.days.,station),factor))%>%  
  mutate(Order, Order = if_else(Order_rep_rel_abund < 0.05, str_c("others <5%"), Order)) %>% 
  mutate(Phylum = fct_relevel(Order,"others <5%", after = Inf))


sky <- Polychrome::sky.colors()
names(sky) <- NULL
polychrome_light <- light.colors()
swatch(sky)
names(polychrome_light) <- NULL
show_col(polychrome_light)
swatch(polychrome_light)

p_o <- ggplot(o, aes(x=polymer_photo, y= Order_rep_rel_abund, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = sky)+
  scale_y_continuous(labels=percent)+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint.days., c('5',"10", "30", "45"))~station)+ 
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=10, face="bold", angle = 0))
p_o


t <- read_csv("./new/tidyPE462_v5.csv")
source("./rel_abund_by_ranks.R")
t1 <-   rel_abund_by_ranks(t, Sample, detail)
