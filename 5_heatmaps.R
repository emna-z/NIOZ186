k <- read_csv("./new/Kingdom_PE462_v5.csv")
p <- read_csv("./new/phyla_PE462_v5.csv")
o <-read_csv("./new/order_PE462_v5.csv")

ggplot(o, aes(x=Material, y= Kingdom_rep_rel_abund, fill=Kingdom))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic2()+  facet_grid (timepoint.days.~treatment)

g_all <- read_csv("./new/genus_PE462_v5.csv")
# coast <- g_all %>% filter(station %in% c("C13")) %>% 
#   mutate(station = str_replace_all(station,c("C13" = "Coastal station")))
# os <- g_all %>% filter(station %in% c("C05")) %>% 
#   mutate(station = str_replace_all(station,c("C05"= "open water station"))) 
g_all$timepoint.days.
unique(g_all$timepoint.days.)
unique(g_all$station)
g_all <- g_all %>% filter(timepoint.days.%nin% c("mock_DNA","Neg_PCR")) %>% 
  filter(station%nin% c("mock_DNA","Neg_PCR"))%>% 
  mutate(station = str_replace_all(station,c("C13" = "Coastal station"))) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station"))) 

top5 <- g_all%>% 
  mutate(across(c(timepoint.days., detail),factor))%>% distinct() %>% 
  group_by(detail) %>% slice_max(order_by = Genus_rep_rel_abund, n = 3)
top5 <- top5 %>% 
  filter(Genus %nin% c("unassigned"))

coast <- g_all %>% filter(station %in% c("Coastal station")) %>% 
  filter (Genus %in% (unique(top5$Genus)))

pe0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
peuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
pet0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
petuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
ps0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
psuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
nylon0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
nylonuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
wood0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="wood")))
coast2 <- coast %>% filter(polymer %nin% c("Glass_fiber")) %>%
  bind_rows (pe0, peuv0, pet0, petuv0, ps0, psuv0, nylon0, nylonuv0, wood0)

x <- seq(0, 1, length.out = 200)
library("ggh4x")

c <- ggplot(data = coast2, mapping = aes(y = fct_relevel(Genus,rev), 
            x = fct_relevel(timepoint.days.,  c("0","5", "10","30","45")), fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance", colours = no_white, limits=c(0,0.75), breaks = x )+
  facet_nested(~ polymer + treatment, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic"),
        axis.text.x = element_text(face = "bold"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#585858"), strip.text.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position="none") +
  ggtitle(label = "Top genera relative abundance in coastal station") 
c

library(rvg)
p_o_v <- dml(ggobj=c)
ppt<- read_pptx(path = "./plots/heat11.pptx") %>% add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with( value = p_o_v, location = ph_location_fullsize() )
print(ppt, target = "./plots/heat11.pptx")



###
os <- g_all %>% filter(station %in% c("open water station")) %>% 
  filter (Genus %in% (unique(top5$Genus)))

pe0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
peuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
pet0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
petuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
ps0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
psuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
nylon0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
nylonuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
wood0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="wood")))
os2 <- os %>% filter(polymer %nin% c("Glass_fiber")) %>%
  bind_rows (pe0, peuv0, pet0, petuv0, ps0, psuv0, nylon0, nylonuv0, wood0)

x <- seq(0, 0.75, length.out = 200)
library("ggh4x")

osp <- ggplot(data = os2, mapping = aes(y = fct_relevel(Genus,rev), 
                                         x = fct_relevel(timepoint.days.,  c("0","5", "10","30","45")), fill = Genus_rep_rel_abund)) + 
  geom_tile()  +scale_fill_gradientn(name = "relative abundance", colours = no_white, limits=c(0,0.75), breaks = x )+
  facet_nested(~ polymer + treatment, scales = "free_x", space = "free_x") + xlab(label = "incubation time (days)") +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic"),
        axis.text.x = element_text(face = "bold"),
        strip.background = element_rect(fill = "#FFFFFF", color = "#585858"), strip.text.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position="none") +
  ggtitle(label = "Top genera relative abundance in open water station") 
osp
