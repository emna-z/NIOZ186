library(tidyverse)
library(Hmisc)
library(ggpubr)

alpha_tab <- read_csv("./new/alpha_div_indexes_microbiome_package.csv")
alpha_tab <- alpha_tab %>%  filter(station %in% c("C05", "C13")) %>% filter(material %nin% c("Glass_fiber")) %>% 
  filter(timepoint.days. %nin% c("15"))
alpha_tab$polymer <- if_else(alpha_tab$material=="wood", str_c(alpha_tab$material), alpha_tab$polymer)
alpha_tab_coast <- alpha_tab %>% filter(station %in% c("C13")) %>% mutate(across(c(treatment),factor)) %>% 
  mutate(across(c(treatment),factor))
alpha_tab_os <- alpha_tab %>% filter(station %in% c("C05")) %>% mutate(across(c(treatment),factor)) %>% 
  mutate(across(c(treatment),factor))

###create graphs

p_coast <-ggplot(alpha_tab_coast,aes(x = reorder(timepoint.days., as.numeric(timepoint.days.)), y = evenness_simpson, fill=treatment))+
  geom_boxplot(position=position_dodge(1))  +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_minimal()+
  facet_wrap(~polymer, nrow = 5)+
  labs(title="Alpha diversity Coastal station",x="incubation time (days)", y = "Simpson eveness index")+ ylim(0,0.35)
p_os <-ggplot(alpha_tab_os,aes(x = reorder(timepoint.days., as.numeric(timepoint.days.)), y = evenness_simpson, fill=treatment))+
  geom_boxplot(position=position_dodge(1)) +scale_y_continuous(limits = c(0,0.35))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_minimal()+
  facet_wrap(~polymer, nrow = 5)+
  labs(title="Alpha diversity open water station",x="incubation time (days)", y = "Simpson eveness index")
fig1 <- ggarrange(p_coast,p_os, common.legend = T) 
#ggexport(fig,filename = "simpson_index.pdf")
#fig1

sh_coast <-ggplot(alpha_tab_coast,aes(x = reorder(timepoint.days., as.numeric(timepoint.days.)), y = diversity_shannon, fill=treatment))+
  geom_boxplot(position=position_dodge(1))+scale_y_continuous(limits = c(0,6))  +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_minimal()+
  facet_wrap(~polymer, nrow = 5)+
  labs(title="Alpha diversity Coastal station",x="incubation time (days)", y = "Shannon diversity index")

sh_os <-ggplot(alpha_tab_os,aes(x = reorder(timepoint.days., as.numeric(timepoint.days.)), y = diversity_shannon, fill=treatment))+
  geom_boxplot(position=position_dodge(1))+scale_y_continuous(limits = c(0,6))  + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_minimal()+
  facet_wrap(~polymer, nrow = 5)+
  labs(title="Alpha diversity open water station",x="incubation time (days)", y = "Shannon diversity index")
fig2 <- ggarrange(sh_coast,sh_os, common.legend = T) 
#ggexport(fig,filename = "shannon_index.pdf")
#fig2

