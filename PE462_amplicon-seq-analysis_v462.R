##%######################################################%##
#                                                          #
####              Exploratory Analysis of               ####
####        16S/18S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####               Project: PE462 films                 ####
#                                                          #
##%######################################################%##

################Packages_init###################
library(tidyverse)
library(phyloseq)
library(grid)
library(scales)
library(vegan)
library(rmarkdown)
library(knitr)
library(DESeq2)
library(microbiome)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(usedist)
library(heatmaply)
library(Hmisc)
library(kableExtra)

#####data_import######
tax1 <- as.matrix(read.delim("./original/representative_seq_set_tax_assignments.txt", row.names = 1, na.strings = "NA"))
tax2 <- read.delim("./original/representative_seq_set_tax_assignments_unpaired.txt", row.names = 1, na.strings = "NA")
rownames(tax2) <- paste0("asv.",(length(rownames(tax1))+1):(length(rownames(tax1))+length(rownames(tax2))))
tax2 ['Species']= NA
tax2 <- as.matrix(tax2)
tax <- rbind(tax1,tax2) %>% tax_table()
otu1 <- as.matrix(read.delim("./original/asv_table.txt", row.names = 1))
otu2 <- read.delim("./original/asv_table_unpaired.txt", row.names = 1)
rownames(otu2) <- rownames(tax2)
samples_missing <- colnames(otu1)[c(which(colnames(otu1) %nin% colnames(otu2)))]
for (i in 1:length(samples_missing)) {otu2[samples_missing[i]] <- 0 }
otu2 <- as.matrix (otu2)
otu <-rbind(otu1,otu2) %>%  otu_table(taxa_are_rows = T)
map <- sample_data(read.delim("./original/mapPE462.txt", row.names = 1, na.strings = c("NA", "")))
physeq_object = merge_phyloseq(otu, tax, map)                 


####basic_info##############
summarize_phyloseq(physeq_object)
ntaxa(physeq_object)
nsamples(physeq_object)
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) 
physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE) #no singletons
min(taxa_sums(physeq_object))
max(sample_sums(physeq_object))
##########getting rid of wonky taxonomy assignments ###########
get_taxa_unique(physeq_object, "Kingdom") 
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c("", "Unassigned"))
get_taxa_unique(physeq_object, "Kingdom")
get_taxa_unique(physeq_object, "Phylum") 

###rid of NA in taxonomy
taxo <- as.data.frame(physeq_object@tax_table)


for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (is.na(taxo[i,y]) || any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome", "Metagenome","unknown", "Unknown","NA")))) {
      taxo[i,y] <- "unassigned" }
  }
}

taxo <- tax_table(as.matrix(taxo))


physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map)

###rid of chloroplast & Mitochodria
any((get_taxa_unique(physeq_object, "Order") == "Chloroplast"))
any((get_taxa_unique(physeq_object, "Family") == "Mitochondria"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c("Chloroplast")) 
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria"))

physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE)#no singletons

############### alpha div ###################

#microbiome package
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
alpha_tab <- alpha_tab %>% cbind(data.frame(physeq_object@sam_data))
#write_csv(alpha_tab, file = "./new/alpha_div_indexes_microbiome_package.csv")
alpha_tab <- read_csv("./new/alpha_div_indexes_microbiome_package.csv")
#phyloseq package
a_div <- phyloseq::estimate_richness(physeq_object)
#write_csv(a_div, file = "./new/alpha_div_phyloseq.csv")

#p <- ggboxplot(metad, x = "material", y = "evenness_simpson",
#               color = "material", palette =c("#5FB233FF" ,"#6A7F93FF" ,"#F57206FF" ,"#EB0F13FF", "#8F2F8BFF", "#1396DBFF"),
#              add = "jitter", shape = "treatment", size = 1) + facet_wrap(~timepoint.days.)
#p

alpha <- read_csv("./new/alpha_div_indexes_microbiome_package.csv")

plot <- plot_richness(physeq_object, "material", "treatment", measures="Chao1")+facet_grid(treatment~timepoint.days.)
plot + geom_boxplot(data=plot$data, aes(material,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Chao1", x =NULL , y = NULL )+theme_light()

plot <- plot_richness(physeq_object, "polymer", "treatment", measures="Simpson")+facet_grid(treatment~station)
plot + geom_boxplot(data=plot$data, aes(polymer,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Simpson", x =NULL , y = NULL )+theme_light()

microbiome::plot_taxa_prevalence(physeq_object, "Kingdom")+ theme(legend.position = "none") #prevalence


#########merge samples per surface all replicates together##########
source("./tidy_psmelt.R")
tidy_physeq_asv <- tidy_psmelt(physeq_object)
tidy_physeq_asv$Species <- if_else(!tidy_physeq_asv$Species=="unassigned", str_c(tidy_physeq_asv$Genus," ",tidy_physeq_asv$Species), tidy_physeq_asv$Species)
any(tidy_physeq_asv$Species == correct_species)
tidy_physeq_asv$detail <- if_else(tidy_physeq_asv$material=="wood", str_c(tidy_physeq_asv$material,"_",tidy_physeq_asv$station,"_",tidy_physeq_asv$timepoint.days.), tidy_physeq_asv$detail)
any(tidy_physeq_asv$detail == wood_detail)

#write_csv(tidy_physeq_asv, "./new/tidyPE462.csv")
tidy_physeq_asv <- read.csv("./new/tidyPE462.csv")

t3 <- tidy_physeq_asv  %>% group_by(Sample) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup() %>%
  group_by(detail) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
  ungroup() %>% 
  #Kingdom_section
  group_by(Sample, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  group_by(detail, Kingdom) %>% 
  mutate(Kingdom_st_dev_abund_samples = sd(Kingdom_rel_abund_Sample)) %>% # standard dev of Kingdom relative abundances between replicates of detail (ployner_timepoint.days._treatment)
  mutate(Kingdom_rep_rel_abund = sum(rep_rel_abund)) %>% #Kingdom relative abundance per samples of desc 
  ungroup() %>%
  #Phylum_section
  group_by(Sample, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Phylum) %>% 
  mutate(st_dev_Phylum_abund = sd(Phylum_rel_abund_Sample)) %>%
  mutate(Phyla_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Class_section
  group_by(Sample, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Class) %>% 
  mutate(st_dev_Class_abund = sd(Class_rel_abund_Sample)) %>%
  mutate(Class_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Order_section
  group_by(Sample, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Order) %>% 
  mutate(st_dev_Order_abund = sd(Order_rel_abund_Sample)) %>%
  mutate(Order_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Family_section
  group_by(Sample, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Family) %>% 
  mutate(st_dev_Family_abund = sd(Family_rel_abund_Sample)) %>%
  mutate(Family_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Genus_section
  group_by(Sample, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Genus) %>% 
  mutate(st_dev_Genus_abund = sd(Genus_rel_abund_Sample)) %>%
  mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Species_section
  group_by(Sample, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Species) %>% 
  mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup()  

polymer_station <-  str_c(t3$polymer, "_", t3$station)
polymer_photo <- str_c(t3$polymer, "_", t3$treatment)
pol_photo_station <- str_c(t3$polymer_photo, "_", t3$station)

t3 <- t3 %>% 
  add_column(polymer_station, .before ="Kingdom") %>% 
  add_column(polymer_photo, .before ="Kingdom") %>% 
  add_column(pol_photo_station, .before ="Kingdom")

t3$polymer_station <- if_else(t3$material=="wood", str_c(t3$material,"_",t3$station), t3$polymer_station)
t3$polymer_photo <- if_else(t3$material=="wood", str_c(t3$material), t3$polymer_photo)
t3$pol_photo_station <- if_else(t3$material=="wood", str_c(t3$material,"_",t3$station), t3$pol_photo_station)

#write_csv(t3, "./new/tidyPE462_abundances_calc_w.csv")  

#t <- read_csv("./new/tidyPE462_abundances_calc_w.csv")

mock <- t3 %>% filter(polymer %in% c("mock_DNA"))
mock$Genus
D10_45 <- t3 %>% filter(timepoint.days.%in% c("10","45")) %>% filter(polymer %nin% c("Neg_PCR","mock_DNA" ))
#write_csv(D10_45, "./analysis/D10_45.csv")
##############tables for each rank#############
t3 <- read_csv("./new/tidyPE462_abundances_calc_w.csv")
t3$polymer <- if_else(t3$material=="wood", str_c(t3$material), t3$polymer)
Kingdom <- t3  %>%  select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,polymer_photo,
                           material, detail,Kingdom, Kingdom_rep_rel_abund,Kingdom_st_dev_abund_samples)%>% 
  distinct() 
#write_csv(Kingdom, "./new/Kingdom_PE462.csv")

#colors_vector_to_personalize#
CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 

ggplot(Kingdom, aes(x=polymer_photo, y= Kingdom_rep_rel_abund, fill=Kingdom))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic2()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid (timepoint.days.~station)


#creating table phyla
Phyla <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                polymer_photo, material, detail, Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
  distinct() 
head(Phyla)
Phyla <- Phyla %>% filter(timepoint.days.%in% c("10","45")) %>% filter(polymer %nin% c("Neg_PCR","mock_DNA" ))
#write_csv(Phyla, "./analysis/phyla_PE462.csv")

#Phyla <- read_csv("./analysis/phyla_PE462.csv")
Phyla_10 <- Phyla %>% filter(timepoint.days.=="10")
p10 <- ggplot(Phyla_10, aes(x=polymer_photo, y= Phyla_rep_rel_abund, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+
  theme_classic2()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")+
  facet_grid (~station)+ ggtitle ("day 10")

Phyla_45 <- Phyla %>% filter(timepoint.days.=="45")
p45 <- ggplot(Phyla, aes(x=polymer_photo, y= Phyla_rep_rel_abund, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+
  theme_classic2()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")+
  facet_grid (timepoint.days.~station)
p45


bact <- Phyla %>% filter(Phylum %in% c("Proteobacteria","Planctomycetota", "Myxococcota","Bacteroidetes","Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota",
                                       "Deinococcota", "Firmicutes", "Actinobacteria", "Cyanobacteria"))


p2 <- ggplot(bact, aes(x=detail, y= Phyla_rep_rel_abund, fill=Phylum)) +
  geom_col(stat="identity", position="dodge")+ scale_fill_manual(values = CPCOLS) +
  theme_pubclean()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  labs(title = "OTU", subtitle = NULL, x =NULL , y = NULL ) + ylim(0,1)+
  geom_errorbar(aes(ymin=Phyla_rep_rel_abund, ymax=Phyla_rep_rel_abund+st_dev_Phylum_abund), width=.2,
                position=position_dodge(.9))
p2

Class <-t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail,
                       Class, Class_rep_rel_abund, st_dev_Class_abund )%>% 
  distinct()  
head(Class)
#write_csv(Class, "./analysis/class_PE462.csv")
#Class <- read_csv("../../Analysis/class.csv")
Class <- Class %>% mutate(Class = ifelse(Class_rep_rel_abund<0.01, "others<0.01", Class))
ggplot(Class, aes(x=pol_photo_station, y= Class_rep_rel_abund, fill=Class))+
  geom_bar(stat="identity", position="stack")+
  theme_classic()+ theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid (timepoint.days.~station)

Order <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail,
                       Order, Order_rep_rel_abund, st_dev_Order_abund )%>% 
  distinct() #
head(Order)
#write_csv(Order, "./analysis/order_PE462.csv")
#Order <- read_csv("../../Analysis/order_v460.csv")
Order <- Order %>% mutate(Order = ifelse(Order_rep_rel_abund<0.01, "others<0.01", Order))

ggplot(Order, aes(x=pol_photo_station, y= Order_rep_rel_abund, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ 
  theme_classic()+  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid (timepoint.days.~station)

Family <- t3 %>% select(station,timepoint.days.,treatment,polymer, polymer_photo, material, detail,
                        Family, Family_rep_rel_abund, st_dev_Family_abund )%>% 
  distinct()
head(Family)
#write_csv(Family, "./analysis/family_PE462.csv")

Family <- Family %>%  mutate(Family = ifelse(Family_rep_rel_abund<0.01, "others<0.01", Family)) %>% 
  filter ( timepoint.days. %in% c("T1", "T6"))

ggplot(Family, aes(x=material, y= Family_rep_rel_abund, fill=Family))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+
  theme_classic()+ +  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid (timepoint.days.~treatment)

Genus <- t3%>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
head(Genus)
#write_csv(Genus, "./analysis/genus_PE462.csv")
#Genus <- read_csv( "../../Analysis/genus_v460.csv")
Genus <- Genus %>% mutate(Genus = ifelse(Genus_rep_rel_abund<0.01, "others<0.01", Genus))
g10_45 <- ggplot(Genus, aes(x=polymer_photo, y= Genus_rep_rel_abund, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")+ facet_grid (timepoint.days.~station)

Species <- t3 %>% select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,
polymer_photo, material, detail,  Species, Species_rep_rel_abund, st_dev_Species_abund)%>% 
  distinct()

head(Species)
#write_csv(Species, "./analysis/species_PE462.csv") 

Sp10_45 <- Species %>%  mutate(ng_species = ifelse(Species_rep_rel_abund<0.01, "others<0.01", ng_species))
 
s10_45 <- ggplot(Species, aes(x=polymer_photo, y= Species_rep_rel_abund, fill=ng_species))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+
  theme_classic()+   theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")+ facet_grid (timepoint.days.~station)

#  theme (axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold") ) 
plotly::ggplotly(s10_45)
#################

Genus_no_nc <- Genus %>% filter(Genus %nin% g) %>% filter ( timepoint.days. %in% c("T1", "T6")) %>% mutate(Genus = ifelse(Genus_rep_rel_abund<0.01, "others<0.01", Genus))

################## heatmap on genus & no NAs #############
t2_no_na_genus <- filter(t2, !Genus==" NA")%>% group_by(detail, Genus) %>% mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>% 
  select(detail, Genus, Genus_rep_rel_abund, treatment, timepoint.days., material)%>% distinct() %>% 
  arrange(desc(Genus_rep_rel_abund))

t2_no_na_genus <- filter(t2_no_na, Genus %in% (unique(t2_no_na$Genus)[1:20]))


pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(t2_no_na_genus,aes(x=detail,y=Genus,fill=Genus_rep_rel_abund))+
  geom_tile(colour="white",size=0.25)+ labs(x="",y="")+
  geom_text(aes(label = round(Genus_rep_rel_abund, 3)), colour = "Black" , size = 3)+ scale_fill_gradient(colours =  pal)+
  theme (axis.text.x = element_text(face="bold", angle=90), axis.text.y = element_text(face="bold") ) 


###############hydrocarbon_degraders###################
oil_deg <- read_lines("../../Documentation/oil_degraders.txt")
oil <-  t3 %>%  filter (Genus %in% oil_deg) %>% distinct() 
#write_csv(oil, "../../Analysis/oil_degraders_tab_v460.csv")
unique(oil$Genus)

oil_ecreme <- oil %>% select(treatment, timepoint.days., material, detail, Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
head(oil_ecreme)
#write_csv(oil_ecreme, "../../Analysis/genus_oil_degraders_v460.csv")

Genus_oil <- oil_ecreme %>%
  filter ( timepoint.days. %in% c("T1", "T6"))

ggplot(Genus_oil, aes(x=material, y= Genus_rep_rel_abund, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ 
  theme_classic()+  facet_grid (timepoint.days..~treatment)


###########beta_div####################




#######################mia############################
# convert phyloseq to TSE
BiocManager::install("mia")

TSE <- makeTreeSummarizedExperimentFromPhyloseq(physeq_object) 
