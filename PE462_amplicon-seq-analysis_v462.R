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
library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(rmarkdown)
library(knitr)
library("DESeq2")
library("microbiome")
library(ggpubr)
library("FactoMineR")
library("factoextra")
library(usedist)
library("heatmaply")
library(Hmisc)


#####data_import######
tax <- as.matrix(read.delim("./original/representative_seq_set_tax_assignments.txt", row.names = 1, na.strings = "NA"))
tax <- tax_table(tax)
otu <- as.matrix(read.delim("./original/asv_table.txt", row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim("./original/mapPE462.txt", row.names = 1, na.strings = c("NA", "")))
physeq_object = merge_phyloseq(otu, tax, map)                 


####basic_info##############
summarize_phyloseq(physeq_object)
ntaxa(physeq_object)
nsamples(physeq_object)  ###there's 71 samples in the OTU table and 74 in the mapping file any idea why? never mind, I found in the report that 3 of them had nothing whatsoever
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #it says 7 sequences here. In the report it says 8, not much of a difference but wondering how to deal wih a sample like that

physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE) #no singletons
min(taxa_sums(physeq_object))
#physeq_object <-  subset_samples(physeq_object,(sample_sums(physeq_object) >= 1000))
max(sample_sums(physeq_object))

#####subset T3 & merge ####
#sub1 <- subset_samples(physeq_object, timepoint.days. %in% c("T1", "T6"))
#sub2 <- subset_samples(physeq_object, surface=="negative_c")
#physeq_object <- merge_phyloseq(sub1, sub2)
#physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE) #no singletons

##########getting rid of wonky taxonomy assignments ###########
get_taxa_unique(physeq_object, "Kingdom") # unassigned in Domains
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c("", "Unassigned")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Kingdom") # all good now
get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA"," NA" )) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum"))
any((get_taxa_unique(physeq_object, "Order") == "Chloroplast"))
any((get_taxa_unique(physeq_object, "Family") == "Mitochondria"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c("Chloroplast")) 
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria"))

physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE)#no singletons

#loops to redefine weird taxonomy to a single common character string "unassigned" 
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



############### alpha div ###################
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
#write_csv(alpha_tab, file = "./analysis/alpha_div_indexes_microbiome_package.csv")
metad <- data.frame(physeq_object@sam_data) 
metad$Shannon <- alpha_tab$diversity_shannon 
metad$evenness_simpson <- alpha_tab$evenness_simpson 
#m <- subset_samples(physeq_object, timepoint.days. %in% c("T1", "T6"))

#p <- ggboxplot(metad, x = "material", y = "evenness_simpson",
#               color = "material", palette =c("#5FB233FF" ,"#6A7F93FF" ,"#F57206FF" ,"#EB0F13FF", "#8F2F8BFF", "#1396DBFF"),
#              add = "jitter", shape = "treatment", size = 1) + facet_wrap(~timepoint.days.)
#p


a <- phyloseq::estimate_richness(physeq_object)
#write_csv(a, file = "./analysis/alpha_div_phyloseq.csv")
plot <- plot_richness(physeq_object, "material", "treatment", measures="Chao1")+facet_grid(treatment~timepoint.days.)
plot + geom_boxplot(data=plot$data, aes(material,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Chao1", x =NULL , y = NULL )+theme_light()

plot <- plot_richness(physeq_object, "polymer", "treatment", measures="Simpson")+facet_grid(treatment~station)
plot + geom_boxplot(data=plot$data, aes(polymer,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Simpson", x =NULL , y = NULL )+theme_light()

microbiome::plot_taxa_prevalence(physeq_object, "Phylum")+ theme(legend.position = "none") #prevalence


#########merge samples per surface all replicates together##########

#getting the phyloseq object as tidy tibble 
tidy_psmelt <- function(physeq) {
  ### INSERT Initial variable and rank name checking and modding from `psmelt`
  # Get the OTU table with taxa as rows
  rankNames = rank_names(physeq, FALSE)
  sampleVars = sample_variables(physeq, FALSE) 
  otutab <- otu_table(physeq)
  if (!taxa_are_rows(otutab)) {
    otutab <- t(otutab)
  }
  # Convert the otu table to a tibble in tidy form
  tb <- otutab %>% 
    as("matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    tidyr::gather("Sample", "Abundance", -OTU)
  # Add the sample data if it exists
  if (!is.null(sampleVars)) {
    sam <- sample_data(physeq) %>%
      as("data.frame") %>% 
      tibble::as_tibble(rownames = "Sample")
    tb <- tb %>%
      dplyr::left_join(sam, by = "Sample")
  }
  # Add the tax table if it exists
  if (!is.null(rankNames)) {
    tax <- tax_table(physeq) %>%
      as("matrix") %>%
      tibble::as_tibble(rownames = "OTU")
    tb <- tb %>%
      dplyr::left_join(tax, by = "OTU")
  }
  tb %>%
    arrange(desc(Abundance))
  # Optional conversion to a data frame doesn't affect the speed/memory usage
  # %>% as.data.frame
}
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
  
  ###correcting species by combinig genus+species#  easy fix idea merge gen/spec in new col & replace in followng paragraph  
  
  group_by(Sample, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Species) %>% 
  mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup()  

polymer_station <-  str_c(t3$polymer, "_", t3$station)
polymer_photo <- str_c(t3$polymer, "_", t3$treatment)
pol_photo_station <- str_c(polymer_photo, "_", t3$station)

t3 <- t3 %>% 
  add_column(polymer_station, .before ="Kingdom") %>% 
  add_column(polymer_photo, .before ="Kingdom") %>% 
  add_column(pol_photo_station, .before ="Kingdom")

#write_csv(t3, "./new/tidyPE462_abundances_calc_w.csv")  

#t3 <- read_csv("./analysis/.csv")

mock <- t3 %>% filter(polymer %in% c("mock_DNA"))
mock$Genus
D10_45 <- t3 %>% filter(timepoint.days.%in% c("10","45")) %>% filter(polymer %nin% c("Neg_PCR","mock_DNA" ))
#write_csv(D10_45, "./analysis/D10_45.csv")
##############tables for each rank#############
#t3 <- read_csv("./analysis/.csv")
Kingdom <- t3  %>%  select(station,timepoint.days.,treatment,polymer,pol_photo_station, polymer_station,polymer_photo,
                           material, detail,Kingdom, Kingdom_rep_rel_abund,Kingdom_st_dev_abund_samples)%>% 
  distinct() 
#write_csv(Kingdom, "./analysis/Kingdom_PE462.csv")

#Kingdom <- Kingdom %>% filter ( timepoint.days. %in% c("T1", "T6")) #%>% mutate(timepoint.days. = ifelse(material =="negative_c", "T1", timepoint.days.))

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
