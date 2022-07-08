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
#summarize_phyloseq(physeq_object)
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
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
alpha_tab <- alpha_tab %>% cbind(data.frame(physeq_object@sam_data))
#write_csv(alpha_tab, file = "./new/alpha_div_indexes_microbiome_package_with_Euk.csv")
alpha_tab <- read_csv("./new/alpha_div_indexes_microbiome_package.csv")
#phyloseq package
a_div <- phyloseq::estimate_richness(physeq_object)
#write_csv(a_div, file = "./new/alpha_div_phyloseq_with_Euk.csv")
alpha <- read_csv("./new/alpha_div_indexes_microbiome_package.csv")

#########melt data, merge replicates and calculate relative abundances ##########
source("./tidy_psmelt.R")
tidy_physeq_asv <- tidy_psmelt(physeq_object)
tidy_physeq_asv <- tidy_physeq_asv %>% filter(!Phylum == "unassigned") %>% filter(timepoint.days. %nin% c("15"))
tidy_physeq_asv$Species <- if_else(!tidy_physeq_asv$Species=="unassigned", str_c(tidy_physeq_asv$Genus," ",tidy_physeq_asv$Species), tidy_physeq_asv$Species)
tidy_physeq_asv$detail <- if_else(tidy_physeq_asv$material=="wood", str_c(tidy_physeq_asv$material,"_",tidy_physeq_asv$station,"_",tidy_physeq_asv$timepoint.days.), tidy_physeq_asv$detail)
any(tidy_physeq_asv$detail == wood_detail)

#write_csv(tidy_physeq_asv, "./new/tidyPE462_v5.csv")
tidy_physeq_asv <- read.csv("./new/tidyPE462_v5.csv")

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


####final complete tidy Table#####
polymer_station <-  str_c(t3$polymer, "_", t3$station)
polymer_photo <- str_c(t3$polymer, "_", t3$treatment)
pol_photo_station <- str_c(t3$polymer, "_", t3$treatment, "_", t3$station)

t3 <- t3 %>% 
  add_column(polymer_station, .before ="Kingdom") %>% 
  add_column(polymer_photo, .before ="Kingdom") %>% 
  add_column(pol_photo_station, .before ="Kingdom")

t3$polymer_station <- if_else(t3$material=="wood", str_c(t3$material,"_",t3$station), t3$polymer_station)
t3$polymer_photo <- if_else(t3$material=="wood", str_c(t3$material), t3$polymer_photo)
t3$pol_photo_station <- if_else(t3$material=="wood", str_c(t3$material,"_",t3$station), t3$pol_photo_station)

#write_csv(t3, "./new/tidyPE462_abundances_calc_v5.csv")  
