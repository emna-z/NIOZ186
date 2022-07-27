rel_abund_by_ranks <- function(tidy_physeq_asv, Sample, detail) {

  tidy_physeq_asv  %>% group_by(Sample) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
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
}