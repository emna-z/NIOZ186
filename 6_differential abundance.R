library("mia")
library("MicrobiomeStat")
library(vegan)
tax <- read.csv("./new/tax_tab_physeq.csv", row.names = 1)
tax$Species <- if_else(!tax$Species=="unassigned", str_c(tax$Genus," ",tax$Species), tax$Species)
tax <- tax_table(as.matrix(tax))
otu <- otu_table(as.matrix(read.csv("./new/asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)

# map <- read.delim("./original/mapPE462.txt", na.strings = c("NA", ""))
# map$detail <- if_else(map$material=="wood", str_c(map$material,"_",map$station,"_",map$timepoint.days.), map$detail)
# map$polymer <- if_else(map$material=="wood", str_c(map$material), map$polymer)
# polymer_station <-  str_c(map$polymer, "_", map$station)
# polymer_photo <- str_c(map$polymer, "_", map$treatment)
# pol_photo_station <- str_c(map$polymer, "_", map$treatment, "_", map$station)
# map <- map %>% add_column(polymer_station) %>% add_column(polymer_photo) %>% add_column(pol_photo_station)
# map$polymer_station <- if_else(map$material=="wood", str_c(map$material,"_",map$station), map$polymer_station)
# map$polymer_photo <- if_else(map$material=="wood", str_c(map$material), map$polymer_photo)
# map$pol_photo_station <- if_else(map$material=="wood", str_c(map$material,"_",map$station), map$pol_photo_station)
# write_csv(map,"./new/map_from_tidy.csv")

map <- read.csv("./new/map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T)
map<- sample_data(map)
pseq<- merge_phyloseq(otu, tax, map)      
pseq <- pseq %>% 
  subset_taxa(!Phylum %in% c("unassigned")) %>% 
  subset_samples(station %in% c("C05", "C13")) %>%
  subset_samples(timepoint.days. %nin% c("15","0"))
  
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)
permanova <- adonis(t(otu) ~ polymer,
                    data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova$aov.tab)["polymer", "Pr(>F)"])


