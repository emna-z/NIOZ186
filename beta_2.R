library("mia")
library("MicrobiomeStat")
# install.packages("vegan")
library(vegan)
library(scater)
library("factoextra")
library(phyloseq)
set.seed(42)

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
pseq_os <- pseq %>% 
  subset_samples(station %in% c("C05"))
pseq_coast <- pseq %>% 
  subset_samples(station %in% c("C13"))
#TreeSE see https://microbiome.github.io/OMA/differential-abundance.html
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
tse_os <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq_os)
tse_coast <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq_coast)


##############CLR transform################
ps_os_clr <- microbiome::transform(pseq_os, "clr")
ps_coast_clr <- microbiome::transform(pseq_coast, "clr")
ps_clr <- microbiome::transform(pseq, "clr")
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
head(ord_clr$CA$eig)                                                  
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     


clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(pseq, ord_clr, type="samples", color="timepoint.days.", shape = "polymer") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = timepoint.days.), linetype = 2)

clr_dist_matrix <- phyloseq::distance(ps_clr, method = "bray") 
dist_matrix <- phyloseq::distance(pseq, method = "bray")
######################################
ps_clr <- microbiome::transform(pseq_coast, "clr")
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
head(ord_clr$CA$eig)                                                  
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     


clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(pseq_coast, ord_clr, type="samples", color="timepoint.days.", shape = "polymer") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = timepoint.days.), linetype = 2)
#############PERMAOVAS################
perm <- vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$c)

perm <- vegan::adonis2(dist_matrix ~ phyloseq::sample_data(pseq)$polymer)

dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(pseq)$polymer)

dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

