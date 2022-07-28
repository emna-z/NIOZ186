library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(knitr)
library(tidyverse)
library(ANCOMBC)

set.seed(42)
##physeq
tax <- read.csv("./new/tax_tab_physeq.csv", row.names = 1)
tax$Species <- if_else(!tax$Species=="unassigned", str_c(tax$Genus," ",tax$Species), tax$Species)
tax <- tax_table(as.matrix(tax))
otu <- otu_table(as.matrix(read.csv("./new/asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)
map <- read.csv("./new/map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T)
map<- sample_data(map)
pseq<- merge_phyloseq(otu, tax, map)      
pseq <- pseq %>% 
  subset_taxa(!Phylum %in% c("unassigned")) %>% 
  subset_samples(station %in% c("C05", "C13")) %>%
  subset_samples(timepoint.days. %nin% c("15","0"))

#TreeSE see https://microbiome.github.io/OMA/differential-abundance.html
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
tse_os <- tse[ ,colData(tse)$station == "C05"]
colData(tse_os)$station <- fct_drop(colData(tse_os)$station, "C13")
tse_os <- subsetByPrevalentTaxa(tse_os, detection = 0, prevalence = 0.1)
tse_coast <- tse[ ,colData(tse)$station == "C13"]
colData(tse_coast)$station <- fct_drop(colData(tse_coast)$station, "C05")
tse_coast <- subsetByPrevalentTaxa(tse_coast, detection = 0, prevalence = 0.1)


#ALDEX2
#Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centred log-ratio transform.
# This is the input for all further analyses.
x <- aldex.clr(
  reads = assay(tse_os),
  conds = colData(tse_os)$timepoint.days., 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 999, 
  denom = "all",
  verbose = FALSE
)
# calculates expected values of the Welch's t-test and Wilcoxon rank test on
# the data returned by aldex.clr
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)
# determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)
