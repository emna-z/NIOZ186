# library("mia")
# library("MicrobiomeStat")
library(vegan)
library(scater)
library(phyloseq)
library(tidyverse)
library(Hmisc)
library("metagMisc")
library(ALDEx2)
# library(finalfit)
# library(knitr)

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

map <- read.csv("./new/map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T) %>% 
  mutate(across(.cols = c("polymer_station","pol_photo_station"), 
                ~str_replace_all(.,c("C05"= "open_water", "C13" = "Coast"))))

map<- sample_data(map)
pseq<- merge_phyloseq(otu, tax, map)      
pseq <- pseq %>% 
  subset_taxa(!Phylum %in% c("unassigned")) %>% 
  subset_samples(station %in% c("C05", "C13")) %>%
  subset_samples(timepoint.days. %nin% c("15","0"))
pseq@sam_data$timepoint.days. <- fct_relevel(pseq@sam_data$timepoint.days. ,c("5", "10", "30", "45"))

pseq@sam_data <- sample_data(data.frame(pseq@sam_data) %>% mutate(stage = if_else(timepoint.days. %in% c(5,10), "early", "late")))

pseq_glom <- tax_glom(pseq, taxrank="Genus")
pseq_c <- subset_samples(pseq, station %in% c("C13"))
pseq_c_g <- subset_samples(pseq_glom, station %in% c("C13"))
pseq_ow <- subset_samples(pseq, station %in% c("C05"))
pseq_ow_g <- subset_samples(pseq_glom, station %in% c("C05"))

##########################ALDEX2################

  aldex_pseq <- function(pseq, ALDEx2, aldex.effect, wi.eBH) {
    otu <- otu_table(pseq_ow_g)
    var <- sample_data(pseq_ow_g)$stage
    taxa_info <- data.frame(tax_table(pseq_ow_g))

    x <- aldex.clr(
      reads = otu,
      conds = var, 
      # 128 recommened for ttest, 1000 for rigorous effect size calculation
      mc.samples = 128, 
      denom = "all",
      verbose = T
    )
    
    # calculates expected values of the Welch's t-test and Wilcoxon rank test on
    # the data returned by aldex.clr
    x_tt <- aldex.ttest(
      x, 
      paired.test = FALSE, 
      verbose = T)
    # determines the median clr abundance of the feature in all samples and in
    # groups, the median difference between the two groups, the median variation
    # within each group and the effect size, which is the median of the ratio
    # of the between group difference and the larger of the variance within groups
    x_effect <- ALDEx2::aldex.effect(x, CI = TRUE, verbose = T)
    # combine all outputs 
    aldex_out <- data.frame(x_tt, x_effect)
    
    
    par(mfrow = c(1, 2))
    aldex.plot(
      aldex_out, 
      type = "MA", 
      test = "welch", 
      xlab = "Log-ratio abundance",
      ylab = "Difference",
      cutoff = 0.05
    )
    aldex.plot(
      aldex_out, 
      type = "MW", 
      test = "welch",
      xlab = "Dispersion",
      ylab = "Difference",
      cutoff = 0.05
    )
    
    res <- rownames_to_column(aldex_out, "Genus") %>%
      filter(wi.eBH <= 0.05) %>%  # here we chose the wilcoxon output rather than tt
      arrange(effect, wi.eBH)
    
    res <- left_join(res, taxa_info, "Genus")
    
    write.csv(res, paste("result_aldex_pseq_genus_glom_ow_time.csv", sep = ""))
  }
  

