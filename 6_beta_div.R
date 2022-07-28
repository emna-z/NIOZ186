library("mia")
library("MicrobiomeStat")
library(vegan)
library(scater)
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
  
# pseq.rel <- microbiome::transform(pseq, "compositional")





tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
tse_rel <- transformSamples(tse, method = "relabundance")
tse_nmds <- runNMDS(tse_rel, FUN = vegan::vegdist, method = "bray", name = "NMDS_BC", exprs_values = "relabundance")
p <- plotReducedDim(tse_nmds, "NMDS_BC", colour_by = "timepoint.days.", shape_by = "polymer", 
                    text_by = "station", point_size = 4) +
  ggtitle(label = "NMDS on overall dataset", subtitle = "distance: Bray-Curtis") 
#p
tse_mds <- runMDS(tse_rel, FUN = vegan::vegdist, name = "MDS_BC", exprs_values = "relabundance")
e <- attr(reducedDim(tse_mds, "MDS_BC"), "eig");
rel_eig <- e/sum(e[e>0])          
p1 <- plotReducedDim(tse_mds, "MDS_BC", colour_by = "timepoint.days.", shape_by = "polymer",
                     text_by = "station", point_size = 3.5) + 
  ggtitle(label = "PCoA - MDS on overall dataset", subtitle = "distance: Bray-Curtis") +
  labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

se <- runMDS(tse_rel, FUN = vegan::vegdist, name = "MDS_ait",method = "robust.aitchison", exprs_values = "counts")
e <- attr(reducedDim(se, "MDS_ait"), "eig");
rel_eig <- e/sum(e[e>0])
p2 <- plotReducedDim(se, "MDS_ait", colour_by = "timepoint.days.", shape_by = "polymer", text_by = "station") + 
  ggtitle(label = "PCoA - MDS on overall dataset", subtitle = "robust.aitchison") +
  labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

se <- runNMDS(tse_rel, FUN = vegan::vegdist, name = "NMDS_ait", method = "robust.aitchison")
p3 <- plotReducedDim(se, "NMDS_ait", colour_by = "timepoint.days.", shape_by = "polymer", text_by = "station") + 
  ggtitle(label = "NMDS on overall dataset", subtitle = "distance: robust.aitchison") 
#ggpubr::ggarrange(p,p1,p2,p3, nrow = 2, ncol = 2, common.legend = TRUE,legend = "top")



# RDA - https://microbiome.github.io/OMA/community-similarity.html
tse
tse_rel <- transformSamples(tse, method = "relabundance")
names(tse@colData@listData)
variable_names <- c("treatment","material","polymer","timepoint.days.","station")
assay <- t(assay(tse_rel, "relabundance"))
coldata <- colData(tse_rel)
formula <- as.formula(paste0("assay ~ ", str_c(variable_names, collapse = " + ")) )
rda <- rda(formula, data = coldata, scale = TRUE, na.action = na.exclude)
rda_info <- list()
variable_name <- "all"
rda_info[[variable_name]] <- c(constrained = rda$CCA$tot.chi, 
                               unconstrainded = rda$CA$tot.chi, 
                               proportion = rda$CCA$tot.chi/rda$CA$tot.chi, 
                               p_value = anova.cca(rda)["Model", "Pr(>F)"] )

# Loop through variables
permutations <- 999
for( variable_name in variable_names ){
  # Create a formula
  formula <- as.formula(paste0("assay ~ ", variable_name) )
  # Perform RDA
  rda_temp <- rda(formula, data = coldata, scale = TRUE, na.action = na.exclude)
  # Add Info to list
  rda_info[[variable_name]] <- c(constrained = rda_temp$CCA$tot.chi, 
                                 unconstrainded = rda_temp$CA$tot.chi, 
                                 proportion = rda_temp$CCA$tot.chi/rda$CA$tot.chi, 
                                 p_value = anova.cca(rda_temp, permutations = permutations
                                 )["Model", "Pr(>F)"] )
}  

# Convert into data.frame
rda_info <- t(as.data.frame(rda_info))
rda_info_clean <- rda_info
# Adjust names
colnames(rda_info_clean) <- 
  c("Explained by variables", "Unexplained by variables", "Proportion expl by vars", 
    paste0("P-value (PERMANOVA ", permutations, " permutations)") )
# Print info
# kable(rda_info_clean)

library("ggord")
library("ggplot2")

coldata <- coldata[ rownames(rda$CCA$wa), ]

# Adjust names
# Get labels of vectors
vec_lab_old <- rownames(rda$CCA$biplot)

# Loop through vector labels
vec_lab <- sapply(vec_lab_old, FUN = function(name){
  # Get the variable name
  variable_name <- variable_names[ str_detect(name, variable_names) ]
  # If the vector label includes also group name
  if( !any(name %in% variable_names) ){
    # Get the group names
    group_name <- unique( coldata[[variable_name]] )[ 
      which( paste0(variable_name, unique( coldata[[variable_name]] )) == name ) ]
    # Modify vector so that group is separated from variable name
    new_name <- paste0(variable_name, " \U2012 ", group_name)
  } else{
    new_name <- name
  }

  new_name <- expr(paste(!!new_name, " (", 
                         !!format(round( rda_info[variable_name, "proportion"]*100, 1), nsmall = 1), 
                         "%, ",italic("P"), " = ", 
                         !!gsub("0\\.","\\.", format(round( rda_info[variable_name, "p_value"], 3), 
                                                     nsmall = 3)), ")"))
  
  return(new_name)
})

names(vec_lab) <- vec_lab_old

# Create labels for axis
xlab <- paste0("RDA1 (", format(round( rda$CCA$eig[[1]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
ylab <- paste0("RDA2 (", format(round( rda$CCA$eig[[2]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")


plot1 <- ggord(rda, grp_in = coldata[["station"]], vec_lab = vec_lab,
              alpha = 0.5,
              size = 4, addsize = -4,
              #ext= 0.7, 
              txt = 3.5, repel = TRUE, 
              #coord_fix = FALSE
) + 
  # Adjust titles and labels
  guides(colour = guide_legend("station"),
         fill = guide_legend("station"),
         group = guide_legend("station"),
         shape = guide_legend("station"),
         x = guide_axis(xlab),
         y = guide_axis(ylab)) +
  theme( axis.title = element_text(size = 12) )
p5 <- plot1 + ggtitle(label = "RDA ordination on overall dataset")

tse <- agglomerateByRank(tse, rank = "Genus")
tse <- transformSamples(tse, method = "relabundance")


permanova1 <- adonis2(t(assay(tse,"relabundance")) ~ timepoint.days.,
                     by = "margin", # each term (here only 'Group') analyzed individually
                     data = colData(tse),
                     method = "euclidean",
                     permutations = 999)

dbrda <- dbrda(t(assay(tse,"relabundance")) ~ timepoint.days., 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 999)
p_values <- c( permanova1["timepoint.days.", "Pr(>F)"], permanova2["timepoint.days.", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values
# polymer
permanova1 <- adonis2(t(assay(tse,"relabundance")) ~ polymer,
                      by = "margin", # each term (here only 'Group') analyzed individually
                      data = colData(tse),
                      method = "euclidean",
                      permutations = 999)

dbrda <- dbrda(t(assay(tse,"relabundance")) ~ polymer, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 999)
p_values <- c( permanova1["polymer", "Pr(>F)"], permanova2["polymer", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values

# treatement
permanova1 <- adonis2(t(assay(tse,"relabundance")) ~ treatment,
                      by = "margin", # each term (here only 'Group') analyzed individually
                      data = colData(tse),
                      method = "robust.aitchison",
                      permutations = 999)

dbrda <- dbrda(t(assay(tse,"relabundance")) ~ treatment, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "robust.aitchison",
                        permutations = 999)
p_values <- c( permanova1["treatment", "Pr(>F)"], permanova2["treatment", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values

# polymer_photo
permanova1 <- adonis2(t(assay(tse,"relabundance")) ~ polymer_photo,
                      by = "margin", # each term (here only 'Group') analyzed individually
                      data = colData(tse),
                      method = "euclidean",
                      permutations = 999)

dbrda <- dbrda(t(assay(tse,"relabundance")) ~ polymer_photo, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 999)
p_values <- c( permanova1["polymer_photo", "Pr(>F)"], permanova2["polymer_photo", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values

# station
permanova1 <- adonis2(t(assay(tse,"relabundance")) ~ station,
                      by = "margin", # each term (here only 'Group') analyzed individually
                      data = colData(tse),
                      method = "euclidean",
                      permutations = 999)

dbrda <- dbrda(t(assay(tse,"relabundance")) ~ station, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 999)
p_values <- c( permanova1["station", "Pr(>F)"], permanova2["station", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values




# otu <- abundances(pseq.rel)
# meta <- meta(pseq.rel)
# permanova <- adonis(t(otu) ~ polymer,
#                     data = meta, permutations=99, method = "bray")
# print(as.data.frame(permanova$aov.tab)["polymer", "Pr(>F)"])
# 
# editable_graph <- dml(ggobj = fig1)
# doc <- read_pptx()
# doc <- add_slide(doc)
# doc <- ph_with(x = doc, editable_graph,
#                location = ph_location_type(type = "body") )
# print(doc, target = "./rvg.pptx")
