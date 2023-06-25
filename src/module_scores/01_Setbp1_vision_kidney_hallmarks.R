# setup libraries and seed
set.seed(2178)
library(Seurat)
library(VISION)
library(GSEABase)
library(dplyr)
library(here)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(tidyr)
library(readr)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# load data
load(here("data/kidney_integrated_celltypes_postSoup.Rdata"))

# set up signatures of interest
signatures <- here("data/unprocessed/mh.all.v2023.1.Mm.symbols.gmt")

# create vision object  (looking at all hallmark sets for kidney since it is less studied in SGS)
vision <- Vision(kidney_int_celltypes, signatures = signatures, pool = FALSE)

# analyze signatures on vision object 
vision <- analyze(vision) #note: this step takes the longest aside from plotting
save(vision, file = here("results/module_scores/kidney_vision_hallmarks.Rdata"))

# # get signatures
# sigscores <- getSignatureScores(vision) %>% as.data.frame(.)
# 
# # create assay object and plot
# kidney_int_celltypes[["signatures"]] <- CreateAssayObject(counts = t(sigscores))
# DefaultAssay(kidney_int_celltypes) <- "signatures"
# Idents(kidney_int_celltypes) <- "cell_type"
# new.cluster.ids <- c("B_cells",
#                      "CDIC_typeA",
#                      "CDIC_typeB",
#                      "CDPC",
#                      "DCT",
#                      "DLH",
#                      "endothelial",
#                      "fibroblasts",
#                      "LOH",
#                      "macrophages",
#                      "PCTS1",
#                      "podocytes",
#                      "PST",
#                      "pericytes",
#                      "PCTS2",
#                      "proximal_tubule",
#                      "smooth_muscle_cells")
# names(new.cluster.ids) <- levels(kidney_int_celltypes)
# kidney_int_celltypes <- RenameIdents(kidney_int_celltypes, new.cluster.ids)
# kidney_int_celltypes@meta.data[["celltype2"]] <- kidney_int_celltypes@active.ident
# #kidney_int_celltypes <- subset(kidney_int_celltypes, idents = c("URO1","URO2"),invert = T)
# kidney_int_celltypes@meta.data$celltype_condition <- paste0(kidney_int_celltypes@meta.data$celltype2, "_", kidney_int_celltypes@meta.data$type)
# Idents(kidney_int_celltypes) <- "celltype_condition"
# 
# levels(kidney_int_celltypes) <- c("B_cells_heterozygous",
#                                   "CDIC_typeA_heterozygous",
#                                   "CDIC_typeB_heterozygous",
#                                   "CDPC_heterozygous",
#                                   "DCT_heterozygous",
#                                   "DLH_heterozygous",
#                                   "endothelial_heterozygous",
#                                   "fibroblasts_heterozygous",
#                                   "LOH_heterozygous",
#                                   "macrophages_heterozygous",
#                                   "PCTS1_heterozygous",
#                                   "podocytes_heterozygous",
#                                   "PST_heterozygous",
#                                   "pericytes_heterozygous",
#                                   "PCTS2_heterozygous",
#                                   "proximal_tubule_heterozygous",
#                                   "smooth_muscle_cells_heterozygous",
#                                   "B_cells_control",
#                                   "CDIC_typeA_control",
#                                   "CDIC_typeB_control",
#                                   "CDPC_control",
#                                   "DCT_control",
#                                   "DLH_control",
#                                   "endothelial_control",
#                                   "fibroblasts_control",
#                                   "LOH_control",
#                                   "macrophages_control",
#                                   "PCTS1_control",
#                                   "podocytes_control",
#                                   "PST_control",
#                                   "pericytes_control",
#                                   "PCTS2_control",
#                                   "proximal_tubule_control",
#                                   "smooth_muscle_cells_control")
# 
# test <- AverageExpression(kidney_int_celltypes, assays = "signatures")
# fig_kid <- pheatmap(test[["signatures"]], scale = "row",cluster_cols = TRUE)
# 
# png(here("results/module_scores/kidney_heatmap_hallmarks_modulescores.png"),
#     res = 250,
#     height = 5000,
#     width = 3000
# )
# fig_kid
# dev.off()
