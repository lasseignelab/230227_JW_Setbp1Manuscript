.libPaths("/data/user/jbarham3/RStudioLibs/")


library(Seurat)
library(dplyr)
library(clustree)
library(reticulate)

use_condaenv("r-reticulate")
py_install("leidenalg")
 
load("/data/user/jbarham3/jw_cellspecific/jw_cellspecific/outputs/220815_integrated_setbp1_kidney.Rdata")
load("/data/user/jbarham3/jw_cellspecific/jw_cellspecific/outputs/220815_integrated_setbp1_cerebral.Rdata")

for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5)){
  kidney.int <- FindClusters(kidney.int, graph.name = "RNA_snn", resolution = res, algorithm = 4, method = "igraph")
}

jpeg(file="/data/user/jbarham3/jw_cellspecific/jw_cellspecific/setbp1_clean/outputs/kidneyint_clustree.jpeg")
clustree(kidney.int@meta.data, prefix = "RNA_snn_res.")
dev.off()

save(kidney.int, file = "/data/user/jbarham3/jw_cellspecific/jw_cellspecific/setbp1_clean/outputs/kidney_leiden.Rdata")



for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5)){
  cerebral.int <- FindClusters(cerebral.int, graph.name = "RNA_snn", resolution = res, algorithm = 4)
}

jpeg(file="/data/user/jbarham3/jw_cellspecific/jw_cellspecific/setbp1_clean/outputs/cerebralint_clustree.jpeg")
clustree(cerebral.int@meta.data, prefix = "RNA_snn_res.")
dev.off()


save(cerebral.int, file = "/data/user/jbarham3/jw_cellspecific/jw_cellspecific/setbp1_clean/outputs/cerebral_leiden.Rdata")