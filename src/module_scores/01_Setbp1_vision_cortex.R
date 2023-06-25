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
load(here("data/setbp1_cerebralintcelltypes.Rdata"))

# set up signatures of interest
signatures <- here("data/unprocessed/mh.all.v2023.1.Mm.symbols.gmt")

# create vision object
vision <- Vision(cerebral_int_celltypes, signatures = signatures, pool = FALSE)


# analyze signatures on vision object 
vision <- analyze(vision) #note: this step takes the longest aside from plotting
save(vision, file = here("results/module_scores/cortex_vision_hallmarks.Rdata"))


