#load in package libraries
library(netZooR)
library(data.table)
#library(here)

getwd() #output wd
setwd("/data/user/jbarham3/230227_JW_Setbp1Manuscript/")
.libPaths() #output libPath

#load in functions
source("/data/user/jbarham3/230227_JW_Setbp1Manuscript/src/functions/functions.R")

#enable args
args <- R.utils::commandArgs(trailingOnly = TRUE)

#load in the input data needed:
motif <- read.table(file = "/data/user/jbarham3/230227_JW_Setbp1Manuscript/data/processed/motif_inputs/mus_motif_all.txt", sep = "\t") #load in motif data
print("motif loaded")

ppi <- read.table(file = "/data/user/jbarham3/230227_JW_Setbp1Manuscript/data/processed/ppi_inputs/mm10_ppi.txt", sep = "\t") #load in ppi data
print("ppi loaded")

expression <- loadRData(args[1]) #load expression data from .Rdata in here from $SAMPLE_LIST
#expression <- readRDS(args[1]) #load expression data from .rds in here from $SAMPLE_LIST
print("expression loaded")

#run panda on multi-omic inputs
set.seed(2178)
print("seed set")

pandaResults <- makePanda(motif, ppi, expression)
name <- sub(".Rdata", "", basename(args[1]))
save(pandaResults, file = paste0("/data/user/jbarham3/230227_JW_Setbp1Manuscript/results/PANDA/", name, "_PANDA.Rdata"))
rm(pandaResults)
print(paste0(name, "_PANDA network made and saved."))