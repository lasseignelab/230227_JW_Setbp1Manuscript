#load in package libraries
library(netZooR)
library(data.table)
library(here)

#set up environment 
setwd(here("results/PANDA/"))

#load in functions
source(here("src/functions/functions.R"))

#load in the input data needed:
motif <- read.table(file = here("data/processed/motif_inputs/mus_motif_all.txt"), sep = "\t") #load in motif data

ppi <- read.table(file = here("data/processed/ppi_inputs/mm10_ppi.txt"), sep = "\t") #load in ppi data


expression <- loadRData(args[1]) #load expression data in here from $SAMPLE_LIST
#run panda on multi-omic inputs
set.seed(2178)
pandaResults <- makePanda(motif, ppi, expression)
name <- sub(".Rdata", "", basename(i))
save(pandaResults, file = paste0(here("results/PANDA/"), name, "_PANDA.Rdata"))
rm(pandaResults)
print(paste0(name, "_PANDA network made and saved."))