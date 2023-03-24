#set up libPath. NOTE: THE USER MUST CHANGE TO THEIR USERNAME
<<<<<<< Updated upstream
##set user for file path below
=======
##set user
>>>>>>> Stashed changes
username <- "jbarham3"

##give libpath
libpath <- paste0("/data/user/home/", username, "/R/x86_64-pc-linux-gnu-library/4.2") #must use this when submitting jobs through SLURM with R version 4.2, if you are not using 4.2 then the libpath needs to be changed

.libPaths(libpath)

##install packages; note this is done once in either the kidney or brain script and then commented out in the other and for all future re-runs (it takes almost 2 hours to install the packages)
#remotes::install_github(repo = "netZoo/netZooR", ref = "master")
#install.packages("data.table")
#install.packages("here")

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

#make list of expression files needed to loop through
files <- list.files(here("data/processed/expression_inputs"), pattern = "expression.Rdata", full.names = TRUE)
files <- files[c(19, 20, 23, 24, 25, 26, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40)]#grab all cerebral samples


#run for loop
for(i in files){
  #load in the expression data files needed and rename:
  fileName = i
  expression <- loadRData(fileName)
  #run panda on multi-omic inputs
  set.seed(2178)
  pandaResults <- makePanda(motif, ppi, expression)
  name <- sub(".Rdata", "", basename(i))
  save(pandaResults, file = paste0(here("results/PANDA/"), name, "_PANDA.Rdata"))
  print(paste0(name, "_PANDA network made and saved."))
}