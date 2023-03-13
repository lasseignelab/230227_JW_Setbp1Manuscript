#set up libPath. NOTE: THE USER MUST UNCOMMENT AND ADD THEIR USERNAME AND COMMENT OUT MINE FOR LIBPATHS
##set user
user <- "~/jbarham3"
#user <- "~/YOUR_USER" #change this

##give libpath
libpath <- "R/x86_64-pc-linux-gnu-library/slurm_4.2/"

.libpaths(paste0(user,libpath, sep = "/")) #must use this when submitting jobs through SLURM with R version 4.2, if you are not using 4.2 then the libpath needs to be changed

##install packages and load
install_github(repo = "netZoo/netZooR", ref = "master")
library(netZooR)

install.packages("data.table")
library(data.table)

install.packages("here")
library(here)

#set up environment 
setwd(here("results/PANDA/"))

#load in functions
source(here("src/functions/functions.R"))

#load in the input data needed:
motif <- read.table(file = here("data/processed/motif_inputs/mus_motif_all.txt"), sep = "\t") #load in motif data

ppi <- read.table(file = here("data/inputs/processed/ppi_inputs/mm10_ppi.txt"), sep = "\t") #load in ppi data

#make list of expression files needed to loop through
files <- list.files(here("data/processed/expression_inputs/"), pattern = "expression.Rdata", full.names = TRUE)
files <- files[c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,20,21,22,25,26,27,28,37,38,39,40,41,42,43,44,45,46)]#grab all kidney samples


#run for loop
for(i in files[2:32]){
  #load in the expression data files needed and rename:
  fileName = i
  expression <- loadRData(fileName)
  #run panda on multi-omic inputs
  set.seed(1235)
  pandaResults <- makePanda(motif, ppi, expression)
  name <- sub(".Rdata", "", basename(i))
  save(pandaResults, file = paste0(name, "_PANDA.Rdata"))
  rm(pandaResults)
  print(paste0(name, "_PANDA network made and saved."))
}