library(here)
library(utils)

#make list of expression files needed for array job:
files <- list.files(here("data/processed/expression_inputs"), full.names = TRUE)
files #save this as a text file to be the array file

write.table(files, file = here("results/array_inputs/Setbp1_PANDA_files.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
