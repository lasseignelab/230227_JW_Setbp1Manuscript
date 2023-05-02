# NOTE: be sure to run this script to set up list of files BEFORE constructing networks with 'PANDA_array_02.sh' bash script.

library(here)
library(utils)

#make list of expression files needed for array job:
files <- list.files(here("data/processed/expression_inputs"), full.names = TRUE) # ensure that ALL exppression inputs for PANDA are in the same directory and that nothing else is in the directory.
files #save this as a text file to be the array file

write.table(files, file = here("results/array_inputs/Setbp1_PANDA_files.txt"), sep = "\t", row.names = FALSE, col.names = FALSE) # for other data user may want to change file name
