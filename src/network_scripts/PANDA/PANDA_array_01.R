# NOTE: be sure to run this script to set up list of files BEFORE constructing networks with 'PANDA_array_02.sh' bash script.

##Load libraries
library(here)
library(utils)

##Make list of .Rdata expression files needed for array job:
files <- list.files(here("data/processed/expression_inputs"), full.names = TRUE) # ensure that ALL exppression inputs for PANDA are in the same directory and that nothing else is in the directory.
files #save this as a text file to be the array file

write.table(files, file = here("results/array_inputs/Setbp1_PANDA_files.txt"), sep = "\t", row.names = FALSE, col.names = FALSE) # for other data user may want to change file name

##Remove " " around each file path

# Specify the input and output file paths
input_file <- here("results/array_inputs/Setbp1_PANDA_files.txt")
output_file <- here("results/array_inputs/Setbp1_PANDA_files_array.txt")

# Read the input file line by line
lines <- readLines(input_file)

# Remove quotes from each line
modified_lines <- gsub("\"", "", lines)

# Write the modified lines to the output file
writeLines(modified_lines, output_file)
