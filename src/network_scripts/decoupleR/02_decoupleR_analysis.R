#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(decoupleR)

# Fetch command line arguments
args <- commandArgs(trailingOnly = TRUE)

rdata_file <- args[1]
prior_net <- args[2]
cell_type <- args[3]
tissue <- args[4]
min_n <- as.numeric(args[5])

print(rdata_file)
print(prior_net)
print(cell_type)
print(tissue)
print(min_n)



load(file = rdata_file)
mat <- as.matrix(data.frame(decoupleR_expression[[cell_type]]))
net <- data.frame(read.csv(prior_net), header = TRUE)

# Define function
tf_activity <- function(mat, net, cell_type, min_n, tissue){
  acts <- run_mlm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = min_n) %>% mutate(cell_type = cell_type)
  #summarize_acts <- acts %>% group_by(source) %>% summarise(mean = mean(score)) %>% mutate(cell_type = cell_type)
  return(acts)
}

# Run function and save results
acts <- tf_activity(mat, net, cell_type, min_n, tissue)

# Specify output path
#output_path <- paste("/data/user/jbarham3/230227_JW_Setbp1Manuscript/results/decoupleR/cortex/", cell_type, "_acts.RData", sep="") # Note: when running cortex be sure to uncomment '/decoupleR/kidney/' for 'decoupleR/cortex/' due to duplicate cell types across tissues 
output_path <- paste("/data/user/jbarham3/230227_JW_Setbp1Manuscript/results/decoupleR/kidney/", cell_type, "_acts.RData", sep="") 
# Save the output
save(acts, file = output_path)
