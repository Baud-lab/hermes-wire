#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(biomformat))
suppressMessages(library(MASS))
suppressMessages(library(rhdf5))
suppressMessages(library(vegan))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(metafolio))
suppressMessages(library(philr))
suppressMessages(library(argparse))
suppressMessages(library(arrow))
suppressMessages(library(xfun))
suppressMessages(library(iNEXT.3D))
suppressMessages(library(MASS))
suppressMessages(library(doParallel))
suppressMessages(library(ggsignif))

# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Matrix processing
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-b", "--beta", type="character", help="File with beta diversity data")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Matrix processing
input_matrix<-args$input_matrix
beta<-args$beta
method=args$method

# -----------------------------------------------------------------------------

### Step 1 ####
print("Step 1: Read input files and arguments and do the custom removal of bad samples")

# Read input matrices
if (input_matrix != "") {
  print("Reading count data")
  load(input_matrix)
} else {
  stop("No count data provided")
}
if (beta != "") {
  print("Reading beta diversity data")
  load(beta)
} else {
  stop("No count data provided")
}

### End of Step 1 ###

### Step 2: Add beta values to the final count data ###
print("Step 2: Add beta values to the final count data")

count_data_beta=count_data_beta[rownames(count_data),]
count_data=cbind(count_data,count_data_beta)

# Save
save(count_data,
     file = 'count_data_with_beta.RData')

### End of Step 3 ###