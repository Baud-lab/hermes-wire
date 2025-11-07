#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(rhdf5))
suppressMessages(library(argparse))

#Create parser object
parser <- ArgumentParser()

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-h5", "--h5_file", type="character", help="H5 input file")
parser$add_argument("-pheno", "--phenotype_version", type="character", help="Phenotype version",default="bacterial_taxa")

args <- parser$parse_args()

phenotype_version<-args$phenotype_version
h5_file<-args$h5_file


h5f = H5Fopen(h5_file)
df_col <- h5read(h5f, name = paste("phenotypes", phenotype_version, "col_header", "phenotype_ID", sep="/"))
len <- length(df_col)
cat(len)
