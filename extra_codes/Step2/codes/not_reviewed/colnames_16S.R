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
# Functions

# 1. Read gzipped or flat files
read_gzipped <- function(input_file) {
  input <- input_file
  if (file_ext(input_file)=="gz") {
    input = gzfile(input_file)  
  }
  return(input)
}

read_tab_file <- function(input_file, header=FALSE) {
  file_content <- read.delim(read_gzipped(input_file), header=header, check.names=FALSE)		
  return(file_content)
}

read_csv_file <- function(input_file) {
  file_content <- read.csv(read_gzipped(input_file), stringsAsFactors = FALSE)	
  return(file_content)
}

# 2. Assign covariates to subsets
assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames <- str_split(x[2], pattern=",", simplify = TRUE) 
  assign(paste('covariates',subset_id,sep='_'), covnames, envir = .GlobalEnv)
}

# 3. Clean ids
clean_ids<-function(oldids) {
  raw_ids<-str_split(oldids, pattern="_", simplify = TRUE)
  ids<-raw_ids[,1]
  if(!all(nchar(ids)>=10)){
    for(i in 1:length(ids)){
      if(!nchar(ids[i])>=10){ids[i] <- paste0("000", ids[i])}
    }
  }
  return(ids)
}

# -----------------------------------------------------------------------------

input_matrix_harm<-"/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/out_table_16S_taxonomic.biom"
counts_table=as.matrix(biom_data(read_biom(input_matrix_harm)))
colnames(counts_table) <- str_sub(colnames(counts_table), - 10, - 1)

input_meta="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt" ### Shallow
sample_identifier="host_subject_id"
# Extract metadata info
if (input_meta != "") {
  print("Extracting metadata info")
  metadata<-read_tab_file(input_meta, TRUE)
  if (ncol(metadata)==1){
    metadata<-read_csv_file(input_meta)
  }
  if (length(which(colnames(metadata)==sample_identifier)) == 0) {
    stop("The metadata file has no column with the same name of the sample identifier informed")
  }
} else {
  stop("No metadata file provided")
}

colnames(counts_table)<-clean_ids(colnames(counts_table))
colnames(counts_table)<-metadata[[sample_identifier]][match(colnames(counts_table),metadata$RFID)]
colnames(counts_table)[1:3]
write_biom(counts_table, "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/out_table_16S_taxonomic_v2.biom")

counts_table=as.matrix(biom_data(read_biom("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/out_table_16S_taxonomic_v2.biom")))
colnames(counts_table)[1:3]