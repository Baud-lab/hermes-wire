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

# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Modules
parser$add_argument("-mod_ana", "--module_analysis", type="character", help="Inclusion of module: Analysis",default="Pooled")
parser$add_argument("-mod_herit", "--module_heritability", type="character", help="Inclusion of module: Heritability Analysis",default="NO")
parser$add_argument("-mod_gwas", "--module_gwas", type="character", help="Inclusion of module: GWAS Analysis",default="NO")
## Matrix processing
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-res", "--residuals", type="character", help="R object with all matrices containing read counts per taxonomic rank and study after regressing out known covariate effects")
parser$add_argument("-sid", "--sample_identifier", type="character", help="Field on the metadata file that identifies the sample",default="host_subject_id")
parser$add_argument("-hid", "--host_identifier", type="character", help="Field on the metadata file that identifies the host",default="RDIF")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
## Taxonomic
parser$add_argument("-tax_ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## Functional
parser$add_argument("-func_ranks", "--functional_ranks", type="character", help="Functional ranks to be analysed",default="EC1,EC2,EC3,EC4")
## Overall Host Genetic Effects Analyses
parser$add_argument("-pheno", "--phenotype_version", type="character", help="Phenotype version",default="bacterial_taxa")
parser$add_argument("-cage", "--cage_version", type="character", help="Name of the cage version",default="all")
parser$add_argument("-dam", "--dam_version", type="character", help="Name of the cage version",default="all")
## Heritability
parser$add_argument("-kin", "--kinship_version", type="character", help="Name of the kinship model for GRM",default="exp_DGA_7_Dec_2021")
## GWAS
parser$add_argument("-dosf", "--dosage_field", type="character", help="Metadata field containing dosage info",default="name_round8_genotypes")
parser$add_argument("-kinloco", "--kinship_version_loco", type="character", help="Name of the kinship model for GRM LOCO",default="pruned_dosages")
## Output file
parser$add_argument("-out", "--output_file", type="character", help="Name of the output h5 file",default="matrix.h5")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Modules
module_analysis<-args$module_analysis
module_heritability<-args$module_heritability
if (module_heritability != "YES" & module_heritability != "NO") {
  stop("Option for module heritability should be either 'YES' or 'NO' in capital letters.", "\n")
}
module_gwas<-args$module_gwas
if (module_gwas != "YES" & module_gwas != "NO") {
  stop("Option for module gwas should be either 'YES' or 'NO' in capital letters.", "\n")
}
## Matrix processing
input_matrix<-args$input_matrix
input_meta<-args$input_meta
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
residuals=args$residuals
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
host_identifier<-args$host_identifier
if (host_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
method=args$method

## Taxonomic
tax_ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(tax_ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}

## Functional
func_ranks<-args$functional_ranks
allowed <- c("EC1", "EC2", "EC3", "EC4")
string_vector <- strsplit(func_ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be:EC1,EC2,EC3 or EC4" ,"\n")
}

##Overall Host Genetic Effects Analyses
phenotype_version<-args$phenotype_version
cage_version<-args$cage_version
dam_version<-args$cage_version

## Heritability
kinship_version<-args$kinship_version

## GWAS
dosage_field<-args$dosage_field
kinship_version_loco<-args$kinship_version_loco

## Output file
output_file<-args$output_file

# -----------------------------------------------------------------------------


### Step 1 ####
print("Step 1: Read input files and arguments and do the custom removal of bad samples")

# Read counts data
if (input_matrix != "") {
  print("Reading counts data")
  load(input_matrix)
} else {
  stop("No counts data provided")
}

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

# Read covariates info
if (list_of_covariates != "") {
  print("Reading covariates and covariate types info")
  covs_data<-read_tab_file(list_of_covariates, TRUE)
  if (ncol(covs_data)==1){
    covs_data<-read_csv_file(list_of_covariates)
  }
  subset_ids <- covs_data$subset_id
  if (nrow(covs_data)==0) {
    stop("No subset was provided. If you only worked with one cohort/Study, please include subset 'ALL' on the input file covariates.txt and assign known covariates to it")
  }
  for (i in 1:nrow(covs_data)) {
    if (covs_data$covariates[i] == "") {
      stop(paste0("Subset ",covs_data$subset_id[i]," has no covariate assigned. Please inform at least one known covariate for this subset on the input file covariates.txt and remember that the name of the covariate should correspond with the name of one of the columns of your matadata file. If more than one covariate are assigned for this subset, remember to list them separated by comma with no spaces.")) 
    }
  }
  for (i in 1:nrow(covs_data)) {
    if (nrow(covs_data) > 1 & covs_data$subset_id[i] == "ALL" & !grepl("Study",covs_data$covariates[i])) {
      covs_data$covariates[i]=paste(covs_data$covariates[i],"Study",sep=",")
      print(paste0("Warning: Subset ",covs_data$subset_id[i]," had no covariate 'Study' assigned to it on the input file covariates.txt. Given that other subsets were listed, this covariate was automaticaly included.")) 
    }
  }
  apply(covs_data, 1,  assignCovariate) ## please check: I changed hard coded 1 for subset_id
  covnames <- unique(c(str_split(covs_data$covariates[covs_data$subset_id=="ALL"], pattern=",", simplify = TRUE)))
  total_covs=c()
  for (subset_id in subset_ids[-which(subset_ids=="ALL")]){
    covariate=get(paste("covariates",subset_id,sep="_"))
    total_covs=unique(c(total_covs,covariate))
  }
  missing_covs <- total_covs[!(total_covs %in% covnames)]
  if (length(missing_covs) == 0) {
    print("All elements in the subsets are present in the subset 'ALL'.")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
  if (list_of_covariate_types != "") {
    covs_types<-read_tab_file(list_of_covariate_types, TRUE)
    if (ncol(covs_types)==1){
      covs_types<-read_csv_file(list_of_covariate_types)
    }
    if (nrow(covs_types)==0) {
      stop("No covariate was provided to the covariate types file. Please, include at least one known covariate on the input file covariates.classification.txt")
    }
    for (i in 1:nrow(covs_types)) {
      if (covs_types$classification[i] == "") {
        stop(paste0("Subset ",covs_types$covariate[i]," has no classification assigned. Please, inform on the input file covariates.classification.txt if the covariate is 'categorical' or 'continuous'.")) 
      }
    }
    for (i in 1:length(covnames)) {
      if (length(which(covs_types$covariate==covnames[i])) == 0) {
        stop(paste0("No classification was assigned to covariate: ",covnames[i],". Please, inform on the input file covariates.classification.txt if this covariate is 'categorical' or 'continuous'"))
      }
    }
    for (subset_id in subset_ids) {
      covariates=get(paste("covariates",subset_id,sep="_"))
      covs=c()
      for (i in 1:length(covariates)){
        cat=covs_types$classification[match(covariates[i],covs_types$covariate)]
        if (cat=="categorical"){
          covs=c(covs,covariates[i])
        }
      }
      assign(paste('covariates',subset_id,sep='_'), covs)
    }
  }  else {
    stop("No covariate types file provided")
  }
} else {
  stop("No covariates file provided")
}

if (module_analysis=="Pooled"){
  subset_ids="ALL"
}

# Read ranks
print("Reading ranks")
if (phenotype_version=="microbiome_taxonomic"){
  ranks<-na.omit(c((sapply(strsplit(tax_ranks, ","),"[",1)),
                   (sapply(strsplit(tax_ranks, ","), "[",2)),
                   (sapply(strsplit(tax_ranks, ","), "[",3)),
                   (sapply(strsplit(tax_ranks, ","), "[",4)),
                   (sapply(strsplit(tax_ranks, ","), "[",5)),
                   (sapply(strsplit(tax_ranks, ","), "[",6)),
                   (sapply(strsplit(tax_ranks, ","), "[",7))))
} else {
  ranks<-na.omit(c((sapply(strsplit(func_ranks, ","),"[",1)),
                   (sapply(strsplit(func_ranks, ","), "[",2)),
                   (sapply(strsplit(func_ranks, ","), "[",3)),
                   (sapply(strsplit(func_ranks, ","), "[",4))))
}

# Read residuals
if (residuals ==""){
  stop("It seems the RData file with the residues was not generated by the previous step or not loaded properly.")
} else {
  print("Reading tables with residual values.")
  load(residuals)
  permutations <- expand.grid(subset_ids,ranks)
  concatenated <- paste(permutations$Var2, permutations$Var1, sep = "_")
  for (i in 1:length(concatenated)){
    residuals_qned_counts=residuals_qned_counts_objs[[i]]
    assign(paste('residuals_qned_counts',concatenated[i],sep='_'), residuals_qned_counts)
  }
}

### End of Step 1 ###


### Step 2: Build the HDF5 phenotype file ###
print("Step 2: : Build the HDF5 phenotype file")

# Create the phenotype group on the H5 file
print("Creating the phenotype group on the H5 file")
h5createGroup(output_file,"phenotypes")
h5createGroup(output_file,paste("phenotypes",phenotype_version,sep="/"))

# Sample IDs

# Replace samples names by host names on the counts table
metadata=metadata[metadata[[sample_identifier]] %in% rownames(count_data),]
rownames(count_data)=metadata[[host_identifier]][match(rownames(count_data),metadata[[sample_identifier]])]
h5createGroup(output_file,paste("phenotypes",phenotype_version,"row_header",sep="/"))
max_size <- max(nchar(rownames(count_data)))
h5createDataset(file=output_file,dataset=paste("phenotypes",phenotype_version,"row_header","sample_ID",sep="/"),dims=dim(count_data)[1],storage.mode='character',size=max_size)
h5write(obj=rownames(count_data),file=output_file,name=paste("phenotypes",phenotype_version,"row_header","sample_ID",sep="/"))

# Phenotype IDs
h5createGroup(output_file,paste("phenotypes",phenotype_version,"col_header",sep="/"))
max_size <- max(nchar(colnames(count_data)))
h5createDataset(file=output_file,dataset=paste("phenotypes",phenotype_version,"col_header","phenotype_ID",sep="/"),dims=dim(count_data)[2],storage.mode='character',size=max_size)
h5write(obj=colnames(count_data),file=output_file,name=paste("phenotypes",phenotype_version,"col_header","phenotype_ID",sep="/"))

# Input the matrix values
h5write(obj=count_data,file=output_file,name=paste("phenotypes",phenotype_version,"matrix",sep="/"))

# Cage IDs
print("Creating the cage group on the H5 file")
## Step 1: Create the samples vector
samples = rownames(count_data)
## Step 2: Create the vector of cages per sample
cages=c()
for (i in 1:length(samples)){
  cat("sample", i, "of", length(samples), "\n")
  m <- match(samples[i],metadata[[host_identifier]], nomatch = 0)
  if (m>0){
    cages[i] <-metadata$cage[m]
  } else {
    cages[i] <- ""
  }
}
names(cages)=samples
## Step 3: Create the "cages" folder on the H5 file
h5createGroup(output_file,"cages")
h5createGroup(output_file,paste("cages",cage_version,sep="/"))
max_size <- max(nchar(cages),na.rm=T)
h5createDataset(file=output_file, dataset=paste("cages",cage_version,"array",sep="/"),
                dims=dim(count_data)[1],storage.mode='character',size=max_size)
h5write(obj=as.character(cages), file=output_file, name=paste("cages",cage_version,"array",sep="/"))
## Step 4: Input cage IDs
max_size <- max(nchar(names(cages)),na.rm=T)
h5createDataset(file=output_file, dataset=paste("cages",cage_version,"sample_ID",sep="/"),
                dims=length(cages), storage.mode='character', size=max_size)
h5write(obj=as.character(names(cages)), file=output_file, name=paste("cages",cage_version,"sample_ID",sep="/"))

# Dam IDs
print("Creating the dam group on the H5 file")
## Step 1: Create the samples vector
samples = rownames(count_data)
## Step 2: Create the vector of dams per sample
dams=c()
for (i in 1:length(samples)){
  cat("sample", i, "of", length(samples), "\n")
  m <- match(samples[i],metadata[[host_identifier]], nomatch = 0)
  if (m>0){
    dams[i] <-metadata$dam[m]
  } else {
    dams[i] <- ""
  }
}
names(dams)=samples
## Step 3: Create the "dams" folder on the H5 file
h5createGroup(output_file,"dams")
h5createGroup(output_file,paste("dams",dam_version,sep="/"))
max_size <- max(nchar(dams),na.rm=T)
h5createDataset(file=output_file, dataset=paste("dams",dam_version,"array",sep="/"),
                dims=dim(count_data)[1],storage.mode='character',size=max_size)
h5write(obj=as.character(dams), file=output_file, name=paste("dams",dam_version,"array",sep="/"))
## Step 4: Input dam IDs
max_size <- max(nchar(names(dams)),na.rm=T)
h5createDataset(file=output_file, dataset=paste("dams",dam_version,"sample_ID",sep="/"),
                dims=length(dams), storage.mode='character', size=max_size)
h5write(obj=as.character(names(dams)), file=output_file, name=paste("dams",dam_version,"sample_ID",sep="/"))

# Subsets
print("Creating the subsets group on the H5 file")
## Step 1: Create the "subsets" folder
h5createGroup(output_file,"subsets")
## Step 2: Input the subsets
for (subset_id in subset_ids) {
  residuals = get(paste('residuals_qned_counts',ranks[1],subset_id,sep='_'))
  motch = match(colnames(residuals),metadata[[sample_identifier]])
  any(is.na(motch))
  this_metadata = metadata[motch,]
  sub_rats = this_metadata[[host_identifier]]
  sub_rats = sub_rats[!is.na(sub_rats)]
  soze = max(nchar(sub_rats))
  h5createDataset(file=output_file,dataset=paste("subsets/",subset_id,sep=''),dims=length(sub_rats),storage.mode='character',size=soze + 2)
  h5write(obj=sub_rats,file=output_file,name=paste("subsets/",subset_id,sep=''))
}

# Input the matrix values
print("Inputting the matrix values")
h5write(obj=count_data,file=output_file,name=paste("phenotypes",phenotype_version,"matrix",sep="/"))

# IMPORTANT: No co-variates are added as it uses residuals throughout

# Check the file
print("Checking the H5 file")
h5ls(output_file)

### End of Step 2 ###