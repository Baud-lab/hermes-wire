#!/usr/bin/env Rscript

# Packages

suppressMessages(library(rhdf5))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(xfun))
suppressMessages(library(qqman))

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

# 2. Assign covariates to subsets
assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames <- str_split(x[2], pattern=",", simplify = TRUE) 
  assign(paste('covariates',subset_id,sep='_'), covnames, envir = .GlobalEnv)
}

# -----------------------------------------------------------------------------

#Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

## Matrix processing
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-res", "--residuals", type="character", help="R object with all matrices containing residuals")
parser$add_argument("-sid", "--sample_identifier", type="character", help="Field on the metadata file that identifies the sample",default="host_subject_id")
parser$add_argument("-hid", "--host_identifier", type="character", help="Field on the metadata file that identifies the host",default="RDIF")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
## Taxonomic
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## GWAS
parser$add_argument("-h5", "--h5_file", type="character", help="Name of the output h5 file",default="matrix.h5")
parser$add_argument("-qtls", "--qtls", type="character", help="File with SNP associations",default="All_QTLs.RData")


#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Matrix processing
input_meta<-args$input_meta
list_of_covariates<-args$list_of_covariates
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
ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}
## GWAS
h5=args$h5_file

# -----------------------------------------------------------------------------

### Step 1 ####
print("Step 1: Read input files and arguments")

# Load QTLS
print("Loading SNP associations")
load(qtls)

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
    print("All elements in the subsets are present in the subset 'ALL'.\n")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
} else {
  stop("No covariates file provided")
}

# Read taxonomic ranks
print("Reading taxonomic ranks")
ranks<-na.omit(c((sapply(strsplit(ranks, ","),"[",1)),
                 (sapply(strsplit(ranks, ","), "[",2)),
                 (sapply(strsplit(ranks, ","), "[",3)),
                 (sapply(strsplit(ranks, ","), "[",4)),
                 (sapply(strsplit(ranks, ","), "[",5)),
                 (sapply(strsplit(ranks, ","), "[",6)),
                 (sapply(strsplit(ranks, ","), "[",7))))

# Read count tables after residuals
if (residuals ==""){
  stop("Please, provide a RData file with the count tables after residuals if you are not using module 'Matrix Processing'.")
} else {
  print("Reading count tables after residuals")
  load(residuals)
  permutations <- expand.grid(subset_ids,ranks)
  concatenated <- paste(permutations$Var2, permutations$Var1, sep = "_")
  for (i in 1:length(concatenated)){
    residuals_qned_counts=residuals_qned_counts_objs[[i]]
    assign(paste('residuals_qned_counts',concatenated[i],sep='_'), residuals_qned_counts)
  }
}

# Read genotypes
print("Opening H5 file with genotypes")
h5f = H5Fopen(h5)

### End of Step 1 ###

### Step 2: Make boxplots ###
print("Step 2: : Make boxplots of residuals by genotypes for significant associations")

# Filter significant associations
print("Filtering out non-significant associations")
filt_QTLS=all_QTLs[grepl("Significant",all_QTLs$Significance),]
chrs=unique(filt_QTLS$chr)
print("Making boxplots")
pdf('genotypes_boxplot.pdf',bg='white')
genotypes=h5read(h5f, paste0("/direct/","/matrix"))
rownames(genotypes)=h5read(h5f, paste0("/direct/","/row_header/sample_ID"))
for (subset_id in subset_ids){
  filt_subset=filt_QTLS[grepl(subset_id,filt_QTLS$subset_id),]
  for (chr in chrs) {
    genotypes_chr=genotypes
    colnames(genotypes_chr)=h5read(h5f, paste0("/direct/","/col_header/chr"))
    m=which(colnames(genotypes_chr) == chr)
    genotypes_pos=genotypes
    colnames(genotypes_pos)=h5read(h5f, paste0("/direct/","/col_header/pos"))
    genotypes_pos=genotypes_pos[,m]
    filt_chrs=filt_subset[which(filt_subset$chr==chr),]
    positions=unique(filt_chrs$pos)
    for (position in positions){
      filt_pos=filt_chrs[which(filt_chrs$pos==position),]
      traits=data.frame(trait=filt_pos$measure,rank=filt_pos$rank)
      for (trait in traits$trait){
        rank=traits$rank[match(trait,traits$trait)]
        matrix=get(paste("residuals_qned_counts",rank,subset_id,sep="_"))
        colnames(matrix)=metadata[[host_identifier]][match(colnames(matrix),metadata[[sample_identifier]])]
        y=data.frame(residuals=matrix[which(rownames(matrix)==trait),])
        y$genotype=genotypes_pos[match(rownames(y),rownames(genotypes_pos)),which(colnames(genotypes_pos)==position)]
        y=na.omit(y)
        y$genotype=round(y$genotype,digit=0)
        genos=sort(unique(y$genotype))
        y$genotype[which(y$genotype==genos[1])]="AA"
        y$genotype[which(y$genotype==genos[2])]="AB"
        y$genotype[which(y$genotype==genos[3])]="BB"
        order=c("AA","AB","BB")
        myplot<-ggplot(y, aes(x=factor(genotype,level = order), y=residuals)) +
          geom_boxplot() +
          xlab(paste0("Genotypes (chr",chr,":",position,")")) +
          ylab(paste0("Relative Abundances (After CLR transformation - residuals): ",trait)) +
          ggtitle(paste0("Subset: ",subset_id)) +
          theme_bw() +
          theme(legend.position="bottom",
                axis.line = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 13),
                axis.title = element_text(size = 13),
                legend.text=element_text(size=10),
                plot.title = element_text(size=15, hjust = 0.5)) +
          geom_signif(comparisons = list(c("AA", "AB"), c("AB", "BB")),
                      map_signif_level = TRUE, textsize = 5)
        print(myplot)
      }
    }
  }
}
dev.off()
# Close genotypes
print("Closing H5 file with genotypes")
h5closeAll()
