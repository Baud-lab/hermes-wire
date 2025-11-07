#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(vegan))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(philr))
suppressMessages(library(ggplot2))
suppressMessages(library(xfun))
suppressMessages(library(argparse))
suppressMessages(library(factoextra))
suppressMessages(library(cluster))

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
clean_ids<-function(oldids,metadata) {
  raw_ids<-str_split(oldids, pattern="_", simplify = TRUE)
  for (i in 1:ncol(raw_ids)){
    w=which(metadata[[sample_identifier]] == raw_ids[1,i])
    if (length(w) != 0){
      ids<-sapply(strsplit(oldids, "_"), "[",i)
    }
  }
  return(ids)
}

# 4. Data transformation (CLR)
calculate_gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}
do_clr_default = function(data){
  transform_col = function(x){
    gMean = calculate_gm_mean(x)
    transformed_col = log(x / gMean)
    return(transformed_col)
  }
  transformed_data = apply(data, MAR = 2, FUN = transform_col)
  return(transformed_data)
}

# 5. Inverse rank normalization and removal of covariate effects by regression
invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}
my_regress = function(row) {
  covariates = c(get(paste('all_covariates',subset_id,sep='_')))
  resids =  rep(NA, length(row))	
  lm1 = lm(as.formula(paste('row ~ 1 + ',paste(covariates, collapse = '+'),sep=''))) 
  resids[!is.na(row)] = residuals(lm1)
  return(resids)
}

# -----------------------------------------------------------------------------


# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Modules
parser$add_argument("-mod_ana", "--module_analysis", type="character", help="Inclusion of module: Analysis",default="Pooled")
parser$add_argument("-mod_comp", "--module_comparisons", type="character", help="Inclusion of module: Comparisons",default="NO")
## Matrix processing
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-residues", "--input_residues", type="character", help="File with input matrix: residues")
parser$add_argument("-counts", "--counts_data", type="character", help="File with counts data")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-lib", "--input_libsize", type="character", help="File with library depth info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the host",default="host_subject_id")
parser$add_argument("-phe", "--pheno", type="character", help="Define the phenotype: microbiome_taxonomic or microbiome_functional", default="microbiome_taxonomic")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-tax_ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## Functional
parser$add_argument("-func_ranks", "--functional_ranks", type="character", help="Functional ranks to be analysed",default="EC1,EC2,EC3,EC4")
## Enterotype
parser$add_argument("-ent_rank", "--enterotype_rank", type="character", help="Define the rank that will define the enterotype names", default="Genus")
parser$add_argument("-ent_kmax", "--enterotype_kmax", type="integer", help="Maximum number of enterotypes", default=5)
parser$add_argument("-ent_nstart", "--enterotype_nstart", type="integer", help="Number of initial random sets", default=25)

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

## Modules
module_analysis<-args$module_analysis
module_comparisons<-args$module_comparisons
if (module_comparisons != "YES" & module_comparisons != "NO") {
  stop("Option for module alpha diversity should be either 'YES' or 'NO' in capital letters.", "\n")
}

## Matrix processing
input_matrix<-args$input_matrix
input_residues<-args$input_residues
counts<-args$counts_data
input_meta<-args$input_meta
input_libsize<-args$input_libsize
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
method=args$method
pheno<-args$pheno

## Taxonomic
input_taxonomy<-args$input_taxonomy
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

## Enterotype
ent_rank<-args$enterotype_rank
ent_kmax<-args$enterotype_kmax
ent_nstart<-args$enterotype_nstart

# -----------------------------------------------------------------------------


### Step 1 ####
print("Step 1: Read input files and arguments")

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

# Read library depth info
print("Reading library depth info")
if (method=="16S"){
  depths<-data.frame(Sample=metadata[[sample_identifier]],Raw_Reads=metadata$Depth_16S)
} else {
  if (input_libsize != "") {
    depths<-read_tab_file(input_libsize, TRUE)
    if (module_comparisons=="YES"){
      depths$Sample<-clean_ids(depths$Sample,metadata)
    }
    rownames(metadata)=metadata[[sample_identifier]]
    metadata=metadata[depths$Sample,]
    metadata$Depth=depths$Raw_Reads
  } else {
    stop("No library depth file provided")
  }
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
      assign(paste('all_covariates',subset_id,sep='_'), c(covariates))
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

# Read taxonomy
print("Reading taxonomy info")
if (input_taxonomy != "") {
  taxonomy<-read_tab_file(input_taxonomy, FALSE)
  if (grepl("GTDB", input_taxonomy)) {
    colnames(taxonomy)=c("Genome","Taxonomy")
    taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
    taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
  } else {
    if (grepl("Greengenes2", input_taxonomy)){
      colnames(taxonomy) <- c("Genome","Taxonomy","Confidence","GTDB")
      taxonomy$Taxonomy<-gsub("\\s(?!.*\\s)", "_", as.character(taxonomy$Taxonomy), perl = TRUE)
      taxonomy$Taxonomy<-gsub(" ", "", as.character(taxonomy$Taxonomy), perl = TRUE)
      taxonomy$ASV=paste0("asv__",rownames(taxonomy))
    } else {
      stop("The taxonomy provided should be GTDB for shotgun data or Greengenes2 for 16S")
    }
  }
  taxonomy$Phylum<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",2)
  taxonomy$Class<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",3)
  taxonomy$Order<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",4)
  taxonomy$Family<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",5)
  taxonomy$Genus<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",6)
  taxonomy$Species<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",7)
} else {
  stop("No taxonomy file provided for the main dataset")
}

# Read ranks
print("Reading ranks")
if (pheno=="microbiome_taxonomic"){
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

# Read residues
if (input_residues != "") {
  print("Reading filtered matrices per subset")
  load(input_residues)
  if (module_comparisons=="YES"){
    permutations <- expand.grid(subset_ids,ranks,methods)
    colnames(permutations)=c("subset_ids","ranks","methods")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, permutations$methods,sep = "_")
  } else {
    permutations <- expand.grid(subset_ids,ranks)
    colnames(permutations)=c("subset_ids","ranks")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
  }
  concatenated=concatenated[length(concatenated)] ####Change
  for (i in 1:length(concatenated)){
    residuals_qned_counts=residuals_qned_counts_ilr_objs[[i]]
    assign(paste('residuals_qned_counts',concatenated[i],sep='_'), residuals_qned_counts)
  }
} else {
  stop("No matrices provided")
}

if (input_matrix != "") {
  print("Reading filtered matrices per subset")
  load(input_matrix)
  if (module_comparisons=="YES"){
    permutations <- expand.grid(subset_ids,ranks,methods)
    colnames(permutations)=c("subset_ids","ranks","methods")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, permutations$methods,sep = "_")
  } else {
    permutations <- expand.grid(subset_ids,ranks)
    colnames(permutations)=c("subset_ids","ranks")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
  }
  for (i in 1:length(concatenated)){
    filtered_prev=filtered_prev_objs[[i]]
    assign(paste('filtered_prev',concatenated[i],sep='_'), filtered_prev)
  }
} else {
  stop("No matrices provided")
}

# Read input matrices
if (counts != "") {
  print("Reading count data")
  load(counts)
  previous_counts=count_data
} else {
  stop("No count data provided")
}

### End of Step 1 ###

### Step 2: Enterotype estimation ###
print("Step 2: Enterotype estimation")

enterotypes_samples_objs=list()
stats_enterotypes_samples_objs=list()
prop_data_objs=list()
pdf('enterotypes.pdf',bg='white')
if (("Species" %in% ranks) || ("EC4" %in% ranks)){
  for (subset_id in subset_ids){
    # Load the matrix
    if (pheno=="microbiome_taxonomic"){
      matrix=t(get(paste('residuals_qned_counts',"Species",subset_id,sep="_")))
    } else {
      matrix=t(get(paste('residuals_qned_counts',"EC4",subset_id,sep="_")))
    }
    # Determine the optimal number of enterotypes using the silhouette analysis
    nbclust_plot <- fviz_nbclust(matrix, kmeans, method = "silhouette", k.max = ent_kmax) +
      ggtitle(paste0("Number of Clusters vs. Silhouette Width: ",subset_id)) +
      theme(plot.title = element_text(hjust = 0.5))
    plot(nbclust_plot)
    silhouette_results <- fviz_nbclust(matrix, kmeans, method = "silhouette", k.max = ent_kmax)
    silhouette=silhouette_results$data
    silhouette$abs_diff_y <- c(NA, abs(diff(silhouette$y)))
    optimal_enterotypes <- which.max(silhouette$abs_diff_y)
    # Define enterotypes
    set.seed(123)
    kmeans_results <- kmeans(matrix, centers = optimal_enterotypes, nstart = ent_nstart)
    # Plot enterotypes
    fviz_cluster_plot <- fviz_cluster(kmeans_results, data = matrix,
                                      label = NULL,
                                      geom = "point",
                                      pointsize = 3,
                                      ellipse = TRUE,
                                      show.clust.cent = TRUE,
                                      ellipse.type = "convex",
                                      ggtheme = theme_minimal()) +
      ggtitle(paste0("Cluster Plot with Ellipses: ",subset_id)) +
      theme(plot.title = element_text(hjust = 0.5))
    plot(fviz_cluster_plot)
    # Make the dataframe with the enterotypes
    enterotypes_samples=data.frame(sample=rownames(matrix),Enterotype=kmeans_results$cluster)
    enterotypes_samples$study=metadata$Study[match(enterotypes_samples$sample,metadata[[sample_identifier]])]
    enterotypes_samples$sex=metadata$Sex[match(enterotypes_samples$sample,metadata[[sample_identifier]])]
    enterotypes_samples$Enterotype=paste("enterotype",enterotypes_samples$Enterotype,sep="__")
    if (subset_id == "ALL"){
      stats_enterotypes_samples=count(enterotypes_samples,study,sex,Enterotype)
    } else {
      stats_enterotypes_samples=count(enterotypes_samples,sex,Enterotype)
    }
    stats_enterotypes_samples$Enterotype=gsub("enterotype__","Enterotype ",stats_enterotypes_samples$Enterotype)
    # Define objects to store the data
    enterotypes_samples_obj <- paste('enterotypes_samples',subset_id,sep='_')
    assign(enterotypes_samples_obj, enterotypes_samples)
    enterotypes_samples_objs<- append(enterotypes_samples_objs, list(enterotypes_samples_obj=enterotypes_samples))
    stats_enterotypes_samples_obj <- paste('stats_enterotypes_samples',subset_id,sep='_')
    assign(stats_enterotypes_samples_obj, stats_enterotypes_samples)
    stats_enterotypes_samples_objs<- append(stats_enterotypes_samples_objs, list(stats_enterotypes_samples_obj=stats_enterotypes_samples))
    
    # Enterotypes statistics (total)
    myplot <- ggplot(stats_enterotypes_samples, aes(x = sex, y = n, fill = Enterotype)) +
      geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
      #scale_fill_manual(values = colours_phylum) +
      labs(fill = "") +
      facet_wrap(~ study) +
      #ggtitle("Proportion of species with significant heritability (FDR<10%)") +
      xlab("") +
      ylab("Number of samples") +
      scale_x_discrete(guide = guide_axis(angle = 60)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=12),legend.title=element_text(size=15),)
    print(myplot)
    
    
    print(subset_id)
  }
} else {
  stop("To estimate Enterotypes, please define it on the Species level.", "\n")
}
dev.off()

# Save results
save(enterotypes_samples_objs,
     stats_enterotypes_samples_objs,
     file ="enterotypes.RData")

### End of Step 2 ###