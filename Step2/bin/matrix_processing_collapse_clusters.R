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
suppressMessages(library(compositions))

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
parser$add_argument("-prev", "--input_prev", type="character", help="File with input matrix after filtering by prevalence")
parser$add_argument("-not_prev", "--input_notprev", type="character", help="File with input matrix before filtering by prevalence")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-lib", "--input_libsize", type="character", help="File with library depth info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the host",default="host_subject_id")
parser$add_argument("-mfd", "--min_final_depth", type="integer", help="Minimum number of reads accepted after filtering noise", default=30000)
parser$add_argument("-mp", "--min_prevalence", type="double", help="Minimum prevalence accepted", default=0.5)
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-mets", "--methods", type="character", help="Define the methods to be compared", default="Shallow,16S")
parser$add_argument("-phe", "--pheno", type="character", help="Define the phenotype: microbiome_taxonomic or microbiome_functional", default="microbiome_taxonomic")
## Clusters
parser$add_argument("-umap", "--umap", type="character", help="File with clusters info")
parser$add_argument("-harm", "--taxa_harmonization", type="character", help="Species harmonization between individual and clustered", default="YES")
## Taxonomic
parser$add_argument("-phy", "--phylotree", type="character", help="File with phylogeny tree info")
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-tax_ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## Functional
parser$add_argument("-func_ranks", "--functional_ranks", type="character", help="Functional ranks to be analysed",default="EC1,EC2,EC3,EC4")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Modules
module_analysis<-args$module_analysis
module_comparisons<-args$module_comparisons
if (module_comparisons != "YES" & module_comparisons != "NO") {
  stop("Option for module alpha diversity should be either 'YES' or 'NO' in capital letters.", "\n")
}

## Matrix processing
input_notprev<-args$input_notprev
input_prev<-args$input_prev
input_meta<-args$input_meta
input_libsize<-args$input_libsize
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
min_final_depth<-args$min_final_depth
min_prevalence<-args$min_prevalence
method=args$method
methods<-args$methods
pheno<-args$pheno

##Clusters
umap_dir<-args$umap
harm=args$taxa_harmonization

## Taxonomic
phylotree<-args$phylotree
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

# Hierarchical tree info
if ((pheno=="microbiome_taxonomic")){
  print("Reading phylogenetic info")
  if (phylotree == ""){
    stop("Please, provide the phylogeny if you want to estimate the phylogenetic diversity (PD)")
  } else {
    tree = ape::read.tree(phylotree)
  } 
} else {
  print("Preparing a functional hierarquical tree")
  ec_numbers=c()
  for (subset_id in subset_id){
    counts_table_EC4=get(paste('filtered_prev_EC4',subset_id,sep='_'))
    numbers <- sapply(strsplit(rownames(counts_table_EC4), "_"),"[",3)
    ec_numbers=unique(c(ec_numbers,numbers))
  }
  ec_components <- data.frame(do.call(rbind, strsplit(ec_numbers, "\\.")))
  colnames(ec_components) <- c("L1", "L2", "L3", "L4")
  ec_components$Level1 <- ec_components$L1
  ec_components$Level2 <- with(ec_components, paste(Level1, L2, sep = "."))
  ec_components$Level3 <- with(ec_components, paste(Level2, L3, sep = "."))
  ec_components$Level4 <- with(ec_components, paste(Level3, L4, sep = "."))
  ec_components=ec_components[,5:8]
  rownames(ec_components)=paste0("EC4__",ec_components$Level4)
  tree <- as.phylo(hclust(dist(ec_components)))
}

# Read methods
print("Reading methods")
methods<-na.omit(c((sapply(strsplit(methods, ","),"[",1)),
                   (sapply(strsplit(methods, ","), "[",2))))
if (methods[1]!=method){
  stop("ERROR: The 1st method pointed on the 'methods' field should be the main method to be used for profiling.")
}

# Read input matrices
## Before filtering by prevalence
if (input_notprev != "") {
  print("Reading filtered matrices per subset")
  load(input_notprev)
  if (module_comparisons=="YES"){
    permutations <- expand.grid(subset_ids,methods)
    colnames(permutations)=c("subset_ids","methods")
    concatenated <- paste(permutations$subset_ids,permutations$methods, sep = "_")
  } else {
    permutations <- expand.grid(subset_ids)
    colnames(permutations)=c("subset_ids")
    concatenated <- paste(permutations$subset_ids, sep = "_")
  }
  for (i in 1:length(concatenated)){
    filtered_trait=filtered_trait_objs[[i]]
    assign(paste('filtered_trait',concatenated[i],sep='_'), filtered_trait)
  }
} else {
  stop("No matrices provided")
}

## After filtering by prevalence
if (input_prev != "") {
  print("Reading filtered matrices per subset")
  load(input_prev)
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

# Read clusters information
if (umap_dir != "") {
  print("Reading clusters information")
  umap_files <- data.frame(file=list.files(umap_dir, pattern=".csv"))
  umap_files$subset=gsub(".csv","",umap_files$file)
  umap_files$subset= sapply(strsplit(umap_files$subset, "_"), "[",4)
  for (i in 1:length(umap_files$file)){
    clusters=read_csv_file(paste(umap_dir,umap_files$file[i],sep="/"))
    colnames(clusters)[1]="Genome"
    if (module_comparisons=="YES"){
      matrix=get(paste("filtered_trait",umap_files$subset[i],method,sep="_"))
    } else {
      matrix=get(paste("filtered_trait",umap_files$subset[i],sep="_"))
    }
    clusters$Genome=rownames(matrix)
    clusters$Cluster=gsub("-1","noise",clusters$Cluster)
    clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")
    assign(paste('clusters',umap_files$subset[i],sep='_'), clusters)
  }
} else {
  stop("No clusters info provided")
}

### End of Step 1 ###

### Step 2: Collapsing tables based on clusters ###
print("Step 2: Collapse tables based on clusters, remove non-clustered species and low-prevalence clusters")

filtered_prev_cluster_objs <- list()
if (harm=="YES"){
  filtered_trait_objs <- list()
}

for (subset_id in subset_ids) {
  
  if (module_comparisons=="YES"){
    matrix = get(paste('filtered_trait',subset_id,method,sep='_'))
  } else {
    matrix = get(paste('filtered_trait',subset_id,sep='_'))
  }
  print(paste0(subset_id,": ",dim(matrix)," (Before clustering)"))
  cluster = get(paste('clusters',subset_id,sep='_'))
  # Sum clusters' counts
  matrix=t(matrix[rownames(matrix) %in% cluster$Genome,])
  colnames(matrix)=as.character(cluster$Cluster)[match(colnames(matrix), cluster$Genome)]
  matrix=t(matrix)
  matrix=t(sapply(by(matrix,rownames(matrix),colSums),identity))
  print(paste0(subset_id,": ",dim(matrix)," (After clustering)"))
  # Remove non-clustered species (noise)
  matrix=matrix[!grepl("noise",rownames(matrix)),]
  print(paste0(subset_id,": ",dim(matrix)," (After removing non-clustered species)"))
  
  ## Remove low depth samples
  #depth <- colSums(matrix)
  #names(depth) <- colnames(matrix)
  #samp.keep <- colnames(matrix)[depth >= min_final_depth]
  #matrix <- matrix[,samp.keep]
  #print(paste0(subset_id,": ",dim(matrix)," (After filtering out low depth samples)"))
  
  # Remove low-prevalence clusters
  spe.nonZero  <- rowSums(matrix>0)/ncol(matrix) # Prevalence
  spe.keep <- rownames(matrix)[spe.nonZero > min_prevalence]
  length(spe.keep)
  matrix <- matrix[spe.keep,]
  print(paste0(subset_id,": ",dim(matrix)," (After filtering out low prevalence clusters)"))
  
  # Define objects to store the data
  filtered_prev_cluster_obj <- paste('filtered_prev_Clusters',subset_id,sep='_')
  assign(filtered_prev_cluster_obj, matrix)
  filtered_prev_cluster_objs<- append(filtered_prev_cluster_objs, list(filtered_prev_cluster_obj=matrix))

  if (pheno=="microbiome_taxonomic" && harm=="YES"){
    
    print("Harmonizing with taxonomic tables")
    # Filter non-clustered species from the taxonomic tables
    cluster=cluster[!grepl("noise",cluster$Cluster),]
    if (module_comparisons=="YES"){
      matrix = get(paste('filtered_trait',subset_id,method,sep='_'))
    } else {
      matrix = get(paste('filtered_trait',subset_id,sep='_'))
    }
    print(paste0(subset_id,": ",dim(matrix)," (Before removing non-clustered species)"))
    matrix=matrix[rownames(matrix) %in% cluster$Genome,]
    print(paste0(subset_id,": ",dim(matrix)," (After removing non-clustered species)"))
    
    # Define objects to store the data
    filtered_trait_obj <- paste('filtered_trait',subset_id,sep='_')
    assign(filtered_trait_obj, matrix)
    filtered_trait_objs<- append(filtered_trait_objs, list(filtered_trait_obj=matrix))
    
    # Collapse ranks
    for (rank in ranks){
      this_matrix = t(matrix)
      colnames(this_matrix)=taxonomy[[rank]][match(colnames(this_matrix), taxonomy$Genome)]
      collapsed_matrix=t(this_matrix)
      collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
      print(paste0(rank,"/", subset_id,": ",dim(collapsed_matrix)," (After collapsing rank)"))
      
      ## Remove low depth samples
      #depth <- colSums(collapsed_matrix)
      #names(depth) <- colnames(collapsed_matrix)
      #samp.keep <- colnames(collapsed_matrix)[depth >= min_final_depth]
      #filt.matrix <- collapsed_matrix[,samp.keep]
      #print(paste0(rank,"/", subset_id,": ",dim(filt.matrix)," (After filtering out low depth samples)"))
      
      # Remove low-prevalence taxa
      spe.nonZero  <- rowSums(collapsed_matrix>0)/ncol(collapsed_matrix) # Prevalence
      spe.keep <- rownames(collapsed_matrix)[spe.nonZero > min_prevalence]
      length(spe.keep)
      collapsed_matrix <- collapsed_matrix[spe.keep,]
      print(paste0(rank,"/", subset_id,": ",dim(collapsed_matrix)," (After filtering out low prevalence taxa)"))     
      
      # Define objects to store the data
      assign(paste('filtered_prev',rank,subset_id,sep='_'), collapsed_matrix)
    }
  }
}

filtered_prev_objs <- list()
for (rank in ranks){
  for (subset_id in subset_ids){
    if (module_comparisons=="YES" & harm!="YES") {
      matrix = get(paste('filtered_prev',rank,subset_id,method,sep='_'))
    } else {
      matrix = get(paste('filtered_prev',rank,subset_id,sep='_'))
    }
    filt.matrix=get(paste('filtered_prev',rank,subset_id,sep='_'))
    filtered_prev_obj <- paste('filtered_prev',rank,subset_id,sep='_')
    assign(filtered_prev_obj, filt.matrix)
    filtered_prev_objs<- append(filtered_prev_objs, list(filtered_prev_obj=filt.matrix))
  }
}

# Save results for posterior analyses

save(filtered_prev_cluster_objs,
     file="filtered_prev_clusters.RData")
save(filtered_prev_objs,
     file="filtered_prev_new.RData")
save(filtered_trait_objs,
     file="filtered_traits_matrix_noprev_new.RData")

### End of Step 2 ###

### Step 3: Relative abundance plots ###
print("Step 3.1: Cluster relative abundance plots")
# Define colors
total_cluster=c()
for (subset_id in subset_ids){
  cluster_matrix <- get(paste("filtered_prev_Clusters",subset_id,sep="_"))
  total_cluster=unique(c(total_cluster,rownames(cluster_matrix)))
}
colours = metafolio::gg_color_hue(n =  length(total_cluster), c = 100, l = 65)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]
# Make the plot
pdf('cluster_stacked_barplots.pdf')
for (subset_id in subset_ids){
  # Make into relative abundance
  cluster_matrix <- apply(get(paste("filtered_prev_Clusters",subset_id,sep="_")), 2, function(i) i/sum(i))
  # Meanwhile, transpose the count table for future wrangling.
  cluster_matrix <- data.frame(t(cluster_matrix))
  # No covariates
  colnames(cluster_matrix)=gsub("cluster__","Guild ",colnames(cluster_matrix))
  cluster_matrix$Sample=rownames(cluster_matrix)
  # Wrangle the data to long format for easy plotting
  barlong = cluster_matrix %>%
    pivot_longer(!Sample, names_to = c("Guild"), values_to = "value") %>%
    mutate(Guild = str_replace(Guild, ".*_or_", ""))
  # Calculate mean values across covariates
  grouped_matrix <- barlong %>%
    group_by(Guild) %>%
    summarise(mean_value_per_Cluster = mean(value, na.rm = TRUE), .groups = "drop")
  myplot <- grouped_matrix %>%
    ggplot(aes(x = factor(1), y = mean_value_per_Cluster, fill = Guild)) +
    geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
    scale_fill_manual(values = colours) +
    labs(color = "") +
    ggtitle(paste0("Subset: ",subset_id)) +
    xlab("") +
    ylab("Mean relative abundance") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text=element_text(size=10),)
  print(myplot)
  # Prepare the data for ggplot by adding in metadata here
  motch = match(rownames(cluster_matrix), metadata[[sample_identifier]])
  any(is.na(motch))
  this_metadata = metadata[motch,]
  covariates=get(paste("covariates",subset_id,sep="_"))
  formula <- c("Sample")
  for (covariate in covariates){
    cluster_matrix[[covariate]]=this_metadata[[covariate]]
    cluster_matrix$Sample=rownames(cluster_matrix)
    # Wrangle the data to long format for easy plotting
    formula <- c(formula, covariate)
    barlong = cluster_matrix %>%
      pivot_longer(!formula, names_to = c("Guild"), values_to = "value") %>%
      mutate(Guild = str_replace(Guild, ".*_or_", ""))
    # Calculate mean values across covariates
    grouped_matrix <- barlong %>%
      group_by(Guild, !!sym(covariate)) %>%
      summarise(mean_value_per_cluster = mean(value))
    myplot <- ggplot(grouped_matrix, aes(x = !!sym(covariate), y = mean_value_per_cluster, fill = Guild)) +
      geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
      scale_fill_manual(values = colours) +
      labs(color = "") +
      ggtitle(paste0("Subset: ",subset_id)) +
      xlab(paste0("Covariate: ",covariate)) +
      ylab("Mean relative abundance") +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),)
    print(myplot)
  }
}
dev.off()

pdf('stacked_barplots_new.pdf')
if (pheno=="microbiome_taxonomic"){
  print("Step 3.2: Phylum relative abundance plots")
  if ("Phylum" %in% ranks) {
    # Define colors
    total_phyla=c()
    for (subset_id in subset_ids){
      phylum_matrix <- get(paste("filtered_prev_Phylum",subset_id,sep="_"))
      total_phyla=unique(c(total_phyla,rownames(phylum_matrix)))
    }
    total_phyla=c(total_phyla,"Rare Taxa")
    colours = metafolio::gg_color_hue(n =  length(total_phyla), c = 100, l = 65)
    shuffled_indices <- sample(length(colours))
    colours <- colours[shuffled_indices]
    # Make the plot
    for (subset_id in subset_ids){
      # Make into relative abundance
      phylum_matrix <- apply(get(paste("filtered_prev_Phylum",subset_id,sep="_")), 2, function(i) i/sum(i))
      colSums(phylum_matrix)[1:5]
      # Define a cutoff for rare taxa
      maxabundances <- apply(phylum_matrix, 1, max)
      # Meanwhile, transpose the count table for future wrangling.
      phylum_matrix <- data.frame(t(phylum_matrix))
      # For every sample, sum up all rare taxa ( < 1% at their highest in this case)
      rare=which(maxabundances < 0.01)
      if (length(rare)>1){
        phylum_matrix$`Rare Taxa` <- rowSums(phylum_matrix[,maxabundances < 0.01], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
        phylum_matrix = phylum_matrix[,c(maxabundances > 0.01, T) ] #`T` to include the `Rare Taxa`.
      } else {
        colnames(phylum_matrix)[which(colnames(phylum_matrix)==names(maxabundances)[rare])]="Rare Taxa"
      }
      # No covariates
      colnames(phylum_matrix)=gsub("p__","",colnames(phylum_matrix))
      phylum_matrix$Sample=rownames(phylum_matrix)
      # Wrangle the data to long format for easy plotting
      barlong = phylum_matrix %>%
        pivot_longer(!Sample, names_to = c("Phylum"), values_to = "value") %>%
        mutate(Phylum = str_replace(Phylum, ".*_or_", ""))
      # Calculate mean values across covariates
      grouped_matrix <- barlong %>%
        group_by(Phylum) %>%
        summarise(mean_value_per_Phylum = mean(value, na.rm = TRUE), .groups = "drop")
      myplot <- grouped_matrix %>%
        ggplot(aes(x = factor(1), y = mean_value_per_Phylum, fill = Phylum)) +
        geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
        scale_fill_manual(values = colours) +
        labs(color = "") +
        ggtitle(paste0("Subset: ",subset_id)) +
        xlab("") +
        ylab("Mean relative abundance") +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text=element_text(size=10),)
      print(myplot)
      # Prepare the data for ggplot by adding in metadata here
      motch = match(rownames(phylum_matrix), metadata[[sample_identifier]])
      any(is.na(motch))
      this_metadata = metadata[motch,]
      covariates=get(paste("covariates",subset_id,sep="_"))
      formula <- c("Sample")
      for (covariate in covariates){
        phylum_matrix[[covariate]]=this_metadata[[covariate]]
        phylum_matrix$Sample=rownames(phylum_matrix)
        # Wrangle the data to long format for easy plotting
        formula <- c(formula, covariate)
        barlong = phylum_matrix %>%
          pivot_longer(!formula, names_to = c("Phylum"), values_to = "value") %>%
          mutate(Phylum = str_replace(Phylum, ".*_or_", ""))
        # Calculate mean values across covariates
        grouped_matrix <- barlong %>%
          group_by(Phylum, !!sym(covariate)) %>%
          summarise(mean_value_per_Phylum = mean(value))
        myplot <- grouped_matrix %>%
          ggplot(aes(x = !!sym(covariate), y = mean_value_per_Phylum, fill = Phylum)) +
          geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
          scale_fill_manual(values = colours) +
          labs(color = "") +
          ggtitle(paste0("Subset: ",subset_id)) +
          xlab(paste0("Covariate: ",covariate)) +
          ylab("Mean relative abundance") +
          scale_x_discrete(guide = guide_axis(angle = 45)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5),
                axis.line = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 15),
                legend.text=element_text(size=10),)
        print(myplot)
      }
    }
  } else {
    print("WARNING: As 'Phylum' was not selected as a rank no Phylum relative abundance plots will be generated")
  }
  dev.off()
} else {
  print("Step 3.2: EC1 relative abundance plots")
  if ("EC1" %in% ranks) {
    # Define colors
    total_functions=c()
    for (subset_id in subset_ids){
      function_matrix <- get(paste("filtered_prev_EC1",subset_id,sep="_"))
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__1")]="Oxidoreductases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__2")]="Transferases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__3")]="Hydrolases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__4")]="Lyases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__5")]="Isomerases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__6")]="Ligases"
      rownames(function_matrix)[which(rownames(function_matrix)=="EC1__7")]="Translocases"
      total_functions=unique(c(total_functions,rownames(function_matrix)))
    }
    total_functions=c(total_functions,"Rare Functions")
    colours = metafolio::gg_color_hue(n =  length(total_functions), c = 100, l = 65)
    shuffled_indices <- sample(length(colours))
    colours <- colours[shuffled_indices]
    # Make the plot
    for (subset_id in subset_ids){
      # Make into relative abundance
      function_matrix <- apply(get(paste("filtered_prev_EC1",subset_id,sep="_")), 2, function(i) i/sum(i))
      colSums(function_matrix)[1:5]
      # Define a cutoff for Rare Functions
      maxabundances <- apply(function_matrix, 1, max)
      # Meanwhile, transpose the count table for future wrangling.
      function_matrix <- data.frame(t(function_matrix))
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__1")]="Oxidoreductases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__2")]="Transferases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__3")]="Hydrolases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__4")]="Lyases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__5")]="Isomerases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__6")]="Ligases"
      colnames(function_matrix)[which(colnames(function_matrix)=="EC1__7")]="Translocases"
      # Prepare the data for ggplot by adding in metadata here
      motch = match(rownames(function_matrix), metadata[[sample_identifier]])
      any(is.na(motch))
      this_metadata = metadata[motch,]
      covariates=get(paste("covariates",subset_id,sep="_"))
      formula <- c("Sample")
      for (covariate in covariates){
        function_matrix[[covariate]]=this_metadata[[covariate]]
        function_matrix$Sample=rownames(function_matrix)
        # Wrangle the data to long format for easy plotting
        formula <- c(formula, covariate)
        barlong = function_matrix %>%
          pivot_longer(!formula, names_to = c("EC1"), values_to = "value") %>%
          mutate(EC1 = str_replace(EC1, ".*_or_", ""))
        # Calculate mean values across covariates
        grouped_matrix <- barlong %>%
          group_by(EC1, !!sym(covariate)) %>%
          summarise(mean_value_per_EC1 = mean(value))
        myplot <- grouped_matrix %>%
          ggplot(aes(x = !!sym(covariate), y = mean_value_per_EC1, fill = EC1)) +
          geom_bar(stat = "identity", position = "stack", col = "white", linewidth = 2, width = 1) +
          scale_fill_manual(values = colours,
                            labels = total_functions,
                            limits = total_functions) +
          labs(color = "EC1") +
          ggtitle(paste0("Subset: ",subset_id)) +
          xlab(paste0("Covariate: ",covariate)) +
          ylab("Mean relative abundance") +
          scale_x_discrete(guide = guide_axis(angle = 45)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5),
                axis.line = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 15),
                legend.text=element_text(size=10),)
        print(myplot)
      }
    }
  } else {
    print("WARNING: As 'EC1' was not selected as a rank no EC1 relative abundance plots will be generated")
  }
  dev.off()
}


### End of Step 3 ###

### Step 4: Centered log-ratio (CLR) transformation (Important for Differential Abundance Analyses) ###
print("Step 4.1: Clusters centered log-ratio (CLR) transformation")

offset = 0.00001
filtered_clr_counts_cluster_objs <- list()
for (subset_id in subset_ids) {
  matrix.filtered = get(paste('filtered_prev_Clusters',subset_id,sep='_'))
  matrix.filtered = matrix.filtered + offset
  clr_counts = do_clr_default(matrix.filtered)
  filtered_clr_counts_cluster_obj <- paste('filtered_clr_counts',subset_id,sep='_')
  assign(filtered_clr_counts_cluster_obj, clr_counts)
  #Define objects to store the data
  filtered_clr_counts_cluster_objs<- append(filtered_clr_counts_cluster_objs, list(filtered_clr_counts_cluster_obj=clr_counts))
  print(paste0("Subset:", subset_id))
}

print("Step 4.2: Ranks centered log-ratio (CLR) transformation")
offset = 0.00001
filtered_clr_counts_objs <- list()
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix.filtered = get(paste('filtered_prev',rank,subset_id,sep='_'))
    matrix.filtered = matrix.filtered + offset
    clr_counts = do_clr_default(matrix.filtered)
    filtered_clr_counts_obj <- paste('filtered_clr_counts',rank,subset_id,sep='_')
    assign(filtered_clr_counts_obj, clr_counts)
    #Define objects to store the data
    filtered_clr_counts_objs<- append(filtered_clr_counts_objs, list(filtered_clr_counts_obj=clr_counts))
    print(paste0(rank,":", subset_id))
  }
}

print("Step 4.3: Ranks isometric log-ratio (ILR) transformation")
offset = 1
filtered_ilr_counts_objs <- list()
for (rank in ranks[length(ranks)]){
  for (subset_id in subset_ids) {
    if (module_comparisons=="YES") {
      matrix.filtered = get(paste('filtered_prev',rank,subset_id,method,sep='_'))
    } else {
      matrix.filtered = get(paste('filtered_prev',rank,subset_id,sep='_'))
    }
    matrix.filtered=(round(matrix.filtered, digits = 0))
    if (pheno=="microbiome_taxonomic"){
      rownames(matrix.filtered)=taxonomy$Genome[match(rownames(matrix.filtered),taxonomy$Species)]
    }
    # Prune the phylogenetic tree
    keep_taxa=rownames(matrix.filtered)
    remove_taxa = setdiff(tree$tip.label, keep_taxa)
    pruned_tree = drop.tip(tree, remove_taxa)
    pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
    matrix.filtered=matrix.filtered[rownames(matrix.filtered) %in% pruned_tree$tip.label,]
    matrix.filtered = matrix.filtered + offset
    ilr_counts = t(as.matrix(philr(t(matrix.filtered), pruned_tree, 
                                   part.weights='enorm.x.gm.counts', 
                                   ilr.weights='blw.sqrt')))
    filtered_ilr_counts_obj <- paste('filtered_ilr_counts',rank,subset_id,sep='_')
    assign(filtered_ilr_counts_obj, ilr_counts)
    #Define objects to store the data
    filtered_ilr_counts_objs<- append(filtered_ilr_counts_objs, list(filtered_ilr_counts_obj=ilr_counts))
    print(paste0(rank,":", subset_id))
  }
}

# Save filtered clr_counts in a list of objects for posterior analyses
save(filtered_clr_counts_cluster_objs,
     file = 'filtered_clr_counts_clusters.RData')
save(filtered_clr_counts_objs,
     filtered_ilr_counts_objs,
     file = 'filtered_clr_counts_new.RData')

### End of Step 4 ###

### Step 5: Rank normalization and removal of known covariate effects by regression ###
print("Step 5.1: Cluster normalization and removal of known covariates effects by regression")
residuals_qned_counts_clusters_objs<-list()
for (subset_id in subset_ids) {
  matrix = get(paste('filtered_clr_counts',subset_id,sep='_'))
  if (any(!colnames(matrix) %in% metadata[[sample_identifier]])) {stop('wherzit')}
  motch = match(colnames(matrix), metadata[[sample_identifier]])
  this_metadata = metadata[motch,]
  Study_cov = this_metadata$Study
  Study_cov = factor(Study_cov)
  Study_cov = droplevels(Study_cov)
  cc = complete.cases(Study_cov)
  for (covname in covnames){
    cov = this_metadata[,covname]
    if(any(is.na(cov))) {stop(covname)}
    cc = complete.cases(cc,cov)
    assign(paste('cov',covname,sep='_'),cov)
  }
  for (covname in covnames){
    cov = get(paste('cov',covname,sep='_'))
    cov = cov[cc]
    assign(covname,cov)
  }
  #Study_cov = Study_cov[cc]
  names=rownames(matrix)
  matrix = matrix[,cc]
  rownames(matrix)=names
  qned_counts = t(apply(matrix, FUN = invrank, MAR = 1))
  residuals_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1))
  colnames(residuals_qned_counts) = colnames(matrix)
  residuals_qned_counts_clusters_obj <- paste('residuals_qned_counts',subset_id,sep='_')
  assign(residuals_qned_counts_clusters_obj, residuals_qned_counts)
  #Define objects to store the data
  residuals_qned_counts_clusters_objs<- append(residuals_qned_counts_clusters_objs, list(residuals_qned_counts_clusters_obj=residuals_qned_counts))
  print(paste0("Subset:", subset_id))
}

print("Step 5.2: Rank normalization and removal of known covariates effects by regression")

print("ILR")
residuals_qned_counts_ilr_objs<-list()
for (rank in ranks[length(ranks)]){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered_ilr_counts',rank,subset_id,sep='_'))
    if (any(!colnames(matrix) %in% metadata[[sample_identifier]])) {stop('wherzit')}
    motch = match(colnames(matrix), metadata[[sample_identifier]])
    this_metadata = metadata[motch,]
    Study_cov = this_metadata$Study
    Study_cov = factor(Study_cov)
    Study_cov = droplevels(Study_cov)
    cc = complete.cases(Study_cov)
    for (covname in covnames){
      cov = this_metadata[,covname]
      if(any(is.na(cov))) {stop(covname)}
      cc = complete.cases(cc,cov)
      assign(paste('cov',covname,sep='_'),cov)
    }
    for (covname in covnames){
      cov = get(paste('cov',covname,sep='_'))
      cov = cov[cc]
      assign(covname,cov)
    }
    #Study_cov = Study_cov[cc]
    names=rownames(matrix)
    matrix = matrix[,cc]
    rownames(matrix)=names
    qned_counts = t(apply(matrix, FUN = invrank, MAR = 1))
    residuals_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1))
    colnames(residuals_qned_counts) = colnames(matrix)
    residuals_qned_counts_obj <- paste('residuals_qned_counts',rank,subset_id,sep='_')
    assign(residuals_qned_counts_obj, residuals_qned_counts)
    #Define objects to store the data
    residuals_qned_counts_ilr_objs<- append(residuals_qned_counts_ilr_objs, list(residuals_qned_counts_obj=residuals_qned_counts))
    print(paste0(rank,":", subset_id))
  }
}

print("CLR")
residuals_qned_counts_objs<-list()
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered_clr_counts',rank,subset_id,sep='_'))
    if (any(!colnames(matrix) %in% metadata[[sample_identifier]])) {stop('wherzit')}
    motch = match(colnames(matrix), metadata[[sample_identifier]])
    this_metadata = metadata[motch,]
    Study_cov = this_metadata$Study
    Study_cov = factor(Study_cov)
    Study_cov = droplevels(Study_cov)
    cc = complete.cases(Study_cov)
    for (covname in covnames){
      cov = this_metadata[,covname]
      if(any(is.na(cov))) {stop(covname)}
      cc = complete.cases(cc,cov)
      assign(paste('cov',covname,sep='_'),cov)
    }
    for (covname in covnames){
      cov = get(paste('cov',covname,sep='_'))
      cov = cov[cc]
      assign(covname,cov)
    }
    #Study_cov = Study_cov[cc]
    names=rownames(matrix)
    matrix = matrix[,cc]
    rownames(matrix)=names
    qned_counts = t(apply(matrix, FUN = invrank, MAR = 1))
    residuals_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1))
    colnames(residuals_qned_counts) = colnames(matrix)
    residuals_qned_counts_obj <- paste('residuals_qned_counts',rank,subset_id,sep='_')
    assign(residuals_qned_counts_obj, residuals_qned_counts)
    #Define objects to store the data
    residuals_qned_counts_objs<- append(residuals_qned_counts_objs, list(residuals_qned_counts_obj=residuals_qned_counts))
    print(paste0(rank,":", subset_id))
  }
}

# Save filtered clr_counts in a list of objects for posterior  correlations with alpha_diversity
save(residuals_qned_counts_clusters_objs,
     file = 'residuals_qned_counts_clusters.RData')
save(residuals_qned_counts_objs,
     residuals_qned_counts_ilr_objs,
     file = 'residuals_qned_counts_new.RData')

### End of Step 5 ###

### Step 6: Build the final counts table ###
print("Step 6.1 : Build the final counts table for clusters")

# Make subsets based on the different studies (Clusters)
all_rats = c()
ncols = 0
coolnames = c()
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('residuals_qned_counts',subset_id,sep='_'))
  all_rats = c(all_rats, colnames(this_residuals_qned_counts))
  ncols = ncols + dim(this_residuals_qned_counts)[1]
  coolnames = c(coolnames, paste(rownames(get(paste('residuals_qned_counts',subset_id,sep='_'))),subset_id,sep='_'))
}
all_rats = unique(all_rats)
count_data = matrix(ncol = ncols, nrow = length(all_rats))
colnames(count_data) = coolnames
rownames(count_data) = all_rats
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('residuals_qned_counts',subset_id,sep='_'))
  motch = match(all_rats, colnames(this_residuals_qned_counts))
  count_data[,paste(rownames(this_residuals_qned_counts),subset_id,sep='_')] = t(this_residuals_qned_counts[,motch])
  print(paste0("Subset:", subset_id))
}
motch = match(rownames(count_data), metadata[[sample_identifier]])
any(is.na(motch))
this_metadata = metadata[motch,]
rownames(count_data) = this_metadata[[sample_identifier]]
count_data[is.na(count_data)] = (-999)
assign(paste('count_data_clusters',sep='_'),count_data)

print("Step 6.2 : Build the final counts table for ranks")

# Make subsets based on the different studies (Ranks)
for (rank in ranks){
  all_rats = c()
  ncols = 0
  coolnames = c()
  for (subset_id in subset_ids) {
    this_residuals_qned_counts = get(paste('residuals_qned_counts',rank,subset_id,sep='_'))
    all_rats = c(all_rats, colnames(this_residuals_qned_counts))
    ncols = ncols + dim(this_residuals_qned_counts)[1]
    coolnames = c(coolnames, paste(rownames(get(paste('residuals_qned_counts',rank,subset_id,sep='_'))),subset_id,sep='_'))
  }
  all_rats = unique(all_rats)
  count_data = matrix(ncol = ncols, nrow = length(all_rats))
  colnames(count_data) = coolnames
  rownames(count_data) = all_rats
  for (subset_id in subset_ids) {
    this_residuals_qned_counts = get(paste('residuals_qned_counts',rank,subset_id,sep='_'))
    motch = match(all_rats, colnames(this_residuals_qned_counts))
    count_data[,paste(rownames(this_residuals_qned_counts),subset_id,sep='_')] = t(this_residuals_qned_counts[,motch])
    print(paste0(rank,":", subset_id))
  }
  motch = match(rownames(count_data), metadata[[sample_identifier]])
  any(is.na(motch))
  this_metadata = metadata[motch,]
  rownames(count_data) = this_metadata[[sample_identifier]]
  count_data[is.na(count_data)] = (-999)
  assign(paste('count_data',rank,sep='_'),count_data)
}

# Concatenate
count_data=t(get(paste("count_data",ranks[1],sep="_")))
if (length(ranks)>1){
  for (rank in ranks[2:length(ranks)]){
    count_rank=get(paste("count_data",rank,sep="_"))
    count_data=rbind(count_data,t(count_rank))
  }
}

count_data=t(count_data)

### Add cluster values to the final count data ###
print("Adding cluster values to the final count data")

count_data_alpha=count_data_clusters[rownames(count_data),]
count_data=cbind(count_data,count_data_clusters)

# Save
save(count_data,
     file = 'count_data_new.RData')

### End of Step 6 ###