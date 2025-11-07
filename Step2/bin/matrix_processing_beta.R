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

# 4. Inverse rank normalization and removal of covariate effects by regression
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
parser$add_argument("-cont", "--controls", type="character", help="File with input matrix")
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-lib", "--input_libsize", type="character", help="File with library depth info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the host",default="host_subject_id")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-phe", "--pheno", type="character", help="Define the phenotype: microbiome_taxonomic or microbiome_functional", default="microbiome_taxonomic")
## Taxonomic
parser$add_argument("-phy", "--phylotree", type="character", help="File with phylogeny tree info")
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-tax_ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## Functional
parser$add_argument("-func_ranks", "--functional_ranks", type="character", help="Functional ranks to be analysed",default="EC1,EC2,EC3,EC4")
## Beta Diversity
parser$add_argument("-divt", "--diversity_type", type="character", help="Alpha diversity level: 'PD' for 'Phylogenetic' and 'TD' for 'Taxonomic'",default="PD")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

## Modules
module_analysis<-args$module_analysis
module_comparisons<-args$module_comparisons
if (module_comparisons != "YES" & module_comparisons != "NO") {
  stop("Option for module alpha diversity should be either 'YES' or 'NO' in capital letters.", "\n")
}
## Matrix processing
controls<-args$controls
input_matrix<-args$input_matrix
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

## Beta Diversity
diversity<-args$diversity_type
allowed <- c("TD","PD")
string_vector <- strsplit(diversity, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All diversity types are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid diversity type(s) found: ", paste(invalid, collapse = ", "),"\n","Diversity types should be: 'TD' for Taxonomic Diversity and 'PD' for Phylogenetic Diversity." ,"\n")
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

# Read controls and input matrices
if (controls != "") {
  print("Reading controls")
  load(controls)
  if (module_comparisons=="YES"){
    matrix.sample_controls=matrix.sample_controls_objs[[1]]
    metadata_controls=metadata_controls_objs[[1]]
  }
} else {
  stop("No controls provided")
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

# Read diversity parameters
print("Processing diversity parameters")
div_types<-na.omit(c(sapply(strsplit(diversity, ","), "[",1),sapply(strsplit(diversity, ","), "[",2)))

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

### End of Step 1 ###

### Step 2: Beta diversity calculation ###
print("Step 2: Calculate beta diversity")

#if (method=="16S"){
#  for (subset_id in subset_ids){
#    # Collapse ASVs in GTDB Genomes (Species level)
#    print(paste0("Collapsing ASVs in GTDB Genomes (Species level): ",subset_id))
#    matrix = get(paste('filtered_prev_ASV',subset_id,sep='_'))
#    matrix = t(matrix)
#    colnames(matrix)=taxonomy$GTDB[match(colnames(matrix), taxonomy$Genome)]
#    collapsed_matrix=t(matrix)
#    collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
#    # Prune the phylogenetic tree and the matrix
#    print(paste0("Pruning the phylogenetic tree and the matrix: ",subset_id))
#    keep_taxa=rownames(collapsed_matrix)
#    remove_taxa = setdiff(tree$tip.label, keep_taxa)
#    pruned_tree = drop.tip(tree, remove_taxa)
#    pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
#    collapsed_matrix=collapsed_matrix[rownames(collapsed_matrix) %in% pruned_tree$tip.label,]
#    assign(paste('filtered_prev_Species',subset_id,sep='_'), collapsed_matrix)
#  }
#}


print("Calculating beta diversity per subset")
gp.ilr_objs=list()
gp.dist_objs=list()
beta_objs=list()
for (subset_id in subset_ids){
  for (i in c("with_controls", "without_controls")){
    if (pheno=="microbiome_taxonomic"){
      matrix = get(paste('filtered_prev_Species',subset_id,sep='_'))
    } else {
      matrix = get(paste('filtered_prev_EC4',subset_id,sep='_'))
    }
    if (i=="with_controls"){
      motch=match(rownames(matrix),rownames(matrix.sample_controls))
      this_controls=matrix.sample_controls[motch,]
      matrix=cbind(matrix,this_controls)
    }
    matrix=(round(matrix, digits = 0))
    if (pheno=="microbiome_taxonomic"){
      rownames(matrix)=taxonomy$Genome[match(rownames(matrix),taxonomy$Species)]
    }
    # Prune the phylogenetic tree
    keep_taxa=rownames(matrix)
    remove_taxa = setdiff(tree$tip.label, keep_taxa)
    pruned_tree = drop.tip(tree, remove_taxa)
    pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
    matrix=matrix[rownames(matrix) %in% pruned_tree$tip.label,]
    # Calculate dissimilarities in microbiome composition
    offset = 1
    filtered_trait_off= matrix + offset
    beta_total=data.frame()
    for (div_type in div_types){
      if (div_type=="PD"){
        gp.ilr <- philr(t(filtered_trait_off), pruned_tree, 
                          part.weights='enorm.x.gm.counts', 
                          ilr.weights='blw.sqrt')
        gp.dist <- dist(gp.ilr, method="euclidean")
      } else {
        gp.philr <- philr(t(filtered_genome_off), pruned_tree, 
                          part.weights='enorm.x.gm.counts', 
                          ilr.weights='uniform')
        gp.dist <- dist(gp.ilr, method="euclidean")
      }
      if (i=="without_controls"){
        # Calculate PCoA and extract the principal coordinates (PCs)
        genome.pcoa <- ape::pcoa(gp.dist)
        genome.pc <- data.frame(genome.pcoa$vectors)
        #Create matrix with dispersion values for host genetic effects analyses
        group <- rep("all", times = length(colnames(as.matrix(gp.dist))))
        betadisper_result <- betadisper(gp.dist,group, type="centroid")
        beta_values=as.data.frame(betadisper_result$distances)
        colnames(beta_values)="values"
        beta_values$samples=colnames(as.matrix(gp.dist))
        beta_values$values=as.numeric(beta_values$values)
        beta <- matrix(beta_values$values, nrow = 1, ncol = length(beta_values$samples), byrow = TRUE)
        colnames(beta)=beta_values$samples
        diversity=paste0("beta__",div_type,"_dispersion")
        rownames(beta)[1]=diversity
        PC1=t(as.matrix(genome.pc$Axis.1))
        colnames(PC1)=colnames(beta)
        rownames(PC1)[1]=paste0("beta__",div_type,"_PC1")
        PC2=t(as.matrix(genome.pc$Axis.2))
        colnames(PC2)=colnames(beta)
        rownames(PC2)[1]=paste0("beta__",div_type,"_PC2")
        beta=rbind(beta,PC1,PC2)
        beta_total=rbind(beta_total,beta)
        beta_obj <- paste('beta',subset_id,sep='_')
        assign(beta_obj, beta_total)
        beta_objs<- append(beta_objs, list(beta_obj=beta_total))
      }
      gp.ilr_obj <- paste('gp.ilr',i,div_type,subset_id,sep='_')
      assign(gp.ilr_obj, gp.ilr)
      gp.ilr_objs<- append(gp.ilr_objs, list(gp.ilr_obj=gp.ilr))
      gp.dist_obj <- paste('gp.dist',i,div_type,subset_id,sep='_')
      assign(gp.dist_obj, gp.dist)
      gp.dist_objs<- append(gp.dist_objs, list(gp.dist_obj=gp.dist))
      print(paste0(subset_id,"/",div_type,"/",i))
    }
  }
}

# Save Beta diversity objects for posterior analyses
save(gp.ilr_objs,
     gp.dist_objs,
     file="beta_diversity.RData")

### End of Step 2 ###

### Step 3: Plot PCoA for all covariates defined per subset ###
print("Step 3: Plot PCoA for all covariates defined per subset")
pdf('beta_diversity.pdf',bg='white')
for (subset_id in subset_ids){
  # Filter covariates and samples by subset
  print(subset_id)
  covariates=get(paste("covariates",subset_id,sep="_"))
  if (pheno=="microbiome_taxonomic"){
    matrix = get(paste('filtered_prev_Species',subset_id,sep='_'))
  } else {
    matrix = get(paste('filtered_prev_EC4',subset_id,sep='_'))
  }
  motch = match(colnames(matrix), metadata[[sample_identifier]])
  any(is.na(motch))
  this_metadata = metadata[motch,]
  # Add controls for Beta diversity analyses #
  motch=match(rownames(matrix),rownames(matrix.sample_controls))
  this_controls=matrix.sample_controls[motch,]
  matrix=cbind(matrix,this_controls)
  if (length(rownames(metadata_controls))>0){
    metadata_controls$Depth=0
    this_metadata=rbind(this_metadata,metadata_controls)
  }
  # Calculate PCoA
  for (covariate in covariates){
    for (div_type in div_types){
      if (div_type=="PD"){
        beta="Phylogenetic"
      } else {
        beta="Taxonomic"
      }
      # Calculate PCoA and extract the principal coordinates (PCs)
      genome.pcoa <- ape::pcoa(get(paste("gp.dist_with_controls",div_type,subset_id,sep="_")))
      genome.pc <- data.frame(genome.pcoa$vectors)
      # Make plot
      myplot<-ggplot(genome.pc, aes(x=Axis.1, y=Axis.2, colour = as.factor(this_metadata[[covariate]]))) +
        geom_point() + 
        ggtitle(paste0("Subset: ",subset_id," / Cov: ",covariate, " / Dist. Type: ",beta)) +
        xlab("PCo1") +
        ylab("PCo2") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 17),
              axis.line = element_line(),
              plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text = element_text(size = 13),
              axis.title = element_text(size = 16),
              legend.text=element_text(size=13),
              legend.title = element_blank(),
              legend.position="bottom")
      print(myplot)
    }
  }
}
dev.off()

### End of Step 3 ###

### Step 4: Calculate PERMANOVA ###
print("Step 4: Calculate PERMANOVA for all covariates defined per subset")
permanova_total=data.frame()
for (subset_id in subset_ids) {
  # Filter covariates and samples by subset
  covariates=get(paste("covariates",subset_id,sep="_"))
  covariates=c(covariates,"all")
  if (pheno=="microbiome_taxonomic"){
    matrix = get(paste('filtered_prev_Species',subset_id,sep='_'))
  } else {
    matrix = get(paste('filtered_prev_EC4',subset_id,sep='_'))
  }
  motch = match(colnames(matrix), metadata[[sample_identifier]]) #check colnames
  any(is.na(motch))
  this_metadata = metadata[motch,]
  for (div_type in div_types){
    gp.dist=get(paste("gp.dist_without_controls",div_type,subset_id,sep="_"))
    for (covariate in covariates){
      if (covariate=="all"){
        formula <- as.formula(paste("gp.dist ~", paste(covariates[-length(covariates)], collapse = " + ")))
      } else {
        formula <- as.formula(paste("gp.dist ~", paste(covariate, collapse = " + ")))
      }
      permanova.beta <- adonis2(formula, data = this_metadata, permutations = 999)
      permanova.beta$subset=subset_id
      permanova.beta$covariate=covariate
      permanova.beta$diversity=div_type
      permanova_total=rbind(permanova_total,permanova.beta)
      print(paste0(subset_id,"/",covariate,"/",div_type))
    }
  }
}

# Save PERMANOVA results for posterior analyses
save(permanova_total,
     file="permanova_beta_diversity.RData")

### End of Step 4 ###

### Step 5: Rank normalization and removal of known covariate effects by regression ###
print("Step 5: Rank normalization and removal of known covariates effects by regression")
residuals_qned_counts_beta_objs<-list()
gp.dist_objs=list()
beta_tot=data.frame()
for (subset_id in subset_ids) {
  for (div_type in div_types){
    #matrix = get(paste("beta",subset_id,sep='_'))
    matrix = as.matrix(get(paste('gp.ilr_without_controls',div_type,subset_id,sep='_')))
    colnames(matrix) <- paste0("n", seq_len(ncol(matrix)))
    matrix=t(matrix)
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
    Study_cov = Study_cov[cc]
    this_matrix=matrix
    names=rownames(this_matrix)
    this_matrix = this_matrix[,cc]
    rownames(this_matrix)=names
    qned_counts = t(apply(this_matrix, FUN = invrank, MAR = 1))
    residuals_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1))
    colnames(residuals_qned_counts) = colnames(this_matrix)
    # Calculate distances
    gp.dist <- dist(t(residuals_qned_counts), method="euclidean")
    # Calculate PCoA and extract the principal coordinates (PCs)
    genome.pcoa <- ape::pcoa(gp.dist)
    genome.pc <- data.frame(genome.pcoa$vectors)
    #Create matrix with dispersion values for host genetic effects analyses
    group <- rep("all", times = length(colnames(as.matrix(gp.dist))))
    betadisper_result <- betadisper(gp.dist,group, type="centroid")
    beta_values=as.data.frame(betadisper_result$distances)
    colnames(beta_values)="values"
    beta_values$samples=colnames(as.matrix(gp.dist))
    beta_values$values=as.numeric(beta_values$values)
    beta <- matrix(beta_values$values, nrow = 1, ncol = length(beta_values$samples), byrow = TRUE)
    colnames(beta)=beta_values$samples
    diversity=paste0("beta__",div_type,"_dispersion")
    rownames(beta)[1]=diversity
    PC1=t(as.matrix(genome.pc$Axis.1))
    colnames(PC1)=colnames(beta)
    rownames(PC1)[1]=paste0("beta__",div_type,"_PC1")
    PC2=t(as.matrix(genome.pc$Axis.2))
    colnames(PC2)=colnames(beta)
    rownames(PC2)[1]=paste0("beta__",div_type,"_PC2")
    beta=rbind(beta,PC1,PC2)
    beta_tot=rbind(beta_tot,beta)
    #Define objects to store the data
    gp.dist_obj <- paste('gp.dist',div_type,subset_id,sep='_')
    assign(gp.dist_obj, gp.dist)
    gp.dist_objs<- append(gp.dist_objs, list(gp.dist_obj=gp.dist))
  }
  residuals_qned_counts_beta_obj <- paste('beta',subset_id,sep='_')
  assign(residuals_qned_counts_beta_obj, beta_tot)
  residuals_qned_counts_beta_objs<- append(residuals_qned_counts_beta_objs, list(residuals_qned_counts_beta_obj=beta_tot))
  print(subset_id)
}

save(gp.dist_objs,
     residuals_qned_counts_beta_objs,
     file = 'residuals_qned_counts_beta.RData')

### End of Step 5 ###

### Step 6: Build the final counts table ###
print("Step 6: : Build the final counts table")

all_rats = c()
ncols = 0
coolnames = c()
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('beta',subset_id,sep='_'))
  all_rats = c(all_rats, colnames(this_residuals_qned_counts))
  ncols = ncols + dim(this_residuals_qned_counts)[1]
  coolnames = c(coolnames, paste(rownames(get(paste('beta',subset_id,sep='_'))),subset_id,sep='_'))
}
all_rats = unique(all_rats)
count_data = matrix(ncol = ncols, nrow = length(all_rats))
colnames(count_data) = coolnames
rownames(count_data) = all_rats
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('beta',subset_id,sep='_'))
  motch = match(all_rats, colnames(this_residuals_qned_counts))
  count_data[,paste(rownames(this_residuals_qned_counts),subset_id,sep='_')] = t(this_residuals_qned_counts[,motch])
  print(subset_id)
}
motch = match(rownames(count_data), metadata$host_subject_id)
this_metadata = metadata[motch,]
rownames(count_data) = this_metadata$host_subject_id
count_data[is.na(count_data)] = (-999)
assign('count_data_beta',count_data)

# Save
save(count_data_beta,
     file = 'count_data_beta.RData')

### End of Step 6 ###