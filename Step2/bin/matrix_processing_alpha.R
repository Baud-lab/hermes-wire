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
## Alpha Diversity
parser$add_argument("-divt", "--diversity_type", type="character", help="Alpha diversity level: 'PD' for 'Phylogenetic' and 'TD' for 'Taxonomic'",default="PD")
parser$add_argument("-alpha_h", "--alpha_hill_number", type="character", help="Types of diversity measurements: '0' for 'Richness', '1' for 'Evenness' and '2' for 'Probability'",default="2")
parser$add_argument("-alpha_n", "--alpha_nboot", type="double", help="Number of bootstrap replications for the estimation of the asymptotic alpha diversity",default=50)
parser$add_argument("-alpha_c", "--alpha_conf", type="double", help="Level of confidence interval for the estimation of the asymptotic alpha diversity",default=0.95)

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

## Alpha Diversity
diversity<-args$diversity_type
allowed <- c("TD","PD")
string_vector <- strsplit(diversity, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All diversity types are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid diversity type(s) found: ", paste(invalid, collapse = ", "),"\n","Diversity types should be: 'TD' for Taxonomic Diversity and 'PD' for Phylogenetic Diversity." ,"\n")
}
alpha_h<-args$alpha_hill_number
allowed <- c("0","1","2")
string_vector <- strsplit(alpha_h, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All Hill numbers are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid Hill number(s) found: ", paste(invalid, collapse = ", "),"\n","Hill number should be: '0' for 'Richness', '1' for 'Evenness' and '2' for 'Probability'." ,"\n")
}
alpha_n<-as.numeric(args$alpha_nboot)
alpha_c<-as.numeric(args$alpha_conf)

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
hill_numbers<-na.omit(c(as.numeric(sapply(strsplit(alpha_h, ","),"[",1)),
                        as.numeric(sapply(strsplit(alpha_h, ","), "[",2)),
                        as.numeric(sapply(strsplit(alpha_h, ","), "[",3))))

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

### Step 2: Alpha diversity) ###
print("Step 2: Calculate alpha diversity")

#subset_ids="ALL" ####### Remove!!!!

#if (method=="16S"){
#  for (subset_id in subset_ids){
#    # Collapse ASVs in GTDB Genomes (Species level)
#    print(paste0("Collapsing ASVs in GTDB Genomes (Species level): ",subset_id))
#    matrix = get(paste('filtered_prev_Species',subset_id,sep='_'))
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

# Calculate the asymptotic diversity for all samples according to the subset
print("Calculating the asymptotic diversity for all samples according to the subset")

# Select the number of cores
cl <- makeCluster(detectCores())
print(cl)
registerDoParallel(cl)

# Calculate Alpha Diversity in parallel
alpha_table_objs=list()
for (subset_id in subset_ids){
  covariates=get(paste("covariates",subset_id,sep="_"))
  if (pheno=="microbiome_taxonomic"){
    matrix = get(paste('filtered_prev_Species',subset_id,sep='_'))
  } else {
    matrix = get(paste('filtered_prev_EC4',subset_id,sep='_'))
  }
  motch = match(colnames(matrix), metadata[[sample_identifier]])
  this_metadata = metadata[motch,]
  matrix=round(matrix, digits = 0)
  if (pheno=="microbiome_taxonomic"){
    rownames(matrix)=taxonomy$Genome[match(rownames(matrix),taxonomy$Species)]
  }
  labels=rownames(matrix)
  remove_taxa = setdiff(tree$tip.label, labels)
  pruned_tree = drop.tip(tree, remove_taxa)
  pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
  for (div_type in div_types){
    for (hill_number in hill_numbers){
      alpha_table <- foreach(i = 1:ncol(matrix), .combine = rbind) %dopar% {
        col_data <- matrix[,i]
        names(col_data)=rownames(matrix)
        iNEXT.3D::AO3D(data=col_data,
                       diversity = div_type,
                       q = as.numeric(hill_number),
                       datatype = "abundance",
                       nboot = alpha_n,
                       conf = alpha_c,
                       method = "Asymptotic",
                       PDtree = pruned_tree)
      }
      alpha_table$Assemblage=colnames(matrix)
      colnames(alpha_table)[3]="qD"
      for (covariate in covariates){
        alpha_table[[covariate]]=this_metadata[[covariate]][match(alpha_table$Assemblage,this_metadata[[sample_identifier]])]
      }
      alpha_table_obj <- paste('alpha_table',div_type,hill_number,subset_id,sep='_')
      assign(alpha_table_obj, alpha_table)
      alpha_table_objs<- append(alpha_table_objs, list(alpha_table_obj=alpha_table))
      print(paste0(subset_id,"/",div_type,"/q_",hill_number))
    }
  }
}

# Stop parallel back end and save matrix
stopImplicitCluster()
save(alpha_table_objs,
     file="alpha_diversity.RData")

### End of Step 2 ###

### Step 3: Plot boxplots for all covariates defined per subset ###
print("Step 3: Plot boxplots for all covariates defined per subset")
pdf('alpha_diversity.pdf',bg='white')
for (subset_id in subset_ids){
  covariates=get(paste("covariates",subset_id,sep="_"))
  for (div_type in div_types){
    for (hill_number in hill_numbers){
      if (div_type=="TD"){
        if (hill_number=="0"){
          axis_text="Richness"
        } else {
          if (hill_number=="1"){
            axis_text="Shannon Entropy"
          } else {
            axis_text="Simpson's Index"
          }
        }
      } else {
        if (hill_number=="0"){
          axis_text="Faith's Phylogenetic Diversity"
        } else {
          if (hill_number=="1"){
            axis_text="Phylogenetic Entropy"
          } else {
            axis_text="Rao's Quadratic Index"
          }
        }
      }
      # Define limits for barplots
      alpha_table <- get(paste('alpha_table',div_type,hill_number,subset_id,sep='_'))
      # Make plot
      for (covariate in covariates){
        myplot<-ggplot(alpha_table, aes(x=reorder(!!as.name(covariate),qD,FUN = median), y=qD)) +
          geom_boxplot() +
          xlab(paste0("Subset: ",subset_id," / Covariate: ",covariate)) +
          #ylab(paste0("Asymptotic Diversity ","(",div_type,": q=",hill_number,")")) +
          ylab(paste0("Asymptotic ",axis_text)) +
          theme_bw() +
          theme(legend.position="bottom",
                axis.line = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.text=element_text(size=10))+
          scale_x_discrete(guide = guide_axis(angle = 45)) +
          geom_signif(comparisons = list(unique(alpha_table[[covariate]])),
                      map_signif_level = TRUE)
        print(myplot)
      }
    }
  }
}
dev.off()

### End of Step 3 ###

### Step 4: Calculate ANOVA ###
print("Step 4: Calculate ANOVA for all covariates defined per subset. NOTE: Inverse rank normalization will be applied")
anova_total=data.frame()
for (subset_id in subset_ids){
  covariates=get(paste("covariates",subset_id,sep="_"))
  covariates=c(covariates,"all")
  for (div_type in div_types){
    for (hill_number in hill_numbers){
      alpha_table <- get(paste('alpha_table',div_type,hill_number,subset_id,sep='_'))
      #Transform data to make it parametric
      alpha_table$norm_qD = invrank(alpha_table$qD)
      for (covariate in covariates){
        if (covariate=="all"){
          formula <- as.formula(paste("norm_qD ~", paste(covariates[-length(covariates)], collapse = " + ")))
        } else {
          formula <- as.formula(paste("norm_qD ~", paste(covariate, collapse = " + ")))
        }
        anova.alpha <- summary(aov(formula, data = alpha_table))
        anova.alpha=anova.alpha[[1]]
        anova.alpha$R2=anova.alpha$`Sum Sq`/sum(anova.alpha$`Sum Sq`)
        anova.alpha$subset=subset_id
        anova.alpha$covariate=covariate
        anova.alpha$diversity=div_type
        anova.alpha$q=hill_number
        anova_total=rbind(anova_total,anova.alpha)
        print(paste0(subset_id,"/",covariate,"/",div_type, "/",hill_number))
      }
    }
  }
}

# Save ANOVA results for posterior analyses
save(anova_total,
     file="anova_alpha_diversity.RData")

### End of Step 4 ###

### Step 5: Create matrix with alpha diversity values for host genetic effects analyses ###
print("Step 5: Create matrix with alpha diversity values for host genetic effects analyses")
alpha_objs=c()
for (subset_id in subset_ids){
  alpha_total=data.frame()
  for (div_type in div_types){
    for (hill_number in hill_numbers){
      alpha_table <- get(paste('alpha_table',div_type,hill_number,subset_id,sep='_'))
      alpha_table=alpha_table[c(3,1)]
      alpha <- matrix(alpha_table$qD, nrow = 1, ncol = length(alpha_table$Assemblage), byrow = TRUE)
      colnames(alpha)=alpha_table$Assemblage
      diversity=paste0("alpha__",div_type,"_q",hill_number)
      rownames(alpha)[1]=diversity
      alpha_total=rbind(alpha_total,alpha)
      alpha_obj <- paste('alpha',subset_id,sep='_')
      assign(alpha_obj, alpha_total)
      alpha_objs<- append(alpha_objs, list(alpha_obj=alpha_total))
      print(paste0(subset_id,"/",div_type, "/",hill_number))
    }
  }
}

### End of Step 5 ###

### Step 6: Rank normalization and removal of known covariate effects by regression ###
print("Step 6: Rank normalization and removal of known covariates effects by regression")
residuals_qned_counts_alpha_objs<-list()
for (subset_id in subset_ids) {
  matrix = get(paste("alpha",subset_id,sep='_'))
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
  residuals_qned_counts_alpha_obj <- paste('residuals_qned_counts',"alpha",subset_id,sep='_')
  assign(residuals_qned_counts_alpha_obj, residuals_qned_counts)
  #Define objects to store the data
  residuals_qned_counts_alpha_objs<- append(residuals_qned_counts_alpha_objs, list(residuals_qned_counts_alpha_obj=residuals_qned_counts))
  print(subset_id)
}

save(residuals_qned_counts_alpha_objs,
     file = 'residuals_qned_counts_alpha.RData')

### End of Step 6 ###

### Step 7: Build the final counts table ###
print("Step 7: : Build the final counts table")

all_rats = c()
ncols = 0
coolnames = c()
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('residuals_qned_counts_alpha',subset_id,sep='_'))
  all_rats = c(all_rats, colnames(this_residuals_qned_counts))
  ncols = ncols + dim(this_residuals_qned_counts)[1]
  coolnames = c(coolnames, paste(rownames(get(paste('residuals_qned_counts_alpha',subset_id,sep='_'))),subset_id,sep='_'))
}
all_rats = unique(all_rats)
count_data = matrix(ncol = ncols, nrow = length(all_rats))
colnames(count_data) = coolnames
rownames(count_data) = all_rats
for (subset_id in subset_ids) {
  this_residuals_qned_counts = get(paste('residuals_qned_counts_alpha',subset_id,sep='_'))
  motch = match(all_rats, colnames(this_residuals_qned_counts))
  count_data[,paste(rownames(this_residuals_qned_counts),subset_id,sep='_')] = t(this_residuals_qned_counts[,motch])
  print(subset_id)
}
motch = match(rownames(count_data), metadata$host_subject_id)
this_metadata = metadata[motch,]
rownames(count_data) = this_metadata$host_subject_id
count_data[is.na(count_data)] = (-999)
assign('count_data_alpha',count_data)

# Save
save(count_data_alpha,
     file = 'count_data_alpha.RData')

### End of Step 7 ###