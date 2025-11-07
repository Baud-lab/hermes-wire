#!/usr/bin/env Rscript

# Packages
suppressMessages(library(qvalue))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(gap))
suppressMessages(library(argparse))
suppressMessages(library(xfun))
suppressMessages(library(ggsignif))
suppressMessages(library(ape))
suppressMessages(library(tidyr))

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

# 3. Inverse rank normalization
invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}

# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define arguments:
parser$add_argument("-mod_ana", "--module_analysis", type="character", help="Inclusion of module: Analysis",default="Pooled")
parser$add_argument("-mod_comp", "--module_comparisons", type="character", help="Inclusion of module: Comparisons",default="NO")
parser$add_argument("-mod_clust", "--module_cluster", type="character", help="Inclusion of module: Clusters",default="NO")
parser$add_argument("-inputN", "--input_null", type="character", help="Input Null file")
parser$add_argument("-inputF", "--input_full", type="character", help="Input Full file")
parser$add_argument("-kin", "--kinship_version", type="character", help="Name of the kinship model",default="exp_DGA_7_Dec_2021")
parser$add_argument("-model", "--model_type", type="character", help="Name of the model",default="uni")
parser$add_argument("-pheno", "--phenotype_version", type="character", help="Phenotype version",default="bacterial_taxa")
parser$add_argument("-effN", "--effect_null", type="character", help="Null effect to be checked")
parser$add_argument("-effF", "--effect_full", type="character", help="Full effect to be checked")
parser$add_argument("-sig_f", "--significance_fdr", type="double", help="Significance level (FDR)",default=0.05)
parser$add_argument("-sig_b", "--significance_bonferroni", type="double", help="Significance level Bonferroni)",default=0.05)
parser$add_argument("-mp", "--min_prevalence", type="double", help="Minimum prevalence accepted", default=0.5)
parser$add_argument("-tax_ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
parser$add_argument("-func_ranks", "--functional_ranks", type="character", help="Functional ranks to be analysed",default="EC1,EC2,EC3,EC4")
parser$add_argument("-prev", "--filtered_prev", type="character", help="R object with all matrices containing read counts per taxonomic rank and study after filtering taxa by prevalence")
parser$add_argument("-stats", "--catalogue_stats", type="character", help="File with statistics about GTDB Reference Genomes")
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-phy", "--phylotree", type="character", help="File with phylogeny tree info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-mets", "--methods", type="character", help="Define the methods to be compared", default="Shallow,16S")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()


## Arguments
module_analysis<-args$module_analysis
module_comparisons<-args$module_comparisons
module_cluster<-args$module_cluster
MODEL=args$model_type
GRM_V=args$kinship_version
INPUTF=args$input_full
INPUTN=args$input_null
EFF_full=args$effect_full
EFF_null=args$effect_null
sig_f=args$significance_fdr
sig_b=args$significance_bonferroni
min_prevalence<-args$min_prevalence
tax_ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(tax_ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}
func_ranks<-args$functional_ranks
allowed <- c("EC1", "EC2", "EC3", "EC4")
string_vector <- strsplit(func_ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be:EC1,EC2,EC3 or EC4" ,"\n")
}
input_prev=args$filtered_prev
catalogue_stats=args$catalogue_stats
PHENO_V=args$phenotype_version
OUT_DIR="./"
input_taxonomy<-args$input_taxonomy
phylotree<-args$phylotree
list_of_covariates<-args$list_of_covariates
method=args$method
methods<-args$methods
pheno<-args$pheno


# -----------------------------------------------------------------------------

### Step 1 ####
print("Step 1: Read input files and arguments")

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
  taxonomy$ASV=paste0("asv__",rownames(taxonomy))
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
if (pheno=="microbiome_taxonomic"){
  print("Reading phylogenetic info")
  if (phylotree == ""){
    stop("Please, provide the phylogeny if you want to estimate the phylogenetic diversity (PD)")
  } else {
    tree = ape::read.tree(phylotree)
  } 
}

# Load model "Full"
load(INPUTF)
all_VCs_full = all_VCs
dim(all_VCs_full)
## Load model "Null"
load(INPUTN)
all_VCs_null = all_VCs
dim(all_VCs_null)

# Read methods
print("Reading methods")
methods<-na.omit(c((sapply(strsplit(methods, ","),"[",1)),
                   (sapply(strsplit(methods, ","), "[",2))))
if (methods[1]!=method){
  stop("ERROR: The 1st method pointed on the 'methods' field should be the main method to be used for profiling.")
}

# Read count tables with taxa filtered by prevalence
if (input_prev != "") {
  print("Reading filtered matrices per subset")
  load(input_prev)
  if (module_comparisons=="YES" & module_cluster!="YES"){
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

# Read catalogue stats
if (pheno=="microbiome_taxonomic"){
  if (catalogue_stats ==""){
    stop("Please, provide a file with GTDB Reference Genomes statistics.")
  } else {
    print("Reading genome catalogue statistics")
    catalogue_stats = read_csv_file(catalogue_stats)
    if (method =="16S"){
      catalogue_stats$ASV=taxonomy$Genome[match(catalogue_stats$Genome,taxonomy$GTDB)]
    } else {
      catalogue_stats$Species=taxonomy$Species[match(rownames(catalogue_stats),taxonomy$Genome)]
    }
  }
}

### End of Step 1 ###

### Step 2 ####
print("Step 2: Build heritability table")

inter=intersect(rownames(all_VCs_null),rownames(all_VCs_full))
length(inter)
uniq_taxa = unique(c(all_VCs_full$trait, all_VCs_null$trait))
all_VCs_full = all_VCs_full[match(uniq_taxa, all_VCs_full$trait),]
all_VCs_null = all_VCs_null[match(uniq_taxa, all_VCs_null$trait),]
all_VCs_full$pvalue_DGE = pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=1,lower.tail = F)
all_VCs_full = all_VCs_full[!is.na(all_VCs_full$pvalue_DGE),]
all_VCs_full$logP=as.numeric(-log10(all_VCs_full$pvalue_DGE))
all_VCs_full$subset_id <- sapply(strsplit(all_VCs_full$trait, "_"), tail, 1)
#all_VCs_full$subset_id[all_VCs_full$subset_id == "behavior" ] <- "TN_behavior" #change
#all_VCs_full$subset_id[all_VCs_full$subset_id == "breeder" ] <- "TN_breeder" #change
for (subset_id in subset_ids){
  all_VCs_full$trait=gsub(paste("\\_",subset_id,"$",sep=""),"",all_VCs_full$trait)
}
all_VCs_full$rank <- sapply(strsplit(all_VCs_full$trait, "_"), head, 1)
for (i in 1:nrow(all_VCs_full)){
  if (all_VCs_full$rank[i]=="p"){
    all_VCs_full$rank[i]="Phylum"
  } else{
    if (all_VCs_full$rank[i]=="c"){
      all_VCs_full$rank[i]="Class"
    } else{
      if (all_VCs_full$rank[i]=="o"){
        all_VCs_full$rank[i]="Order"
      } else{
        if (all_VCs_full$rank[i]=="f"){
          all_VCs_full$rank[i]="Family"
        } else{
          if (all_VCs_full$rank[i]=="g"){
            all_VCs_full$rank[i]="Genus"
          } else{
            if (all_VCs_full$rank[i]=="s"){
              all_VCs_full$rank[i]="Species"
            } else{
              if (all_VCs_full$rank[i]=="alpha"){
                all_VCs_full$rank[i]="Alpha Diversity"
              } else{
                if (all_VCs_full$rank[i]=="beta"){
                  all_VCs_full$rank[i]="Beta Diversity"
                } else {
                  if (all_VCs_full$rank[i]=="asv"){
                    all_VCs_full$rank[i]="ASV"
                  } else {
                    if (all_VCs_full$rank[i]=="cluster"){
                      all_VCs_full$rank[i]="Clusters"
                    } else {
                      if (all_VCs_full$rank[i]=="enterotype") {
                        all_VCs_full$rank[i]="Enterotypes"
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

# Bonferroni and FDR correction
print("Processing Bonferroni and FDR correction")

#subset_id="ALL" ### Remove
#all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove

all_VCs_full_subsets=all_VCs_full
if (module_analysis=="Pooled"){
  all_VCs_full=all_VCs_full[all_VCs_full$subset_id=="ALL",]
} else {
  all_VCs_full=all_VCs_full[all_VCs_full$subset_id!="ALL",]
}

all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
all_VCs_full <- all_VCs_full %>% 
  mutate( Significance = case_when(
    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Significant (Bonferroni)",
    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Significant (FDR)",
    TRUE ~ "Non-significant traits"))
all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

# Save results
save(all_VCs_full,
     all_VCs_full_subsets,
     file ="heritability.RData")

### End of Step 2 ###

### Step 3 ####
print("Step 3: Make plots")

# Plots
## Define maximum values to be used for the axes of all plots
max_logP = max(-log10(all_VCs_full[,'pvalue_DGE']), na.rm = T)
max_var_Ad = max(100*(all_VCs_full[,'var_Ad']), na.rm = T)

## 1. Q-Q plots
print("Q-Q plots")
pdf('heritability_QQplot.pdf',bg='white', width = 12, height = 4)
par(mfrow = c(1,length(subset_ids)), mar = c(2,2,3,1))
for (subset_id in subset_ids){
  all_VCs_plot=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  r=qqunif(all_VCs_plot[,'pvalue_DGE'],ci=T,col = 'black',las=1, main = paste0("Subset = ",subset_id), ylim = c(0,max_logP))
}
dev.off()

## 2. Heritability vs -LogP
print("Heritability vs -LogP plots")
pdf('heritability_vs_logp.pdf',bg='white')
color_groups <- c("Non-significant traits", "Significant (FDR)", "Significant (Bonferroni)")
color_labels <- c("black", "firebrick3", "blue")
for (subset_id in subset_ids){
  all_VCs_plot=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  myplot<-ggplot(all_VCs_plot, aes(x=var_Ad*100, y=-log10(pvalue_DGE))) + 
    geom_point(aes(color = Significance), size = 4/5) +
    xlab(paste0("Heritability (%) ","(Subset = ",subset_id,")")) +
    ylab(expression("-log"[10]*"P-Value Heritability")) +
    theme_bw() +
    theme(legend.position="bottom",
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13)) +
    coord_cartesian(xlim = c(0, max_var_Ad), ylim = c(0, max_logP)) +
    scale_y_continuous(breaks = seq(0, max_logP, by = 2)) +
    scale_x_continuous(breaks = seq(0, max_var_Ad, by = 5)) +
    scale_color_manual(values = color_labels,
                       labels = color_groups,
                       limits = color_groups) +
    labs(color = "")
  print(myplot)
}
dev.off()

## 3. Boxplots
print("Boxplots")
pdf('heritability_boxplot.pdf',bg='white')

if (module_analysis!="Pooled"){
  #Comparisons by subset
  all_VCs_plot=all_VCs_full[!grepl("ALL",all_VCs_full$subset_id),]
  myplot<-ggplot(all_VCs_plot, aes(x=reorder(subset_id,var_Ad*100,FUN = median), y=var_Ad*100)) +
    geom_boxplot() +
    xlab("") +
    ylab("Heritability (%)") +
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
    geom_signif(comparisons = list(unique(all_VCs_plot$subset_id)),
                map_signif_level = TRUE)
  print(myplot)
}

#Comparisons by taxonomic rank
order=c()
if ("Alpha Diversity" %in% unique(all_VCs_full$rank)){
  order=c(order,"Alpha Diversity")
}
if ("Beta Diversity" %in% unique(all_VCs_full$rank)){
  order=c(order,"Beta Diversity")
}
if ("Enterotypes" %in% unique(all_VCs_full$rank)){
  order=c(order,"Enterotypes")
}
if ("Clusters" %in% unique(all_VCs_full$rank)){
  order=c(order,"Clusters")
}
order=c(order,ranks)
for (subset_id in subset_ids){
  all_VCs_plot=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  myplot<-ggplot(all_VCs_plot, aes(x=factor(rank,level = order), y=var_Ad*100)) +
    geom_boxplot() +
    xlab("") +
    ylab(paste0("Heritability (%) ","(Subset = ",subset_id,")")) +
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
    scale_x_discrete(guide = guide_axis(angle = 45))
    #geom_signif(comparisons = list(c("Genus","Species")),
    #            map_signif_level = TRUE)
  print(myplot)
}
dev.off()

if (pheno=="microbiome_taxonomic"){
  if (method!="16S"){
    ## 4. Phylogenetic trees
    if (phylotree != "" && "Species" %in% ranks){
      print("Phylogenetic trees")
      all_VCs_species=all_VCs_full[grepl("Species",all_VCs_full$rank),]
      all_VCs_species$genome=taxonomy$Genome[match(all_VCs_species$trait,taxonomy$Species)]
      pdf('heritability_phylotrees.pdf',bg='white', width = 12, height = 4)
      par(mfrow = c(1,3), mar = c(4,2,2,2))
      for (subset_id in subset_ids){
        all_VCs_plot=all_VCs_species[grepl(subset_id,all_VCs_species$subset_id),]
        keep_taxa=all_VCs_plot$genome
        any(is.na(keep_taxa))
        remove_taxa = setdiff(tree$tip.label, keep_taxa)
        any(is.na(remove_taxa))
        pruned_tree = drop.tip(tree, remove_taxa)
        egde_labels <- all_VCs_plot$Significance
        pruned_tree$egde.label <- egde_labels
        pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
        order=all_VCs_plot[match(pruned_tree$tip.label,all_VCs_plot$genome),]
        significance_colors <- ifelse(order$Significance == "Significant (FDR)", "red",
                                      ifelse(order$Significance == "Significant (Bonferroni)", "blue", "lightgrey"))
        plot.phylo(pruned_tree, cex = 0.8, type = "fan", edge.color = significance_colors, show.tip.label = FALSE, show.node.label = FALSE, main = paste0("Subset: ",subset_id))
        edge_labels=c("FDR<10%","Bonf<5%")
        significance_colors=c("red","blue")
        legend("bottom", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
      }
      dev.off()
    }
  }
}


## 5. Abundance vs prevalence
print("Abundance vs prevalence")
edge_labels=c(paste0("FDR<",100*sig_f,"%"),paste0("Bonf<",100*sig_b,"%"))
significance_colors=c("red","blue")
pdf('heritability_abundance_vs_prevalence.pdf',bg='white', width = 12, height = 7)
par(mar = c(5,5,3,1))
for (subset_id in subset_ids){
  for (rank in ranks){
    if (module_comparisons=="YES") {
      matrix = get(paste('filtered_prev',rank,subset_id,method,sep='_'))
    } else {
      matrix = get(paste('filtered_prev',rank,subset_id,sep='_'))
    }
    ab_vs_prev=data.frame(trait=rownames(matrix),
                          abundance=100*apply(sweep(matrix,2,colSums(matrix),"/"), 1, mean, na.rm=TRUE),
                          prevalence=100*rowSums(matrix>0)/ncol(matrix))
    all_VCs=all_VCs_full[grepl(rank,all_VCs_full$rank) & grepl(subset_id,all_VCs_full$subset_id),]
    motch=match(all_VCs$trait,ab_vs_prev$trait)
    ab_vs_prev=ab_vs_prev[motch,]
    ab_vs_prev$Significance <- all_VCs$Significance[match(ab_vs_prev$trait,all_VCs$trait)]
    if (rank == "Genus" || rank == "Species"){
      max_ab=20
    } else {
      if (rank == "Family" || rank == "Order"){
        max_ab=50
      } else {
        max_ab=100
      }
    }
    plot(x=1, type = "n", ylab=paste0("Mean Abundance (%) (Subset: ",subset_id," / Rank: ",rank,")"), xlab=paste0("Prevalence (%) (Subset: ",subset_id," / Rank: ",rank,")"), ylim=c(0,max_ab), xlim=c(100*min_prevalence,100))
    points(x = ab_vs_prev$prevalence[ab_vs_prev$Significance == "Non-significant traits"],
           y = ab_vs_prev$abundance[ab_vs_prev$Significance == "Non-significant traits"],
           pch = 19,
           col = "grey")
    points(x = ab_vs_prev$prevalence[ab_vs_prev$Significance == "Significant (FDR)"],
           y = ab_vs_prev$abundance[ab_vs_prev$Significance == "Significant (FDR)"],
           pch = 19,
           col = "red")
    points(x = ab_vs_prev$prevalence[ab_vs_prev$Significance == "Significant (Bonferroni)"],
           y = ab_vs_prev$abundance[ab_vs_prev$Significance == "Significant (Bonferroni)"],
           pch = 19,
           col = "blue")
    legend("top", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
    assign(paste('ab_vs_prev',rank,subset_id,sep='_'), ab_vs_prev)
    print(paste0(subset_id,"/",rank))
  }
}
dev.off()

# 6. Overall statistics
print("Overall statistics")
pdf('heritability_statistics.pdf',bg='white', width = 8, height = 6)
par(mar = c(2,2,4,2))

## 1. Per subsets
statistics=data.frame(subset=subset_ids)
### Samples
statistics$samples=all_VCs_full$sample_size[match(statistics[[1]],all_VCs_full$subset_id)]
### Traits
traits=count(all_VCs_full,subset_id)
statistics$traits=traits$n[match(statistics[[1]],traits$subset_id)]
### Significant (FDR)
sig_fdr=all_VCs_full[!all_VCs_full$Significance=="Non-significant traits",]
stats_sig_herit=sig_fdr %>%                            
  group_by(subset_id) %>%
  summarize(min = min(var_Ad*100), median=median(var_Ad*100), max = max(var_Ad*100))
sig_fdr=count(sig_fdr, subset_id)
statistics$sig_fdr=sig_fdr$n[match(statistics[[1]],sig_fdr$subset_id)]
statistics$percentage_fdr=100*statistics$sig_fdr/statistics$traits
statistics$fdr_min_herit=stats_sig_herit$min[match(statistics[[1]],stats_sig_herit$subset_id)]
statistics$fdr_median_herit=stats_sig_herit$median[match(statistics[[1]],stats_sig_herit$subset_id)]
statistics$fdr_max_herit=stats_sig_herit$max[match(statistics[[1]],stats_sig_herit$subset_id)]
### Significant (Bonferroni)
sig_bonf=all_VCs_full[all_VCs_full$Significance=="Significant (Bonferroni)",]
stats_sig_herit=sig_bonf %>%                            
  group_by(subset_id) %>%
  summarize(min = min(var_Ad*100), median=median(var_Ad*100), max = max(var_Ad*100))
sig_bonf=count(sig_bonf, subset_id)
statistics$sig_bonf=sig_bonf$n[match(statistics[[1]],sig_bonf$subset_id)]
statistics$percentage_bonf=100*statistics$sig_bonf/statistics$traits
statistics$bonf_min_herit=stats_sig_herit$min[match(statistics[[1]],stats_sig_herit$subset_id)]
statistics$bonf_median_herit=stats_sig_herit$median[match(statistics[[1]],stats_sig_herit$subset_id)]
statistics$bonf_max_herit=stats_sig_herit$max[match(statistics[[1]],stats_sig_herit$subset_id)]
statistics <- replace(statistics, is.na(statistics), 0)
assign("statistics_subsets",statistics)
### Make bar plot
statistics_plot <- tidyr::pivot_longer(statistics_subsets, cols = c(traits, sig_fdr, sig_bonf), names_to = "Traits", values_to = "value")
statistics_plot$value[is.na(statistics_plot$value)] <- 0
max_value=max(statistics_plot$value)
myplot <- ggplot(statistics_plot, aes(x = subset, y = value, fill = Traits)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, max_value+100)) +
  labs(x = "", y = "Number of traits") +
  scale_fill_manual(values = c("traits" = "azure4", "sig_fdr" = "red", "sig_bonf" = "blue"),
                    labels = c("Significant (Bonferroni)", "Significant (FDR)", "All traits"),
                    name = "") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5)
print(myplot)

## 2. Per rank
order=c()
if ("Alpha Diversity" %in% unique(all_VCs_full$rank)){
  order=c(order,"Alpha Diversity")
}
if ("Beta Diversity" %in% unique(all_VCs_full$rank)){
  order=c(order,"Beta Diversity")
}
if ("Enterotypes" %in% unique(all_VCs_full$rank)){
  order=c(order,"Enterotypes")
}
if ("Clusters" %in% unique(all_VCs_full$rank)){
  order=c(order,"Clusters")
}
order=c(order,ranks)
for (subset_id in subset_ids){
  statistics=data.frame(rank=order)
  ### Traits
  filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  traits=count(filtered_VCs,rank)
  statistics$traits=traits$n[match(statistics[[1]],traits$rank)]
  ### Significant (FDR)
  sig_fdr=filtered_VCs[!filtered_VCs$Significance=="Non-significant traits",]
  stats_sig_herit=sig_fdr %>%                            
    group_by(rank) %>%
    summarize(min = min(var_Ad*100), median=median(var_Ad*100), max = max(var_Ad*100))
  sig_fdr=count(sig_fdr, rank)
  statistics$sig_fdr=sig_fdr$n[match(statistics[[1]],sig_fdr$rank)]
  statistics$percentage_fdr=100*statistics$sig_fdr/statistics$traits
  statistics$fdr_min_herit=stats_sig_herit$min[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$fdr_median_herit=stats_sig_herit$median[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$fdr_max_herit=stats_sig_herit$max[match(statistics[[1]],stats_sig_herit$rank)]
  ### Significant (Bonferroni)
  sig_bonf=filtered_VCs[filtered_VCs$Significance=="Significant (Bonferroni)",]
  stats_sig_herit=sig_bonf %>%                            
    group_by(rank) %>%
    summarize(min = min(var_Ad*100), median=median(var_Ad*100), max = max(var_Ad*100))
  sig_bonf=count(sig_bonf, rank)
  statistics$sig_bonf=sig_bonf$n[match(statistics[[1]],sig_bonf$rank)]
  statistics$percentage_bonf=100*statistics$sig_bonf/statistics$traits
  statistics$bonf_min_herit=stats_sig_herit$min[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$bonf_median_herit=stats_sig_herit$median[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$bonf_max_herit=stats_sig_herit$max[match(statistics[[1]],stats_sig_herit$rank)]
  statistics <- replace(statistics, is.na(statistics), 0)
  ### Make bar plot
  statistics_plot <- tidyr::pivot_longer(statistics, cols = c(traits, sig_fdr, sig_bonf), names_to = "Traits", values_to = "value")
  max_value=max(statistics_plot$value)
  myplot <- ggplot(statistics_plot, aes(x = factor(rank, levels = order), y = value, fill = Traits)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13)) +
    coord_cartesian(ylim = c(0, max_value+100)) +
    labs(x= "", y = paste0("Number of traits (Subset: ",subset_id,")")) +
    scale_fill_manual(values = c("traits" = "azure4", "sig_fdr" = "red", "sig_bonf" = "blue"),
                      labels = c("Significant (Bonferroni)", "Significant (FDR)", "All traits"),
                      name = "") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5)
  print(myplot)
  assign(paste("statistics_ranks",subset_id,sep="_"),statistics)
  print(subset_id)
}
dev.off()

# 7. Combined Effects
print("Combined effects on significant traits")
pdf('heritability_combined_effects.pdf',bg='white', width = 8, height = 6)
par(mar = c(2,2,4,2))
cbPalette <- c("lightblue","lightgrey", "lightcoral", "lightgreen")

for (subset_id in subset_ids){
  for (rank in ranks){
    filtered_VCs=all_VCs_full[!grepl("Non-significant",all_VCs_full$Significance),]
    filtered_VCs=filtered_VCs[grepl(subset_id,filtered_VCs$subset_id),]
    filtered_VCs=filtered_VCs[grepl(rank,filtered_VCs$rank),]
    filtered_VCs <- filtered_VCs[order(filtered_VCs$var_Ad,decreasing=T),]
    if (length(filtered_VCs$trait)>50){
      filtered_VCs=filtered_VCs[c(1:50),]
    }
    filtered_VCs <- filtered_VCs[order(filtered_VCs$var_Ad),]
    filtered_VCs$LBL[filtered_VCs$Significance == "Significant (Bonferroni)"] <- "**"
    filtered_VCs$LBL[filtered_VCs$Significance == "Significant (FDR)"] <- "*"
    filtered_VCs$trait <- factor(as.character(filtered_VCs$trait),levels=filtered_VCs$trait[order(filtered_VCs$var_Ad)])
    filtered_VCs_plot <- filtered_VCs[,c("trait","LBL","var_Ad","var_C", "var_M","var_Ed","STE_Ad")]
    filtered_VCs_plotL <- gather(filtered_VCs_plot,"Var.Exp","Var.Exp.NR", var_Ad:var_Ed,factor_key = T)
    filtered_VCs_plotL$Var.Exp <- as.character(filtered_VCs_plotL$Var.Exp)
    filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_Ed"] <- "Environment"
    filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_C"] <- "Cohousing"
    filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_M"] <- "Maternal Effects"
    filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_Ad"] <- "Host Aggregate Genetics"
    filtered_VCs_plotL$Var.Exp <- factor(as.character(filtered_VCs_plotL$Var.Exp),level = c("Environment","Cohousing","Maternal Effects","Host Aggregate Genetics"))
    myplot<-ggplot(filtered_VCs_plotL,aes(x=trait,y=Var.Exp.NR,fill=Var.Exp)) + 
      scale_fill_manual(values = cbPalette,name="Effect") +
      geom_col(col="grey", width=1,size=0.5) + 
      geom_errorbar(data = subset(filtered_VCs_plotL, Var.Exp == "Host Aggregate Genetics"),
                    mapping=aes(ymin=as.numeric(Var.Exp.NR)-as.numeric(STE_Ad),
                                ymax=as.numeric(Var.Exp.NR)+as.numeric(STE_Ad)),
                    width=0.25, linetype='solid') + 
      geom_text(data = filtered_VCs_plotL,
                aes(x = trait, y=STE_Ad,
                    label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
                color="black", vjust=+0.75, angle = 0, hjust=1,size=5) + 
      scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
      labs(title = paste0("Top 50 traits with significant heritability (Study: ",subset_id," / Rank: ",rank,")")) +
      ylab('Variance Explained') + xlab('') + 
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 7),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12),
            legend.position="bottom",
            plot.title = element_text(hjust = 0.5)) +
      coord_flip()
    print(myplot)
  }
}
dev.off()

if (method=="16S"){
  if ("ASV" %in% ranks){
    # 8. Check the importance of the taxonomic rank for heritability values on the Species level
    print("Checking the importance of the taxonomic rank for heritability values on the ASV level")
    anova_total=data.frame()
    pdf('heritability_asvs_by_rank.pdf',bg='white')
    test_ranks=c("Phylum",  "Class",   "Order",   "Family",  "Genus", "Species")
    for (subset_id in subset_ids){
      ## Filter ASV by study and assign higher taxonomic ranks
      filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
      filtered_VCs=filtered_VCs[grepl("ASV",filtered_VCs$rank),]
      filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$ASV)]
      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$ASV)]
      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$ASV)]
      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$ASV)]
      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$ASV)]
      filtered_VCs$Species=taxonomy$Species[match(filtered_VCs$trait,taxonomy$ASV)]
      ## ANOVA test on ranks
      filtered_VCs$norm_heritability = invrank(filtered_VCs$var_Ad)
      anova.herit <- summary(aov(norm_heritability ~ Phylum + Class + Order + Family + Genus + Species, data = filtered_VCs))
      anova.herit=anova.herit[[1]]
      anova.herit$R2=anova.herit$`Sum Sq`/sum(anova.herit$`Sum Sq`)
      anova.herit$subset_id=subset_id
      anova_total=rbind(anova_total,anova.herit)
      for (test_rank in test_ranks){
        ## Plot the proportion of ASVs with significant heritability by rank
        plot_data=filtered_VCs
        colnames(plot_data)[which(colnames(plot_data)==test_rank)]="test_rank"
        significance_by_rank <- plot_data %>%
          count(test_rank, Significance) %>%
          complete(test_rank, Significance, fill = list(n = 0)) %>%
          group_by(test_rank) %>%
          mutate(percent = n / sum(n) * 100)
        colnames(significance_by_rank)[1]="trait"
        order=significance_by_rank[grepl("Non-significant",significance_by_rank$Significance),]
        order = order[order(order$percent, decreasing = F),]
        if (length( order$trait)>50){
          order= order[1:50,]
        }
        order = order[order(order$percent, decreasing = T),]
        order=order$trait
        significance_by_rank = significance_by_rank[significance_by_rank$trait %in% order,]
        myplot<-ggplot(significance_by_rank, aes(factor(trait, levels = order), y = n, fill = Significance)) +
          geom_bar(stat = "identity", position = "stack") +
          #scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position="bottom",
                plot.title = element_text(hjust = 0.5)) +
          labs(title = paste0("Subset: ",subset_id," / Rank: ",test_rank),
               x = "Top 50 traits with the highest proportion of ASVs with significant heritability",
               y = "Number of ASVs",
               fill = "Significance") +
          coord_flip()
        print(myplot)
        ## Plot the comparisons between the median heritability of all ASVs on the rank with the own rank heritability
        median_var_Ad <- aggregate(var_Ad ~ ., data = filtered_VCs[, c("var_Ad", test_rank)], FUN = median)
        median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
        colnames(median_var_Ad)[1]="trait"
        if (length(median_var_Ad$trait)>10){
          median_var_Ad=median_var_Ad[1:10,]
        }
        order=median_var_Ad$trait
        filtered_VCs_rank=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
        filtered_VCs_rank=filtered_VCs_rank[grepl(test_rank,filtered_VCs_rank$rank),]
        filtered_VCs_rank=filtered_VCs_rank[filtered_VCs_rank$trait %in% median_var_Ad$trait,]
        filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (Bonferroni)" ] <- "**"
        filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (FDR)"] <- "*"
        filtered_VCs_rank$LBL[is.na(filtered_VCs_rank$LBL)] <- ""
        filtered_VCs_rank$heritability=paste(round(filtered_VCs_rank$var_Ad*100,digit=2),filtered_VCs_rank$LBL,sep="")
        filtered_VCs_rank$median=median_var_Ad$var_Ad[match(filtered_VCs_rank$trait,median_var_Ad$trait)]
        filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad > filtered_VCs_rank$median] <- "Rank Heritability > ASVs Median"
        filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad < filtered_VCs_rank$median] <- "Rank Heritability < ASVs Median"
        filtered_VCs_rank$comp[is.na(filtered_VCs_rank$comp)] <- ""
        plot_data=filtered_VCs[,c(which(colnames(filtered_VCs)=="var_Ad"),which(colnames(filtered_VCs)==test_rank))]
        colnames(plot_data)[which(colnames(plot_data)==test_rank)]="test_rank"
        plot_data=plot_data[plot_data$test_rank %in% median_var_Ad$trait,]
        # Plot
        myplot<-ggplot(plot_data, aes(x = factor(test_rank, levels = order), y = var_Ad*100)) +
          geom_boxplot() +
          geom_label(data = filtered_VCs_rank,  
                     aes(x = factor(trait, levels = order),
                         y = var_Ad*100, 
                         label = heritability, fill = comp),size = 3.5) +
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position="bottom") +
          guides(fill = guide_legend(override.aes = list(color = NA)), 
                 color = "none", 
                 shape = "none") +
          scale_fill_manual(values = c("Rank Heritability > ASVs Median" = "lightgreen", "Rank Heritability < ASVs Median" = "lightcoral"),
                            name = "") +
          labs(x = "Top 10 traits with higher ASVs median heritability", y=paste0("Heritability (%) (Subset: ",subset_id," / Rank: ",test_rank,")"))
        print(myplot)
      }
    }
    dev.off()
    save(anova_total,
         file="anova_rank_effect_on_heritability.RData")
  }
} else {
  if ("Species" %in% ranks){
    # 8. Check the importance of the taxonomic rank for heritability values on the Species level
    print("Checking the importance of the taxonomic rank for heritability values on the Species level")
    anova_total=data.frame()
    pdf('heritability_species_by_rank.pdf',bg='white')
    test_ranks=c("Phylum",  "Class",   "Order",   "Family",  "Genus")
    for (subset_id in subset_ids){
      ## Filter Species by study and assign higher taxonomic ranks
      filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
      filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
      filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
      ## ANOVA test on ranks
      filtered_VCs$norm_heritability = invrank(filtered_VCs$var_Ad)
      anova.herit <- summary(aov(norm_heritability ~ Phylum + Class + Order + Family + Genus, data = filtered_VCs))
      anova.herit=anova.herit[[1]]
      anova.herit$R2=anova.herit$`Sum Sq`/sum(anova.herit$`Sum Sq`)
      anova.herit$subset_id=subset_id
      anova_total=rbind(anova_total,anova.herit)
      for (test_rank in test_ranks){
        ## Plot the proportion of species with significant heritability by rank
        plot_data=filtered_VCs
        colnames(plot_data)[which(colnames(plot_data)==test_rank)]="test_rank"
        significance_by_rank <- plot_data %>%
          count(test_rank, Significance) %>%
          complete(test_rank, Significance, fill = list(n = 0)) %>%
          group_by(test_rank) %>%
          mutate(percent = n / sum(n) * 100)
        colnames(significance_by_rank)[1]="trait"
        order=significance_by_rank[grepl("Non-significant",significance_by_rank$Significance),]
        order = order[order(order$percent, decreasing = F),]
        if (length( order$trait)>50){
          order= order[1:50,]
        }
        order = order[order(order$percent, decreasing = T),]
        order=order$trait
        significance_by_rank = significance_by_rank[significance_by_rank$trait %in% order,]
        myplot<-ggplot(significance_by_rank, aes(factor(trait, levels = order), y = n, fill = Significance)) +
          geom_bar(stat = "identity", position = "stack") +
          #scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position="bottom",
                plot.title = element_text(hjust = 0.5)) +
          labs(title = paste0("Subset: ",subset_id," / Rank: ",test_rank),
               x = "Top 50 traits with the highest proportion of species with significant heritability",
               y = "Number of species",
               fill = "Significance") +
          coord_flip()
        print(myplot)
        ## Plot the comparisons between the median heritability of all species on the rank with the own rank heritability
        median_var_Ad <- aggregate(var_Ad ~ ., data = filtered_VCs[, c("var_Ad", test_rank)], FUN = median)
        median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
        colnames(median_var_Ad)[1]="trait"
        if (length(median_var_Ad$trait)>10){
          median_var_Ad=median_var_Ad[1:10,]
        }
        order=median_var_Ad$trait
        filtered_VCs_rank=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
        filtered_VCs_rank=filtered_VCs_rank[grepl(test_rank,filtered_VCs_rank$rank),]
        filtered_VCs_rank=filtered_VCs_rank[filtered_VCs_rank$trait %in% median_var_Ad$trait,]
        filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (Bonferroni)" ] <- "**"
        filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (FDR)"] <- "*"
        filtered_VCs_rank$LBL[is.na(filtered_VCs_rank$LBL)] <- ""
        filtered_VCs_rank$heritability=paste(round(filtered_VCs_rank$var_Ad*100,digit=2),filtered_VCs_rank$LBL,sep="")
        filtered_VCs_rank$median=median_var_Ad$var_Ad[match(filtered_VCs_rank$trait,median_var_Ad$trait)]
        filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad > filtered_VCs_rank$median] <- "Rank Heritability > Species Median"
        filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad < filtered_VCs_rank$median] <- "Rank Heritability < Species Median"
        filtered_VCs_rank$comp[is.na(filtered_VCs_rank$comp)] <- ""
        plot_data=filtered_VCs[,c(which(colnames(filtered_VCs)=="var_Ad"),which(colnames(filtered_VCs)==test_rank))]
        colnames(plot_data)[which(colnames(plot_data)==test_rank)]="test_rank"
        plot_data=plot_data[plot_data$test_rank %in% median_var_Ad$trait,]
        # Plot
        myplot<-ggplot(plot_data, aes(x = factor(test_rank, levels = order), y = var_Ad*100)) +
          geom_boxplot() +
          geom_label(data = filtered_VCs_rank,  
                     aes(x = factor(trait, levels = order),
                         y = var_Ad*100, 
                         label = heritability, fill = comp),size = 3.5) +
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position="bottom") +
          guides(fill = guide_legend(override.aes = list(color = NA)), 
                 color = "none", 
                 shape = "none") +
          scale_fill_manual(values = c("Rank Heritability > Species Median" = "lightgreen", "Rank Heritability < Species Median" = "lightcoral"),
                            name = "") +
          labs(x = "Top 10 traits with higher species median heritability", y=paste0("Heritability (%) (Subset: ",subset_id," / Rank: ",test_rank,")"))
        print(myplot)
      }
    }
    dev.off()
    save(anova_total,
         file="anova_rank_effect_on_heritability.RData")
  }
}

## 9. Comparisons between significant vs non-significant species or ASVs
#
#pdf('heritability_significant_vs_non_significant.pdf',bg='white')
#for (subset_id in subset_ids){
#  filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
#  filtered_VCs <- filtered_VCs %>% 
#    mutate( All_Sig = case_when(
#      Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Significant",
#      TRUE ~ "Non-significant"))
#  if (method=="16S"){
#    if ("ASV" %in% ranks){
#      print("Comparisons between significant and non-significant ASVs")
#      filtered_VCs=filtered_VCs[grepl("ASV",filtered_VCs$rank),]
#      filtered_VCs$Genome_size=catalogue_stats$Size[match(filtered_VCs$trait,catalogue_stats$ASV)]
#      filtered_VCs$Number_of_sequences=catalogue_stats$Number_of_sequences[match(filtered_VCs$trait,catalogue_stats$ASV)]
#      filtered_VCs$Mean_gene_size=catalogue_stats$mean[match(filtered_VCs$trait,catalogue_stats$ASV)]
#      filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Genome)]
#      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Genome)]
#      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Genome)]
#      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Genome)]
#      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Genome)]
#      rank="ASV"
#    }
#  } else {
#      if ("Species" %in% ranks){
#        print("Comparisons between significant and non-significant Species")
#        filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
#        filtered_VCs$Genome_size=catalogue_stats$Size[match(filtered_VCs$trait,catalogue_stats$Species)]
#        filtered_VCs$Number_of_sequences=catalogue_stats$Number_of_sequences[match(filtered_VCs$trait,catalogue_stats$Species)]
#        filtered_VCs$Mean_gene_size=catalogue_stats$mean[match(filtered_VCs$trait,catalogue_stats$Species)]
#        filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
#        filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
#        filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
#        filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
#        filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
#        rank="Species"
#      }
#    }
#  
#  # Calculate the number of significant (FDR) in each rank
#  stats_phylum=count(filtered_VCs,All_Sig,Phylum)
#  stats_class=count(filtered_VCs,All_Sig,Class)
#  stats_order=count(filtered_VCs,All_Sig,Order)
#  stats_family=count(filtered_VCs,All_Sig,Family)
#  stats_genus=count(filtered_VCs,All_Sig,Genus)
#  save(stats_phylum,
#       stats_class,
#       stats_order,
#       stats_family,
#       stats_genus, file="statistics_per_rank.RData")
#  
#  if (length(filtered_VCs[filtered_VCs$All_Sig=="Significant"])>0 ){
#    # Add relative abundance and prevalence
#    matrix=get(paste("filtered_prev",rank,subset_id,sep="_"))
#    ab_vs_prev=data.frame(trait=rownames(matrix),
#                          abundance=100*apply(sweep(matrix,2,colSums(matrix),"/"), 1, mean, na.rm=TRUE),
#                          prevalence=100*rowSums(matrix>0)/ncol(matrix))
#    filtered_VCs$Abundance=ab_vs_prev$abundance[match(filtered_VCs$trait,ab_vs_prev$trait)]
#    filtered_VCs$Prevalence=ab_vs_prev$prevalence[match(filtered_VCs$trait,ab_vs_prev$trait)]
#    # Make the plots
#    ## Comparing all types of Significance
#    order=c("Non-significant traits","Significant (FDR)","Significant (Bonferroni)")
#    ### Relative Abundance
#    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Abundance)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Relative Abundance (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Prevalence
#    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Prevalence)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Prevalence (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Genome Size
#    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Genome_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Genome Size (Total AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Number of genes
#    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Number_of_sequences)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Number of Genes") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Mean Gene Size
#    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Mean_gene_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Mean Gene Size (AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    
#    ## Comparing Significant (FDR) vs Non-Significant
#    order=c("Non-significant","Significant")
#    ### Relative Abundance
#    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=Abundance)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Relative Abundance (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Prevalence
#    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=Prevalence)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Prevalence (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Genome Size
#    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=Genome_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Genome Size (Total AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Number of genes
#    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=Number_of_sequences)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Number of Genes") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Mean Gene Size
#    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=Mean_gene_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Mean Gene Size (AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    
#    ## Taxonomy: Comparing the two genera with a higher number of species with significant heritability
#    top_2_genera <- filtered_VCs %>%
#      filter(All_Sig == "Significant") %>%
#      group_by(Genus) %>%
#      summarize(num_significant = n()) %>%
#      top_n(2, num_significant)
#    filtered_genus=filtered_VCs[filtered_VCs$Genus %in% top_2_genera$Genus,]
#    
#    ## Comparing Significant (FDR) vs Non-Significant (Only Top 2 Genera)
#    order=c("Non-significant","Significant")
#    ### Relative Abundance
#    max_Abundance=max(filtered_genus$Abundance)
#    myplot<-ggplot(filtered_genus, aes(x = factor(All_Sig, levels = order), y=Abundance)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Relative Abundance (%)") +
#      ggtitle(paste0("Subset: ",subset_id)) +
#      facet_wrap(~ Genus, scales = "free") +
#      ylim(0, max_Abundance) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Prevalence
#    max_Prevalence=max(filtered_genus$Prevalence)
#    myplot<-ggplot(filtered_genus, aes(x = factor(All_Sig, levels = order), y=Prevalence)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Prevalence (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      facet_wrap(~ Genus, scales = "free") +
#      ylim(0, max_Prevalence) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Genome Size
#    max_Genome_size=max(filtered_genus$Genome_size)
#    myplot<-ggplot(filtered_genus, aes(x = factor(All_Sig, levels = order), y=Genome_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Genome Size (Total AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      facet_wrap(~ Genus, scales = "free") +
#      ylim(0, max_Genome_size) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Number of genes
#    max_Number_of_sequences=max(filtered_genus$Number_of_sequences)
#    myplot<-ggplot(filtered_genus, aes(x = factor(All_Sig, levels = order), y=Number_of_sequences)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Number of Genes") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      facet_wrap(~ Genus, scales = "free") +
#      ylim(0, max_Number_of_sequences) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Mean Gene Size
#    max_Mean_gene_size=max(filtered_genus$Mean_gene_size)
#    myplot<-ggplot(filtered_genus, aes(x = factor(All_Sig, levels = order), y=Mean_gene_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Mean Gene Size (AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      facet_wrap(~ Genus, scales = "free") +
#      ylim(0, max_Mean_gene_size) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c("Non-significant","Significant")),
#                  map_signif_level = TRUE)
#    print(myplot)
#    
#    ## Comparing Only Top Two Genera
#    ### Heritability
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,100*var_Ad,FUN = median), y=100*var_Ad)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Heritability (%)") +
#      ggtitle(paste0("Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Relative Abundance
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,Abundance,FUN = median), y=Abundance)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Relative Abundance (%)") +
#      ggtitle(paste0("Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Prevalence
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,Prevalence,FUN = median), y=Prevalence)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Prevalence (%)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Genome Size
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,Genome_size,FUN = median), y=Genome_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Genome Size (Total AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Number of genes
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,Number_of_sequences,FUN = median), y=Number_of_sequences)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Number of Genes") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#    ### Mean Gene Size
#    myplot<-ggplot(filtered_genus, aes(x=reorder(Genus,Mean_gene_size,FUN = median), y=Mean_gene_size)) +
#      geom_boxplot() +
#      xlab("") +
#      ylab("Mean Gene Size (AAs)") +
#      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
#      theme_bw() +
#      theme(legend.position="bottom",
#            axis.line = element_line(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.border = element_blank(),
#            panel.background = element_blank(),
#            axis.text = element_text(size = 15),
#            axis.title = element_text(size = 15),
#            legend.text=element_text(size=10),
#            plot.title = element_text(hjust = 0.5))+
#      scale_x_discrete(guide = guide_axis(angle = 45)) +
#      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
#                  map_signif_level = TRUE)
#    print(myplot)
#  }
#}
#dev.off()

### End of Step 3 ###

### Step 4 ####

if (module_analysis!="Pooled"){
  print("Step 4: Calculate common significant taxa across subsets")
  
  ## 1. Common across studies - FDR
  stats_sig_fdr=all_VCs_full[!all_VCs_full$Significance=="Non-significant traits",]
  stats_sig_subset_id=stats_sig_fdr[!grepl("ALL",stats_sig_fdr$subset_id),]
  filtered_stats_sig_subset_id_fdr <- stats_sig_subset_id %>%
    group_by(trait) %>%
    filter(n() > 1)
  common_fdr=as.data.frame(unique(filtered_stats_sig_subset_id_fdr$trait))
  colnames(common_fdr)="Common traits across studies"
  
  ## 2. Common across studies - Bonferroni
  stats_sig_bonf=all_VCs_full[all_VCs_full$Significance=="Significant (Bonferroni)",]
  stats_sig_subset_id=stats_sig_bonf[!grepl("ALL",stats_sig_bonf$subset_id),]
  filtered_stats_sig_subset_id_bonf <- stats_sig_subset_id %>%
    group_by(trait) %>%
    filter(n() > 1)
  common_bonf=as.data.frame(unique(filtered_stats_sig_subset_id_bonf$trait))
  colnames(common_bonf)="Common traits across studies"
} else {
  common_fdr=c()
  common_bonf=c()
}
save(common_fdr,common_bonf,file="common_traits.RData")

### End of Step 4 ####