library(xfun)
library(compositions)
library(iNEXT.3D)
suppressMessages(library(vegan))
suppressMessages(library(ape))
library(doParallel)
library(tidyr)
suppressMessages(library(ggplot2))
suppressMessages(library(philr))
library(ggsignif)

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

ranks=c("Phylum","Class","Order","Family","Genus","Species")
subset_ids=c("ALL")
module_comparisons="YES"
methods=c("16S","Shallow")
harmonization=c("sample_harmonization","sample_and_taxa_harmonization")

sample_identifier= "host_subject_id"
input_meta="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt"

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

#Import taxonomy and phylogeny
taxonomy<- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt", header=FALSE)
colnames(taxonomy)=c("Genome","Taxonomy")
taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
taxonomy$Phylum<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",2)
taxonomy$Class<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",3)
taxonomy$Order<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",4)
taxonomy$Family<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",5)
taxonomy$Genus<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",6)
taxonomy$Species<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",7)
tree = ape::read.tree("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/phylotree_GTDB207.tree")
# Replace species name by genome
# Load count tables after sample and tax harmonization
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Harmonization/harmonized_samples_and_taxa.RData")
#species_all_total=filtered_prev_objs[[12]]
#rownames(species_all_total)=taxonomy$Genome[match(rownames(species_all_total),taxonomy$Species)]
#m=match(rownames(species_all_total),tree$tip.label)
#any(is.na(m))


harmonization="samples_and_taxa_harmonization"
harm=harmonization

#pdf("/users/abaud/fmorillo/paper_figures/16S_SS_diversity_correlations.pdf")
#for (harm in harmonization){
  
  # Load counts
  if (harm=="sample_harmonization") {
    load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Harmonization/harmonized_samples.RData")
    title="Only samples harmonized"
  } else {
    load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Harmonization/harmonized_samples_and_taxa.RData")
    title="Samples and taxa harmonized"
  }
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
  
  # Parameters for alpha and beta
  subset_id="ALL"
  covariate="Study"
  div_type="PD"
  
  # Beta Diversity
  offset = 1
  #matrix=filtered_prev_Species_ALL_Shallow
  ranks=c("Species")
  for (rank in ranks){
    permanova_total=data.frame()
    #beta_cor=data.frame(samples=colnames(matrix))
    for (method in methods){
      matrix=get(paste('filtered_prev',rank,subset_id,method,sep="_"))
      rownames(matrix)=taxonomy$Genome[match(rownames(matrix),taxonomy$Species)]
      matrix=(round(matrix, digits = 0))
      # Prune the phylogenetic tree
      keep_taxa=rownames(matrix)
      remove_taxa = setdiff(tree$tip.label, keep_taxa)
      pruned_tree = drop.tip(tree, remove_taxa)
      pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
      matrix=matrix[rownames(matrix) %in% pruned_tree$tip.label,]
      # Calculate dissimilarities in microbiome composition
      colnames(matrix)=metadata$host_subject_id[match(colnames(matrix),metadata$RFID)]
      this_metadata=metadata[match(colnames(matrix),metadata$host_subject_id),]
      # Calculate Beta Diversity
      matrix=matrix+offset
      gp.philr <- philr(t(matrix), pruned_tree, 
                        part.weights='enorm.x.gm.counts', 
                        ilr.weights='blw.sqrt')
      gp.dist <- dist(gp.philr, method="euclidean")
      genome.pcoa <- ape::pcoa(gp.dist)
      #if (method==methods[1]) {
      #  prop_var1 <- genome.pcoa$values$Relative_eig
      #} else {
      #  prop_var2 <- genome.pcoa$values$Relative_eig
      #  plot(prop_var1, type = "b", col = "blue", xlab = "Principal Coordinate", ylab = "Proportion of Variance", ylim=c(0, 0.60), main = paste0("Variance Explained by PCs (Rank: ",rank,")"))
      #  points(prop_var2, type = "b", col = "red")
      #  legend("topright", legend = c("16S","Shallow Shotgun"), col = c("blue", "red"), pch = 1)
      #}
      genome.pc <- data.frame(genome.pcoa$vectors)
      #beta_cor[[method]]=genome.pc$Axis.1
      ## Correlations
      #if (method==methods[2]){
      #  plot(as.numeric(beta_cor$`16S`),as.numeric(beta_cor$Shallow),
      #       xlab="16S",
      #       ylab="Shallow Shotgun",
      #       cex.lab = 1.3,
      #       main=paste0("Beta Diversity (PCo1): ",title," / ",rank))
      #  cor_result=cor.test(as.numeric(beta_cor$`16S`),as.numeric(beta_cor$Shallow),method="spearman")
      #  r_value <- cor_result$estimate
      #  p_value <- cor_result$p.value
      #  if (p_value < 0.01){
      #    P="P < 0.01"
      #  } else {
      #    if (p_value > 0.01 && p_value < 0.05){
      #      P="P < 0.05"
      #    } else {
      #      P="P > 0.05"
      #    }
      #  }
      #  legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
      #  new_df <- beta_cor %>%
      #    pivot_longer(cols = c(`16S`, Shallow), 
      #                 names_to = "method", 
      #                 values_to = "beta")
      #  order=c("16S", "Shallow")
      #  myplot<-ggplot(new_df, aes(x=factor(method, levels = order), y=beta)) +
      #    geom_boxplot() +
      #    xlab("") +
      #    ylab("Beta Diversity (PCo1)") +
      #    ggtitle(paste0(title," / ",rank)) +
      #    theme_bw() +
      #    theme(legend.position="bottom",
      #          axis.line = element_line(),
      #          panel.grid.major = element_blank(),
      #          panel.grid.minor = element_blank(),
      #          panel.border = element_blank(),
      #          panel.background = element_blank(),
      #          axis.text = element_text(size = 13),
      #          axis.title = element_text(size = 15),
      #          legend.text=element_text(size=12),legend.title=element_text(size=15),
      #          plot.title = element_text(hjust = 0.5))+
      #    scale_x_discrete(guide = guide_axis(angle = 60)) +
      #    geom_signif(comparisons = list(c("16S", "Shallow")),
      #                map_signif_level = TRUE)
      #  print(myplot)
      #}
      # Run PERMANOVA
      formula <- as.formula(paste0("gp.dist ~ ", covariate))
      permanova.beta <- adonis2(formula, data = this_metadata, permutations = 999)
      permanova.beta$method=method
      permanova.beta$subset=subset_id
      permanova.beta$covariate=covariate
      permanova.beta$diversity=div_type
      permanova.beta$number_of_taxa=nrow(matrix)
      permanova_total=rbind(permanova_total,permanova.beta)
      assign(paste("permanova",rank,harm,sep="_"),permanova_total)
      print(paste0("PERMANOVA: ",harm," / ",rank," / ",method))
    }
  }
  
  ## Alpha Diversity
  #hill_number=1
  #alpha_n= 50
  #alpha_c= 0.95
  #matrix=filtered_prev_Species_ALL_Shallow
  #for (rank in ranks[-length(ranks)]){
  #  anova_total=data.frame()
  #  alpha_cor=data.frame(samples=colnames(matrix))
  #  for (method in methods){
  #    matrix=get(paste('filtered_prev',rank,subset_id,method,sep="_"))
  #    matrix=(round(matrix, digits = 0))
  #    colnames(matrix)=metadata$host_subject_id[match(colnames(matrix),metadata$RFID)]
  #    this_metadata=metadata[match(colnames(matrix),metadata$host_subject_id),]
      ## Calculate Beta Diversity
      #cl <- makeCluster(detectCores())
      #print(cl)
      #registerDoParallel(cl)
      #alpha_table <- foreach(i = 1:ncol(matrix), .combine = rbind) %dopar% {
      #  col_data <- matrix[,i]
      #  names(col_data)=rownames(matrix)
      #  iNEXT.3D::AO3D(data=col_data,
      #                 diversity = div_type,
      #                 q = as.numeric(hill_number),
      #                 datatype = "abundance",
      #                 nboot = alpha_n,
      #                 conf = alpha_c,
      #                 method = "Asymptotic")
      #}
      #stopImplicitCluster()
      #alpha_table$Assemblage=colnames(matrix)
      ##colnames(alpha_table)[3]="qD"
      #alpha_table$Sex=this_metadata$Sex[match(alpha_table$Assemblage,this_metadata$host_subject_id)]
      #alpha_table$Batch=this_metadata$Batch[match(alpha_table$Assemblage,this_metadata$host_subject_id)]
      #alpha_table$Study=this_metadata$Batch[match(alpha_table$Assemblage,this_metadata$host_subject_id)]
      #alpha_cor[[method]]=alpha_table$qD
      ## Correlations
      #if (method==methods[2]){
      #  plot(as.numeric(alpha_cor$`16S`),as.numeric(alpha_cor$Shallow),
      #       xlab="16S",
      #       ylab="Shallow Shotgun",
      #       cex.lab = 1.3,
      #       main=paste0("Alpha diversity (Shannon Index): ",title," / ",rank))
      #  cor_result=cor.test(as.numeric(alpha_cor$`16S`),as.numeric(alpha_cor$Shallow),method="spearman")
      #  r_value <- cor_result$estimate
      #  p_value <- cor_result$p.value
      #  if (p_value < 0.01){
      #    P="P < 0.01"
      #  } else {
      #    if (p_value > 0.01 && p_value < 0.05){
      #      P="P < 0.05"
      #    } else {
      #      P="P > 0.05"
      #    }
      #  }
      #  legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
      #  new_df <- alpha_cor %>%
      #    pivot_longer(cols = c(`16S`, Shallow), 
      #                 names_to = "method", 
      #                 values_to = "alpha")
      #  order=c("16S", "Shallow")
      #  myplot<-ggplot(new_df, aes(x=factor(method, levels = order), y=alpha)) +
      #    geom_boxplot() +
      #    xlab("") +
      #    ylab("Alpha diversity (Shannon Index)") +
      #    ggtitle(paste0(title," / ",rank)) +
      #    theme_bw() +
      #    theme(legend.position="bottom",
      #          axis.line = element_line(),
      #          panel.grid.major = element_blank(),
      #          panel.grid.minor = element_blank(),
      #          panel.border = element_blank(),
      #          panel.background = element_blank(),
      #          axis.text = element_text(size = 13),
      #          axis.title = element_text(size = 15),
      #          legend.text=element_text(size=12),legend.title=element_text(size=15),
      #          plot.title = element_text(hjust = 0.5))+
      #    scale_x_discrete(guide = guide_axis(angle = 60)) +
      #    geom_signif(comparisons = list(c("16S", "Shallow")),
      #                map_signif_level = TRUE)
      #  print(myplot)
      #}
      # Run ANOVA
     #alpha_table$norm_qD = invrank(alpha_table$qD)
     #formula <- as.formula(paste0("norm_qD ~ ", "Sex + Batch + Study"))
     #anova.alpha <- summary(aov(formula, data = alpha_table))
     #anova.alpha=anova.alpha[[1]]
     #anova.alpha$R2=anova.alpha$`Sum Sq`/sum(anova.alpha$`Sum Sq`)
     #anova.alpha$method=method
     #anova.alpha$subset=subset_id
     #anova.alpha$covariate="all"
     #anova.alpha$diversity=div_type
     #anova.alpha$q=hill_number
     #anova.alpha$number_of_taxa=nrow(matrix)
     #anova_total=rbind(anova_total,anova.alpha)
     #assign(paste("anova",rank,harm,sep="_"),anova_total)
     #print(paste0("ANOVA: ",harm," / ",rank," / ",method))
#}
dev.off()

permanova_filt=data.frame()
anova_filt=data.frame()
for (harm in harmonization){
  if (harm=="sample_harmonization") {
    title="Only samples harmonized"
  } else {
    title="Samples and taxa harmonized"
  }
  for (rank in ranks[-length(ranks)]){
    # PERMANOVA
    permanova=get(paste('permanova',rank,harm,sep="_"))
    permanova=permanova[c(1,4), c(3,6,ncol(permanova))]
    permanova$rank=rank
    permanova$harmonization=title
    permanova_filt=rbind(permanova_filt,permanova)
    # ANOVA
    anova=get(paste('anova',rank,harm,sep="_"))
    anova=anova[c(1,3), c(6,7,ncol(anova))]
    anova$rank=rank
    anova$harmonization=title
    anova_filt=rbind(anova_filt,anova)
    print(paste0(harm," / ",rank))
  }
}

save(permanova_Phylum_sample_harmonization,
     permanova_Class_sample_harmonization,
     permanova_Order_sample_harmonization,
     permanova_Family_sample_harmonization,
     permanova_Genus_sample_harmonization,
     permanova_Species_sample_harmonization,
     permanova_Phylum_sample_and_taxa_harmonization,
     permanova_Class_sample_and_taxa_harmonization,
     permanova_Order_sample_and_taxa_harmonization,
     permanova_Family_sample_and_taxa_harmonization,
     permanova_Genus_sample_and_taxa_harmonization,
     permanova_Species_sample_and_taxa_harmonization,
     anova_Phylum_sample_harmonization,
     anova_Class_sample_harmonization,
     anova_Order_sample_harmonization,
     anova_Family_sample_harmonization,
     anova_Genus_sample_harmonization,
     anova_Species_sample_harmonization,
     anova_Phylum_sample_and_taxa_harmonization,
     anova_Class_sample_and_taxa_harmonization,
     anova_Order_sample_and_taxa_harmonization,
     anova_Family_sample_and_taxa_harmonization,
     anova_Genus_sample_and_taxa_harmonization,
     anova_Species_sample_and_taxa_harmonization,
     permanova_filt,
     anova_filt,
  file="/users/abaud/fmorillo/paper_figures/alpha_and_beta.RData"
)
