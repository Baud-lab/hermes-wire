suppressMessages(library(vegan))
suppressMessages(library(ape))
#suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(philr))

# Load count tables after sample and tax harmonization
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Harmonization/harmonized_samples_and_taxa.RData")
species_all_Shallow=filtered_prev_objs[[12]]
colnames(species_all_Shallow)=paste(colnames(species_all_Shallow),"Shallow",sep="_")
species_all_16S=filtered_prev_objs[[6]]
colnames(species_all_16S)=paste(colnames(species_all_16S),"16S",sep="_")
# Build combined matrix
species_all_total=cbind(species_all_Shallow,species_all_16S)
# Build metadata
metadata=data.frame(Sample=c(colnames(species_all_Shallow),colnames(species_all_16S)))
metadata$Method=sapply(strsplit(metadata$Sample, "_"), "[",2)
metadata$RFID=sapply(strsplit(metadata$Sample, "_"), "[",1)
metadata_Shallow <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")
metadata$Study = metadata_Shallow$Study[match(metadata$RFID,metadata_Shallow$RFID)]
any(is.na(metadata))
# Import taxonomy and phylogeny
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
rownames(species_all_total)=taxonomy$Genome[match(rownames(species_all_total),taxonomy$Species)]
m=match(rownames(species_all_total),tree$tip.label)
any(is.na(m))

# Calculate Beta Diversity
gp.dist_objs=list()
beta_objs=list()
matrix=(round(species_all_total, digits = 0))
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
div_types=c("PD")
for (div_type in div_types){
  if (div_type=="PD"){
    gp.philr <- philr(t(filtered_trait_off), pruned_tree, 
                      part.weights='enorm.x.gm.counts', 
                      ilr.weights='blw.sqrt')
    gp.dist <- dist(gp.philr, method="euclidean")
  } else {
    gp.philr <- philr(t(filtered_trait_off), pruned_tree, 
                      part.weights='enorm.x.gm.counts', 
                      ilr.weights='uniform')
    gp.dist <- dist(gp.philr, method="euclidean")
  }
  # Calculate PCoA and extract the principal coordinates (PCs)
  gp.dist_obj <- paste('gp.dist',div_type,sep='_')
  assign(gp.dist_obj, gp.dist)
  gp.dist_objs<- append(gp.dist_objs, list(gp.dist_obj=gp.dist))
  print(paste0(div_type))
}
### Plot PCoA
pdf('/users/abaud/fmorillo/paper_figures/16S_Shallow_beta.pdf',bg='white')
# Filter covariates and samples by subset
covariates=c("Method","Study")
matrix=(round(species_all_total, digits = 0))
motch = match(colnames(matrix), metadata$Sample)
any(is.na(motch))
# Calculate PCoA
for (covariate in covariates){
  for (div_type in div_types){
    if (div_type=="PD"){
      beta="Phylogenetic"
    } else {
      beta="Taxonomic"
    }
    # Calculate PCoA and extract the principal coordinates (PCs)
    genome.pcoa <- ape::pcoa(get(paste("gp.dist",div_type,sep="_")))
    genome.pc <- data.frame(genome.pcoa$vectors)
    # Make plot
    myplot<-ggplot(genome.pc, aes(x=Axis.1, y=Axis.2, colour = as.factor(metadata[[covariate]]))) +
      geom_point() + 
      ggtitle(paste0("Cov: ",covariate, " / Dist. Type: ",beta)) +
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
dev.off()

### Calculate PERMANOVA
permanova_total=data.frame()
# Filter covariates and samples by subset
covariates=c("Method","Study")
covariates=c(covariates,"all")
matrix=(round(species_all_total, digits = 0))
motch = match(colnames(matrix), metadata$Sample) #check colnames
any(is.na(motch))
for (div_type in div_types){
  gp.dist=get(paste("gp.dist",div_type,sep="_"))
  for (covariate in covariates){
    if (covariate=="all"){
      formula <- as.formula(paste("gp.dist ~", paste(covariates[-length(covariates)], collapse = " + ")))
    } else {
      formula <- as.formula(paste("gp.dist ~", paste(covariate, collapse = " + ")))
    }
    permanova.beta <- adonis2(formula, data = metadata, permutations = 999)
    permanova.beta$covariate=covariate
    permanova.beta$diversity=div_type
    permanova_total=rbind(permanova_total,permanova.beta)
    print(paste0(covariate,"/",div_type))
  }
}
save(permanova_total,file="/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/Rnor6/taxonomic_profile/comp_shallow_16S/samples_and_taxa_harm/permanova.RData")