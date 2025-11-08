heatmap_result <- pheatmap(
  presence_absence,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  #main=paste0(func_type,": Mean number per ",genus," species"),
  legend = TRUE,
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  annotation_name_row = FALSE)


load("~/NextFlow/Microbiome_Profiling/step2/Output_Shallow/Heritability/heritability.RData")
all_VCs_full$Heritability=100*all_VCs_full$var_Ad

subset_id="ALL"
filtered_VCs_ALL=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs_ALL=filtered_VCs_ALL[grepl("Species",filtered_VCs_ALL$rank),]

subset_id="MI"
filtered_VCs_MI=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs_MI=filtered_VCs_MI[grepl("Species",filtered_VCs_MI$rank),]
filtered_VCs_MI_Prevotella=filtered_VCs_MI[grepl("s__Prevotella",filtered_VCs_MI$trait),]
filtered_VCs_MI_Bacteroides=filtered_VCs_MI[grepl("s__Bacteroides",filtered_VCs_MI$trait),]
filtered_VCs_MI_Bacteroides=filtered_VCs_MI_Bacteroides[!grepl("Bacteroides_F",filtered_VCs_MI_Bacteroides$trait),]

subset_id="NY"
filtered_VCs_NY=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs_NY=filtered_VCs_NY[grepl("Species",filtered_VCs_NY$rank),]
filtered_VCs_NY_Prevotella=filtered_VCs_NY[grepl("s__Prevotella",filtered_VCs_NY$trait),]
filtered_VCs_NY_Bacteroides=filtered_VCs_NY[grepl("s__Bacteroides",filtered_VCs_NY$trait),]
filtered_VCs_NY_Bacteroides=filtered_VCs_NY_Bacteroides[!grepl("Bacteroides_F",filtered_VCs_NY_Bacteroides$trait),]

threshold=6
filtered_VCs_MI_Prevotella <- filtered_VCs_MI_Prevotella %>%
  arrange(desc(Heritability)) %>%
  mutate(
    Heritability_Rank = case_when(
      row_number() <= threshold ~ "High Heritability",
      row_number() > n() - threshold ~ "Low Heritability",
      TRUE ~ "others"
    )
  )

common_bac=intersect(filtered_VCs_MI_Bacteroides$trait,filtered_VCs_NY_Bacteroides$trait)
common_prev=intersect(filtered_VCs_MI_Prevotella$trait,filtered_VCs_NY_Prevotella$trait)

common_NY_Prevotella=filtered_VCs_NY_Prevotella[filtered_VCs_NY_Prevotella$trait %in% common_prev,]
common_MI_Prevotella=filtered_VCs_MI_Prevotella[filtered_VCs_MI_Prevotella$trait %in% common_prev,]
common_NY_Bacteroides=filtered_VCs_NY_Bacteroides[filtered_VCs_NY_Bacteroides$trait %in% common_bac,]
common_MI_Bacteroides=filtered_VCs_MI_Bacteroides[filtered_VCs_MI_Bacteroides$trait %in% common_bac,]
boxplot(common_MI_Bacteroides$Heritability,common_NY_Bacteroides$Heritability,names=c("MI_Bacteroides","NY_Bacteroides"),ylab="Heritability(%)")
boxplot(common_MI_Prevotella$Heritability,common_NY_Prevotella$Heritability,names=c("MI_Prevotella","NY_Prevotella"),ylab="Heritability(%)")
rownames(common_MI_Prevotella)=common_MI_Prevotella$trait
rownames(common_NY_Prevotella)=common_NY_Prevotella$trait
common_MI_Prevotella=common_MI_Prevotella[rownames(common_NY_Prevotella),]

Common_Prevotella=data.frame(trait=common_MI_Prevotella$trait)
Common_Prevotella$Heritability_MI=common_MI_Prevotella$Heritability
Common_Prevotella$Heritability_NY=common_NY_Prevotella$Heritability
Common_Prevotella$Difference=Common_Prevotella$Heritability_NY-Common_Prevotella$Heritability_MI

Common_Bacteroides=data.frame(trait=common_MI_Bacteroides$trait)
Common_Bacteroides$Heritability_MI=common_MI_Bacteroides$Heritability
Common_Bacteroides$Heritability_NY=common_NY_Bacteroides$Heritability
Common_Bacteroides$Difference=Common_Bacteroides$Heritability_NY-Common_Bacteroides$Heritability_MI

threshold=6
Common_Prevotella <- Common_Prevotella %>%
  arrange(desc(Difference)) %>%
  mutate(
    Heritability_Rank = case_when(
      row_number() <= threshold ~ "Highest Increase",
      row_number() > n() - threshold ~ "Lowest Increase",
      TRUE ~ "others"
    )
  )


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

method="Shallow"
subset_id="ALL"

# Read taxonomy
print("Reading taxonomy info")
input_taxonomy="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt"
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

# Read phylogenetic tree info
phylotree="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/phylotree_GTDB207.tree"
if (phylotree == ""){
  stop("No phylogeny provided. Provide it if you want to calculate phylogenetic alpha diversity or any type of beta diversity (Phylr has to load a phylogenetic tree, even for the calculation of the taxonomic beta diversity.)")
} else {
  tree = ape::read.tree(phylotree)
}

load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_CAZy.RData")
matrix=func_matrix_CAZy
colnames(matrix)=taxonomy$Species[match(colnames(matrix),taxonomy$Genome)]
matrix=matrix[,colnames(matrix) %in% Common_Prevotella$trait]
Common_Prevotella$Genome=taxonomy$Genome[match(Common_Prevotella$trait,taxonomy$Species)]
this_Common_Prev <- Common_Prevotella[grepl("Highest Increase", Common_Prevotella$Heritability_Rank) | grepl("Lowest Increase", Common_Prevotella$Heritability_Rank), ]
matrix=matrix[,this_Common_Prev$trait]

matrix=t(matrix)

rownames(matrix)=this_Common_Prev$Heritability_Rank[match(rownames(matrix),this_Common_Prev$trait)]
collapsed_matrix=t(sapply(by(matrix,rownames(matrix),colSums),identity))
collapsed_matrix=t(collapsed_matrix)

matrix_relative <- sweep(collapsed_matrix,2,threshold,"/")
differences=data.frame(Function=rownames(matrix_relative),Difference=matrix_relative[,1]-matrix_relative[,2])
differences<- differences[differences$Difference != 0, ]
matrix_relative_filt=matrix_relative[rownames(matrix_relative) %in% differences$Function,]

presence_absence <- ifelse(collapsed_matrix > 0, 1, 0)

keep_taxa=Common_Prevotella$Genome
any(is.na(keep_taxa))
remove_taxa = setdiff(tree$tip.label, keep_taxa)
any(is.na(remove_taxa))
pruned_tree = drop.tip(tree, remove_taxa)
pruned_tree$tip.label <- taxonomy$Species[match(pruned_tree$tip.label,taxonomy$Genome)]
pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
significance_colors <- ifelse(Common_Prevotella$Heritability_Rank == "Highest Increase", "blue",
                              ifelse(Common_Prevotella$Heritability_Rank == "Lowest Increase","red","black"))
plot.phylo(pruned_tree, cex = 0.6, type = "phylogram", tip.color = significance_colors, show.tip.label = TRUE, show.node.label = FALSE, main = "")
tip_labels=c("Highest increase in Heritability","Lowest increase in Heritability")
significance_colors=c("blue","red")
legend("bottom", legend = tip_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))

annotation=data.frame(Increase_in_Heritability = as.factor(this_Common_Prev$Heritability_Rank))
rownames(annotation)=rownames(matrix)



# Step 5: Create the heatmap
pheatmap(
  matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  #scale = "column",
  annotation_row = annotation,
  annotation_legend = TRUE,
  #annotation_colors = list(Significance = sig_colors),
  #clustering_distance_rows = "correlation",  # Use correlation-based distance
  #correlation_matrix = filtered_SPARCC,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 7,
  legend = TRUE,  # Show legend for color scale
  #legend_title = "Cluster",  # Title for cluster legend
  #legend_labels = list(levels(cluster_colors), levels(sig_colors)),  # Labels for cluster legend
  #legend_colors = list(cluster_cols, sig_cols)  # Specify the colors for the cluster legend
)
dev.off()

######

subset_id="MI"
clusters_MI <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Cluster_Analyses/Cluster_UMAP_HDBSCAN_",subset_id,".csv"))
colnames(clusters_MI)[1]="Genome"
MI_RES <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Matrix_Processing/matrices_residues/",subset_id,".csv"))
clusters_MI$Genome=MI_RES$Genome_id
clusters_MI$Cluster=gsub("-1","noise",clusters_MI$Cluster)
clusters_MI$Cluster=paste("cluster",clusters_MI$Cluster,sep="__")
clusters_MI$Species=taxonomy$Species[match(clusters_MI$Genome,taxonomy$Genome)]
stats=count(clusters_MI,Cluster)

subset_id="NY"
clusters_NY <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Cluster_Analyses/Cluster_UMAP_HDBSCAN_",subset_id,".csv"))
colnames(clusters_NY)[1]="Genome"
MI_RES <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Matrix_Processing/matrices_residues/",subset_id,".csv"))
clusters_NY$Genome=MI_RES$Genome_id
clusters_NY$Cluster=gsub("-1","noise",clusters_NY$Cluster)
clusters_NY$Cluster=paste("cluster",clusters_NY$Cluster,sep="__")
clusters_NY$Species=taxonomy$Species[match(clusters_NY$Genome,taxonomy$Genome)]
stats=count(clusters_NY,Cluster)


Common_Prevotella$Clusters_MI=clusters_MI$Cluster[match(Common_Prevotella$trait,clusters_MI$Species)]
Common_Prevotella$Clusters_NY=clusters_NY$Cluster[match(Common_Prevotella$trait,clusters_NY$Species)]

Common_Bacteroides$Clusters_MI=clusters_MI$Cluster[match(Common_Bacteroides$trait,clusters_MI$Species)]
Common_Bacteroides$Clusters_NY=clusters_NY$Cluster[match(Common_Bacteroides$trait,clusters_NY$Species)]
