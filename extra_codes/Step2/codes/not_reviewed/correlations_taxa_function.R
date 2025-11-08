subset_id="ALL"
method="Shallow"

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

# Load Clusters
Cluster_UMAP_HDBSCAN_ALL <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Cluster_Analyses/Cluster_UMAP_HDBSCAN_ALL.csv")
clusters=Cluster_UMAP_HDBSCAN_ALL
colnames(clusters)[1]="Genome"
ALL_RES <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Matrix_Processing/matrices_residues/ALL.csv")
clusters$Genome=ALL_RES$Genome_id
clusters$Cluster=gsub("-1","noise",clusters$Cluster)
clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")

# Load heritability
heritability="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Heritability/heritability.RData"
load(heritability)
filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
filtered_VCs <- filtered_VCs %>% 
  mutate( All_Sig = case_when(
    Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Significant",
    TRUE ~ "Not significant"))
filtered_VCs$Genome=taxonomy$Genome[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Cluster=clusters$Cluster[match(filtered_VCs$Genome,clusters$Genome)]

filtered_VCs$Sig_Phy=paste(filtered_VCs$Phylum,filtered_VCs$All_Sig,sep="_")

# Filter
filtered_VCs=filtered_VCs[grepl("p__Bacteroidota",filtered_VCs$Phylum) | grepl("p__Firmicutes_A",filtered_VCs$Phylum),]

# Filter by Species
filtered_VCs=filtered_VCs[grepl("s__Prevotella",filtered_VCs$trait),]

filtered_VCs=filtered_VCs[grepl("s__Bacteroides",filtered_VCs$trait),]

filtered_VCs=filtered_VCs[grepl("s__Prevotella",filtered_VCs$trait) | grepl("s__Bacteroides",filtered_VCs$trait),]
filtered_VCs=filtered_VCs[!grepl("s__Bacteroides_F",filtered_VCs$trait),]

# Load taxonomy
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/Cluster_Analyses/residuals_qned_counts.RData")
#phylum_ALL=residuals_qned_counts_objs[[3]]
#family_ALL=residuals_qned_counts_objs[[3]]
#genus_ALL=residuals_qned_counts_objs[[3]]

species_ALL=residuals_qned_counts_objs[[18]]

# Filter
species_ALL=species_ALL[rownames(species_ALL) %in% filtered_VCs$trait,]

# Filter by species
species_ALL=species_ALL[grepl("s__Prevotella",rownames(species_ALL)),]

species_ALL=species_ALL[grepl("s__Bacteroides",rownames(species_ALL)),]

species_ALL=species_ALL[grepl("s__Prevotella",rownames(species_ALL)) | grepl("s__Bacteroides",rownames(species_ALL)),]
species_ALL=species_ALL[!grepl("s__Bacteroides_F",rownames(species_ALL)),]

# Load clusters
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Cluster_Analyses/filtered_prev_clusters.RData")
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
clusters_ALL=filtered_prev_cluster_objs[[number]]

# Load ECs
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input_func/EC_ALL_RES.RData")
residuals_ECs=residuals_ECs

species_ALL=species_ALL[,colnames(species_ALL) %in% colnames(residuals_ECs)]
species_ALL=species_ALL[,colnames(residuals_ECs)]

clusters_ALL=clusters_ALL[,colnames(clusters_ALL) %in% colnames(residuals_ECs)]
clusters_ALL=clusters_ALL[,colnames(residuals_ECs)]

matrix2=clusters_ALL
matrix1=clusters_ALL

matrix1=species_ALL
matrix2=residuals_ECs

matrix1=clusters_ALL
matrix2=residuals_ECs

matrix1=clusters_ALL
matrix2=species_ALL



correlation_matrix <- cor(t(matrix1), t(matrix2), method="spearman")

rownames(filtered_VCs)=filtered_VCs$trait
filtered_VCs=filtered_VCs[rownames(species_ALL),]



annotation=data.frame(
  Phylum = as.factor(filtered_VCs$Phylum),
  Cluster = as.factor(filtered_VCs$Cluster),
  Significance = as.factor(filtered_VCs$All_Sig))
# Define custom color palettes for each factor variable
#significance_colors <- c("Significant" = "blue", "Non-significant" = "red")
#cluster_colors <- rainbow(length(levels(filtered_VCs$Cluster)))
#phylum_colors <- rainbow(length(levels(filtered_VCs$Phylum)))

# Assign custom colors to the factor variables in the annotation data frame
#annotation$Significance <- factor(annotation$Significance, levels = levels(filtered_VCs$Significance), labels = significance_colors)
#annotation$Cluster <- factor(annotation$Cluster, levels = levels(filtered_VCs$Cluster), labels = cluster_colors)
#annotation$Phylum <- factor(annotation$Phylum, levels = levels(filtered_VCs$Phylum), labels = phylum_colors)

#annotation=data.frame(Significance = as.factor(filtered_VCs$All_Sig))

rownames(annotation)=rownames(correlation_matrix)

folder="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Clust_Harm/"
pdf(paste0(folder,'heatmaps_ALL_Cluster_EC4.pdf'),bg='white')
heatmap_result <- pheatmap(
  correlation_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation,
  annotation_legend = TRUE,
  show_rownames = F,
  show_colnames = F,
  legend = TRUE,
  fontsize_row = 7,
  fontsize_legend = 5,
  annotation_name_row = FALSE)
print(heatmap_result)
dev.off()

# Extract row dendrogram
row_dendrogram <- heatmap_result$tree_row
dendogram=data.frame(position=row_dendrogram$order,trait=row_dendrogram$labels)
dendogram$Genus=taxonomy$Genus[match(dendogram$trait,taxonomy$Species)]
dendogram$Family=taxonomy$Family[match(dendogram$trait,taxonomy$Species)]
dendogram$Phylum=taxonomy$Phylum[match(dendogram$trait,taxonomy$Species)]
#dendogram$Cluster=clusters$Cluster[match(dendogram$trait,clusters$trait)]
dendogram$All_Sig=filtered_VCs$All_Sig[match(dendogram$trait,filtered_VCs$trait)]





get_dendrogram_positions <- function(dendrogram) {
  if (is.leaf(dendrogram)) {
    return(list(name = dendrogram$order, position = dendrogram$height))
  } else {
    left <- get_dendrogram_positions(dendrogram$left)
    right <- get_dendrogram_positions(dendrogram$right)
    position <- (left$position + right$position) / 2
    return(list(name = c(left$name, right$name), position = position))
  }
}

# Get positions
row_positions <- get_dendrogram_positions(row_dendrogram)

# Create dataframe
row_dendrogram_df <- data.frame(
  row_name = row_positions$name,
  position = row_positions$position
)

# Print dataframe
print(row_dendrogram_df)




# Collapse species
matrix = t(correlation_matrix)
colnames(matrix)=filtered_VCs$All_Sig[match(colnames(matrix), filtered_VCs$trait)]
collapsed_matrix=t(matrix)

# Step 1: Find the indices of rows "cluster__11" and "cluster__12"
row_indices <- which(rownames(correlation_matrix) %in% c("cluster__11", "cluster__12"))
row_indices <- which(rownames(collapsed_matrix) == "Significant")
row_indices <- which(rownames(collapsed_matrix) != "Significant")

# Step 2: Compute the absolute differences with other rows

median_values_sig <- apply(collapsed_matrix[row_indices, ], 2, median)
median_values_notsig <- apply(collapsed_matrix[row_indices, ], 2, median)

differences_sig_vs_notsig <- abs(median_values_sig - median_values_notsig)

median_values_all <- apply(correlation_matrix, 2, median)
median_values_sig <- apply(correlation_matrix[row_indices[c(1,2)], ], 2, median)

differences_sig_vs_all <- abs(median_values_all - median_values_sig)



differences_11_vs_12 <- abs(correlation_matrix[row_indices[1], ] - correlation_matrix[row_indices[2], ])

dif_table=data.frame(EC=colnames(correlation_matrix),
                     differences_sig_vs_all=differences_sig_vs_all,
                     differences_11_vs_12=differences_11_vs_12,
                     corr_11=correlation_matrix[row_indices[1], ],
                     corr_12=correlation_matrix[row_indices[2], ])

dif_table=data.frame(EC=colnames(correlation_matrix),
                     differences_sig_vs_notsig=differences_sig_vs_notsig,
                     median_values_notsig=median_values_notsig,
                     median_values_sig=median_values_sig)






# Step 3: Identify column names where the difference is minimal
min_difference_cols <- colnames(correlation_matrix)[which(differences == min(differences))]
print(min_difference_cols)
# [1] "EC4__2.3.1.184" (acylhomoserine lactone synthase / AHL synthase)


# Step 3: Identify column names where the difference is maximal
max_difference_cols <- colnames(correlation_matrix)[which(differences == max(differences))]
print(max_difference_cols)
# [1] "EC4__2.3.1.184" (acylhomoserine lactone synthase / AHL synthase)