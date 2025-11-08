##h5="/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/Rnor6/taxonomic_profile/gtdb_profiles/08_02_24/Output_Shallow/Matrix_Processing/exp_DGA_7_Dec_2021.h5"
#h5="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/P50_rats_Rn7.h5"
#h5f = H5Fopen(h5)
#
#
## Read count tables after CLR transformation
#residuals="~/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Matrix_Processing/residuals_qned_counts.RData"
#load(residuals)
#residuals_MI=residuals_qned_counts_objs[[1]]
#residuals_NY=residuals_qned_counts_objs[[2]]
#residuals_ALL=residuals_qned_counts_objs[[3]]
#colnames(residuals_MI)=sapply(strsplit(colnames(residuals_MI), "_"),"[",1)
#colnames(residuals_NY)=sapply(strsplit(colnames(residuals_NY), "_"),"[",1)
#colnames(residuals_ALL)=sapply(strsplit(colnames(residuals_ALL), "_"),"[",1)
#samples=unique(c(colnames(residuals_ALL),colnames(residuals_MI),colnames(residuals_NY)))
#
#grm=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/","/matrix"))
#ids=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/row_header/","/sample_ID"))
#H5Fclose(h5f)
#rownames(grm)=ids
#colnames(grm)=ids
#grm[1:3,1:3]
#grm=grm[rownames(grm) %in% samples,colnames(grm) %in% samples]
#
#
#
#dist_matrix <- 1 - grm
#
## Step 4: Create a distance object
## Use `as.dist` to convert to a distance object
#distance_obj <- as.dist(dist_matrix)
#
## Step 5: Perform PCoA using the `pcoa` function
#pcoa_results <- pcoa(distance_obj)
#
## Extract coordinates
#coordinates <- as.data.frame(pcoa_results$vectors)
#
#metadata_Shallow <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")
#View(metadata_Shallow)
#rownames(metadata_Shallow)=metadata_Shallow$RFID
#this_metadata=metadata_Shallow[rownames(metadata_Shallow) %in% rownames(coordinates),]
#
#myplot<-ggplot(coordinates, aes(x=Axis.1, y=Axis.2, colour = as.factor(this_metadata$Study))) +
#  geom_point() + 
#  xlab("PCo1") +
#  ylab("PCo2") +
#  theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5, size = 17),
#        axis.line = element_line(),
#        plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 13),
#        axis.title = element_text(size = 16),
#        legend.text=element_text(size=13),
#        legend.title = element_blank(),
#        legend.position="bottom")
#print(myplot)
#
#permanova.beta <- adonis2(distance_obj ~ Study, data = this_metadata, permutations = 999)
#
#
###################
#
#relatedeness=data.frame(samples=rownames(grm),distance=pcoa_results$values$Eigenvalues)
#
#load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Enterotypes/enterotypes.RData")
#enterotypes=enterotypes_samples_objs[[3]]
#enterotypes$sample=sapply(strsplit(enterotypes$sample, "_"), "[",1)
#rownames(enterotypes)=enterotypes$sample
#enterotypes$relatedeness=relatedeness$distance[match(enterotypes$sample,relatedeness$samples)]
#
##hist(enterotypes$relatedeness)
##invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}
##enterotypes$relatedeness = invrank(enterotypes$relatedeness)
##t_test_result <- t.test(relatedeness ~ Enterotype, data = enterotypes)
#
#dist_matrix <- 1 - grm
#distance_obj <- as.dist(dist_matrix)
#perm_result <- adonis2(distance_obj ~ Enterotype, data = enterotypes)
#
##enterotypes$enterotype_name=gsub("g__","",enterotypes$enterotype_name)
#
#myplot<-ggplot(enterotypes, aes(x=Enterotype, y=relatedeness)) +
#  geom_boxplot() +
#  xlab("") +
#  ylab("Genetic Distance (Eigenvalues)") +
#  #ggtitle(paste0("Subset: ",subset_id)) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 10),
#        axis.title = element_text(size = 15),
#        legend.text=element_text(size=10),
#        plot.title = element_text(hjust = 0.5))+
#  scale_x_discrete(guide = guide_axis(angle = 60))
#print(myplot)



###################


# Load required libraries
#library(rhdf5, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(vegan)     # For Mantel test and distance matrix handling
library(Matrix)    # For matrix operations (if needed)

# ------------------------------------------
# STEP 1: Load your GRM and Aitchison distance matrices
# ------------------------------------------
# Replace these with your actual file paths or objects
# Example if matrices are stored as text files:
# Read count tables after CLR transformation
#residuals="~/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Matrix_Processing/residuals_qned_counts.RData"
#load(residuals)
#residuals_ALL=residuals_qned_counts_objs[[3]]
#colnames(residuals_ALL)=sapply(strsplit(colnames(residuals_ALL), "_"),"[",1)
#samples=unique(c(colnames(residuals_ALL)))

#h5="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/P50_rats_Rn7.h5"
#h5f = H5Fopen(h5)
#grm=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/","/matrix"))
#ids=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/row_header/","/sample_ID"))
#H5Fclose(h5f)
#rownames(grm)=ids
#colnames(grm)=ids
##grm[1:3,1:3]
#grm=grm[rownames(grm) %in% samples,colnames(grm) %in% samples]
#save(grm,file="/users/abaud/fmorillo/paper_figures/grm.RData")
load("/users/abaud/fmorillo/paper_figures/grm.RData")
grm=na.omit(grm)

#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Diversity/beta_diversity.RData")
#matrix=as.matrix(gp.ilr_objs[[6]])
#save(gp.dist,file="/users/abaud/fmorillo/paper_figures/beta.RData")

load("/users/abaud/fmorillo/paper_figures/beta.RData")
aitchison_dist <- as.matrix(gp.dist)
colnames(aitchison_dist)=sapply(strsplit(colnames(aitchison_dist), "_"),"[",1)
rownames(aitchison_dist)=sapply(strsplit(rownames(aitchison_dist), "_"),"[",1)
aitchison_dist=na.omit(aitchison_dist)


#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Diversity/residuals_qned_counts_beta.RData")
#number=3
#subset_id="ALL"
#if (subset_id=="MI"){
#  number=number-2
#} else {
#  if (subset_id=="NY"){
#    number=number-1
#  }
#}
#beta=residuals_qned_counts_beta_objs[[number]]
#genome.pc <- beta[2,]
#aitchison_dist <- genome.pc - min(genome.pc, na.rm = TRUE)
#names(aitchison_dist)=sapply(strsplit(names(aitchison_dist), "_"),"[",1)

# Step 1: Find common individuals between GRM and Aitchison matrix
common_ids <- intersect(rownames(grm), colnames(aitchison_dist))
grm <- grm[common_ids, common_ids]
# Step 2: Subset and reorder Aitchison matrix to match GRM
#aitchison_dist <- aitchison_dist[colnames(aitchison_dist) %in% common_ids]
aitchison_dist <- aitchison_dist[common_ids, common_ids]

grm [1:5,1:5]
aitchison_dist[1:5,1:5]

# Ensure matrices are in the same order (same individuals, same row/column order)
# If necessary, reorder matrices here (e.g., using match() or rownames)

# ------------------------------------------
# STEP 2: Convert GRM to a distance matrix
# ------------------------------------------
# GRM is a *similarity* matrix; convert it to a *distance* matrix for comparison.
# Common method: Distance = 1 - GRM (or sqrt(1 - GRM) for Euclidean-like properties)
grm_dist <- as.dist(1 - grm)  # Use as.dist() to discard redundant upper triangle

# Convert Aitchison distance to a dist object
aitchison_dist <- as.dist(aitchison_dist)

# ------------------------------------------
# STEP 3: Perform Mantel test
# ------------------------------------------
# Test correlation between GRM-derived distance and Aitchison distance
mantel_result <- mantel(
  grm_dist, 
  aitchison_dist, 
  method = "spearman",   # Use "spearman" for rank-based correlation
  permutations = 9999   # Number of permutations for p-value
)

# Print results
print(mantel_result)
save(mantel_result,file="/users/abaud/fmorillo/paper_figures/mantel_result.RData")

# ------------------------------------------
# STEP 4 (Optional): Visualization
# ------------------------------------------
# Scatterplot of matrix values (upper triangles only)
# Extract upper triangles (to avoid redundancy)

grm_values <- grm_dist[upper.tri(grm_dist)]
aitchison_values <- aitchison_dist[upper.tri(aitchison_dist)]

pdf("/users/abaud/fmorillo/paper_figures/grm_pcoa.pdf")
plot(
  grm_values, 
  aitchison_values,
  pch = 16,
  col = "blue",
  xlab = "GRM-derived distance (1 - GRM)",
  ylab = "Aitchison distance",
  main = "Correlation between GRM and Aitchison Distance"
)
abline(lm(aitchison_values ~ grm_values), col = "red")
dev.off()

# Perform PCoA on Aitchison distance matrix
pcoa_result <- cmdscale(aitchison_dist, k = 2, eig = TRUE)
# Extract PCoA scores (sample coordinates in PCoA space)
pcoa_scores <- pcoa_result$points  # Samples × PCoA axes
#matrix <- gp.ilr_objs[[6]]
colnames(residuals_qned_counts)=sapply(strsplit(colnames(residuals_qned_counts), "_"),"[",1)
#rownames(residuals_qned_counts)=sapply(strsplit(rownames(residuals_qned_counts), "_"),"[",1)
residuals_qned_counts=na.omit(residuals_qned_counts)
residuals_qned_counts_filt <- residuals_qned_counts[, common_ids]
projected_loadings <- cor(t(residuals_qned_counts_filt), pcoa_scores)

# Contributions of taxa to PCoA1
pcoa1_loadings <- projected_loadings[, 1]  # First column corresponds to PCoA1

# Sort taxa by their absolute contributions to PCoA1
sorted_loadings <- sort(abs(pcoa1_loadings), decreasing = TRUE)
top_taxa <- names(sorted_loadings)[1:10]  # Top 10 taxa contributing to PCoA1
print(top_taxa)
