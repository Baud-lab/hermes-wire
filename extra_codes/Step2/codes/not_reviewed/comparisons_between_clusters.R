
genus="g__Bacteroides"
genus="g__Prevotella"
genera=c("g__Bacteroides","g__Prevotella")

for (genus in genera){
  prevotella=filtered_VCs[filtered_VCs$Genus==genus,]
  unique(prevotella$Cluster)
  if (genus=="g__Bacteroides"){
    cluster="cluster__3"
  } else {cluster="cluster__11"}
  clust11=prevotella[prevotella$Cluster==cluster,]
  clust11_size=nrow(clust11)
  clust10=prevotella[prevotella$Cluster=="cluster__10",]
  clust10_size=nrow(clust10)
  
  ####
  #clust11=clust10[clust10$Hub=="Hubs",]
  #clust11_size=nrow(clust11)
  #clust10=clust10[clust10$Hub=="Others",]
  #clust10_size=nrow(clust10)
  ####
  
  this_filtered_VCs2=rbind(clust11,clust10)
  comparison="heritability"
  
  # Plot the histogram
  myplot=ggplot(this_filtered_VCs2, aes(x = Heritability, fill = Cluster)) +
    geom_density(position = "identity", alpha = 0.5) +
    labs(title = paste0("Heritability distribution of ",gsub("g__","",genus)),
         x = "Heritability (%)",
         y = "Frequency") +
    scale_fill_manual(values = c("blue", "red")) +
    theme_minimal()
  plot(myplot)
}
dev.off()


## Plot the phylogenetic tree
#if (comparison=="heritability"){
#  keep_taxa=taxonomy$Genome[match(this_filtered_VCs2$trait,taxonomy$Species)]
#} else {
#  keep_taxa=taxonomy$Genome[match(rownames(chr10_filt),taxonomy$Species)]
#}
#any(is.na(keep_taxa))
#remove_taxa = setdiff(tree$tip.label, keep_taxa)
#any(is.na(remove_taxa))
#pruned_tree = drop.tip(tree, remove_taxa)
#tips=taxonomy$Species[match(pruned_tree$tip.label,taxonomy$Genome)]
#tips=gsub("s__","",tips)
#pruned_tree$tip.label <- tips
#pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
#if (comparison=="heritability"){
#  order=this_filtered_VCs2[match(pruned_tree$tip.label,gsub("s__","",this_filtered_VCs2$trait)),]
#  significance_colors <- ifelse(order$Cluster == "cluster__11", "red",
#                                ifelse(order$Cluster == "cluster__10", "blue", "darkgray"))
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", tip.color = significance_colors, 
#             show.tip.label = TRUE, show.node.label = FALSE, main = "")
#  edge_labels=c("Cluster 11","Cluster 10")
#  significance_colors=c("red","blue")
#  legend("bottom", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
#} else {
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", show.tip.label = TRUE, show.node.label = FALSE, main = "")
#}

func="CAZy"
func="eggNOG_OGs"
func="KEGG_Module"
func="KEGG_ko"
# Filter function
matrix=get(paste("func_matrix",func,sep="_"))
colnames(matrix)=taxonomy$Species[match(colnames(matrix),taxonomy$Genome)]
matrix=t(matrix)
if (comparison=="heritability"){
  matrix <- matrix[match(this_filtered_VCs2$trait,rownames(matrix)),]
} else {
  matrix <- matrix[match(rownames(chr10_filt),rownames(matrix)),]
}
# Filter out non-mapped functions
column_sums <- colSums(matrix)
non_zero_columns <- column_sums != 0
matrix <- matrix[, non_zero_columns]
# Filter out functions with no variation
is_constant_column <- function(col) {
  length(unique(col)) == 1
}
non_constant_columns <- apply(matrix, 2, function(col) !is_constant_column(col))
matrix <- matrix[, non_constant_columns]



# Define Presence/Absence
presence_absence <- ifelse(matrix > 0, 1, 0)
# Filter out functions with no variation
non_constant_columns <- apply(presence_absence, 2, function(col) !is_constant_column(col))
presence_absence <- presence_absence[, non_constant_columns]

matrix_filt=as.data.frame(t(matrix))
matrix_filt$median_high=apply(matrix_filt[, 1:(ncol(matrix_filt)-clust10_size)], 1, median)
matrix_filt$median_low=apply(matrix_filt[, (ncol(matrix_filt)-clust10_size+1):ncol(matrix_filt)], 1, median)
matrix_filt$diff=matrix_filt$median_low-matrix_filt$median_high
matrix_filt$func=rownames(matrix_filt)
matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = FALSE), ]
top=matrix_filt[1:10,]
#matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = TRUE), ]
#bottom=matrix_filt[1:5,]
#matrix_filt=rbind(top,bottom)
matrix_filt=top
matrix=matrix[,colnames(matrix) %in% matrix_filt$func, drop = FALSE]

presence_absence_filt=as.data.frame(t(presence_absence))
presence_absence_anova=as.data.frame(presence_absence)
differences=data.frame(counts_high=rowSums(presence_absence_filt[,1:(ncol(presence_absence_filt)-clust10_size)]),counts_low=rowSums(presence_absence_filt[,(ncol(presence_absence_filt)-clust10_size+1):ncol(presence_absence_filt)]))
presence_absence_filt=differences[(differences$counts_high==clust11_size & differences$counts_low==0)
                                  | (differences$counts_high==0 & differences$counts_low==clust10_size)
                                  ,]
if (length(rownames(presence_absence_filt))!=0){
  presence_absence=presence_absence[,colnames(presence_absence) %in% rownames(presence_absence_filt), drop = FALSE]
}

# Make plots
if (func == "CAZy"){
  operator=TRUE
} else {operator=FALSE}
# Counts
if (genus!="mixed"){
  genus=sapply(strsplit(genus, "__"), "[",2)
}
heatmap_result <- pheatmap(
  matrix,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function counts (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  annotation_name_row = FALSE)
# Presence/Absence
binary_colors <- c("lightgrey", "black")
binary_breaks <- c(-0.1, 0.5, 1.1)
heatmap_result <- pheatmap(
  presence_absence,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function presence (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  color = binary_colors,
  breaks = binary_breaks,
  legend_breaks = c(0, 1),
  legend_labels = c("0", "1"),
  annotation_name_row = FALSE)


# Bacteroides/Prevotella vs CAG-485

prevotella=filtered_VCs[filtered_VCs$Genus=="g__Bacteroides" | filtered_VCs$Genus=="g__Prevotella" | filtered_VCs$Genus=="g__CAG-485",]
unique(prevotella$Genus)
clust11=prevotella[prevotella$Genus=="g__Prevotella" | prevotella$Genus=="g__Bacteroides",]
clust11_size=nrow(clust11)
clust10=prevotella[prevotella$Genus=="g__CAG-485",]
clust10_size=nrow(clust10)
this_filtered_VCs2=rbind(clust11,clust10)
comparison="heritability"

## Plot the phylogenetic tree
#if (comparison=="heritability"){
#  keep_taxa=taxonomy$Genome[match(this_filtered_VCs2$trait,taxonomy$Species)]
#} else {
#  keep_taxa=taxonomy$Genome[match(rownames(chr10_filt),taxonomy$Species)]
#}
#any(is.na(keep_taxa))
#remove_taxa = setdiff(tree$tip.label, keep_taxa)
#any(is.na(remove_taxa))
#pruned_tree = drop.tip(tree, remove_taxa)
#tips=taxonomy$Species[match(pruned_tree$tip.label,taxonomy$Genome)]
#tips=gsub("s__","",tips)
#pruned_tree$tip.label <- tips
#pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
#if (comparison=="heritability"){
#  order=this_filtered_VCs2[match(pruned_tree$tip.label,gsub("s__","",this_filtered_VCs2$trait)),]
#  significance_colors <- ifelse(order$Genus == "g__Bacteroides", "red",
#                                ifelse(order$Genus == "g__Prevotella", "blue", "green"))
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", tip.color = significance_colors, 
#             show.tip.label = TRUE, show.node.label = FALSE, main = "")
#  edge_labels=c("Bacterodes","Prevotella", "CAG-485")
#  significance_colors=c("red","blue","green")
#  legend("bottom", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
#} else {
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", show.tip.label = TRUE, show.node.label = FALSE, main = "")
#}

func="eggNOG_OGs"
# Filter function
matrix=get(paste("func_matrix",func,sep="_"))
colnames(matrix)=taxonomy$Species[match(colnames(matrix),taxonomy$Genome)]
matrix=t(matrix)
if (comparison=="heritability"){
  matrix <- matrix[match(this_filtered_VCs2$trait,rownames(matrix)),]
} else {
  matrix <- matrix[match(rownames(chr10_filt),rownames(matrix)),]
}
# Filter out non-mapped functions
column_sums <- colSums(matrix)
non_zero_columns <- column_sums != 0
matrix <- matrix[, non_zero_columns]
# Filter out functions with no variation
is_constant_column <- function(col) {
  length(unique(col)) == 1
}
non_constant_columns <- apply(matrix, 2, function(col) !is_constant_column(col))
matrix <- matrix[, non_constant_columns]
# Define Presence/Absence
presence_absence <- ifelse(matrix > 0, 1, 0)
# Filter out functions with no variation
non_constant_columns <- apply(presence_absence, 2, function(col) !is_constant_column(col))
presence_absence <- presence_absence[, non_constant_columns]

matrix_filt=as.data.frame(t(matrix))
matrix_filt$median_high=apply(matrix_filt[, 1:(ncol(matrix_filt)-clust10_size)], 1, median)
matrix_filt$median_low=apply(matrix_filt[, (ncol(matrix_filt)-clust10_size+1):ncol(matrix_filt)], 1, median)
matrix_filt$diff=matrix_filt$median_low-matrix_filt$median_high
matrix_filt$func=rownames(matrix_filt)
matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = FALSE), ]
top=matrix_filt[1:10,]
#matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = TRUE), ]
#bottom=matrix_filt[1:5,]
#matrix_filt=rbind(top,bottom)
matrix_filt=top
matrix=matrix[,colnames(matrix) %in% matrix_filt$func, drop = FALSE]

presence_absence_filt=as.data.frame(t(presence_absence))
presence_absence_anova=as.data.frame(presence_absence)
differences=data.frame(counts_high=rowSums(presence_absence_filt[,1:(ncol(presence_absence_filt)-clust10_size)]),counts_low=rowSums(presence_absence_filt[,(ncol(presence_absence_filt)-clust10_size+1):ncol(presence_absence_filt)]))
presence_absence_filt=differences[(differences$counts_high==clust11_size & differences$counts_low==0)
                                  | (differences$counts_high==0 & differences$counts_low==clust10_size)
                                  ,]
if (length(rownames(presence_absence_filt))!=0){
  presence_absence=presence_absence[,colnames(presence_absence) %in% rownames(presence_absence_filt), drop = FALSE]
}

# Make plots
if (func == "CAZy"){
  operator=TRUE
} else {operator=FALSE}
# Counts
if (genus!="mixed"){
  genus=sapply(strsplit(genus, "__"), "[",2)
}
heatmap_result <- pheatmap(
  matrix,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function counts (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  annotation_name_row = FALSE)
# Presence/Absence
binary_colors <- c("lightgrey", "black")
binary_breaks <- c(-0.1, 0.5, 1.1)
heatmap_result <- pheatmap(
  presence_absence,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function presence (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  color = binary_colors,
  breaks = binary_breaks,
  legend_breaks = c(0, 1),
  legend_labels = c("0", "1"),
  annotation_name_row = FALSE)

# Bacteroides

prevotella=filtered_VCs[filtered_VCs$Genus=="g__Bacteroides",]
unique(prevotella$Cluster)
clust11=prevotella[prevotella$Cluster=="cluster__3",]
clust11_size=nrow(clust11)
clust10=prevotella[prevotella$Cluster=="cluster__10",]
clust10_size=nrow(clust10)
this_filtered_VCs2=rbind(clust11,clust10)
comparison="heritability"

## Plot the phylogenetic tree
#if (comparison=="heritability"){
#  keep_taxa=taxonomy$Genome[match(this_filtered_VCs2$trait,taxonomy$Species)]
#} else {
#  keep_taxa=taxonomy$Genome[match(rownames(chr10_filt),taxonomy$Species)]
#}
#any(is.na(keep_taxa))
#remove_taxa = setdiff(tree$tip.label, keep_taxa)
#any(is.na(remove_taxa))
#pruned_tree = drop.tip(tree, remove_taxa)
#tips=taxonomy$Species[match(pruned_tree$tip.label,taxonomy$Genome)]
#tips=gsub("s__","",tips)
#pruned_tree$tip.label <- tips
#pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
#if (comparison=="heritability"){
#  order=this_filtered_VCs2[match(pruned_tree$tip.label,gsub("s__","",this_filtered_VCs2$trait)),]
#  significance_colors <- ifelse(order$Cluster == "cluster__3", "red",
#                                ifelse(order$Cluster == "cluster__10", "blue", "darkgray"))
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", tip.color = significance_colors, 
#             show.tip.label = TRUE, show.node.label = FALSE, main = "")
#  edge_labels=c("Cluster 3","Cluster 10")
#  significance_colors=c("red","blue")
#  legend("bottom", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
#} else {
#  plot.phylo(pruned_tree, cex = 0.5, type = "phylogram", show.tip.label = TRUE, show.node.label = FALSE, main = "")
#}

func="eggNOG_OGs"
# Filter function
matrix=get(paste("func_matrix",func,sep="_"))
colnames(matrix)=taxonomy$Species[match(colnames(matrix),taxonomy$Genome)]
matrix=t(matrix)
if (comparison=="heritability"){
  matrix <- matrix[match(this_filtered_VCs2$trait,rownames(matrix)),]
} else {
  matrix <- matrix[match(rownames(chr10_filt),rownames(matrix)),]
}
# Filter out non-mapped functions
column_sums <- colSums(matrix)
non_zero_columns <- column_sums != 0
matrix <- matrix[, non_zero_columns]
# Filter out functions with no variation
is_constant_column <- function(col) {
  length(unique(col)) == 1
}
non_constant_columns <- apply(matrix, 2, function(col) !is_constant_column(col))
matrix <- matrix[, non_constant_columns]
# Define Presence/Absence
presence_absence <- ifelse(matrix > 0, 1, 0)
# Filter out functions with no variation
non_constant_columns <- apply(presence_absence, 2, function(col) !is_constant_column(col))
presence_absence <- presence_absence[, non_constant_columns]

matrix_filt=as.data.frame(t(matrix))
matrix_filt$median_high=apply(matrix_filt[, 1:(ncol(matrix_filt)-clust10_size)], 1, median)
matrix_filt$median_low=apply(matrix_filt[, (ncol(matrix_filt)-clust10_size+1):ncol(matrix_filt)], 1, median)
matrix_filt$diff=matrix_filt$median_low-matrix_filt$median_high
matrix_filt$func=rownames(matrix_filt)
matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = FALSE), ]
top=matrix_filt[1:10,]
#matrix_filt <- matrix_filt[order(matrix_filt$diff, decreasing = TRUE), ]
#bottom=matrix_filt[1:5,]
#matrix_filt=rbind(top,bottom)
matrix_filt=top
matrix=matrix[,colnames(matrix) %in% matrix_filt$func, drop = FALSE]

presence_absence_filt=as.data.frame(t(presence_absence))
presence_absence_anova=as.data.frame(presence_absence)
differences=data.frame(counts_high=rowSums(presence_absence_filt[,1:(ncol(presence_absence_filt)-clust10_size)]),counts_low=rowSums(presence_absence_filt[,(ncol(presence_absence_filt)-clust10_size+1):ncol(presence_absence_filt)]))
presence_absence_filt=differences[(differences$counts_high==clust11_size & differences$counts_low==0)
                                  | (differences$counts_high==0 & differences$counts_low==clust11_size)
                                  ,]
if (length(rownames(presence_absence_filt))!=0){
  presence_absence=presence_absence[,colnames(presence_absence) %in% rownames(presence_absence_filt), drop = FALSE]
}

# Make plots
if (func == "CAZy"){
  operator=TRUE
} else {operator=FALSE}
# Counts
if (genus!="mixed"){
  genus=sapply(strsplit(genus, "__"), "[",2)
}
heatmap_result <- pheatmap(
  matrix,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function counts (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  annotation_name_row = FALSE)
# Presence/Absence
binary_colors <- c("lightgrey", "black")
binary_breaks <- c(-0.1, 0.5, 1.1)
heatmap_result <- pheatmap(
  presence_absence,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  #main=paste0("Function presence (Functions: ",func," / Genus: ",genus,")"),
  legend = TRUE,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  fontsize_legend = 5,
  #angle_col = "45",
  color = binary_colors,
  breaks = binary_breaks,
  legend_breaks = c(0, 1),
  legend_labels = c("0", "1"),
  annotation_name_row = FALSE)