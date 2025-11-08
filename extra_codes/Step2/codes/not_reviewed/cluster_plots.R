input_taxonomy="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt"
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

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Cluster_Harmonized/Heritability/heritability.RData")

subset_id="ALL"
all_VCs_plot=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
#all_VCs_plot=all_VCs_full[grepl("Clusters",all_VCs_full$rank) | grepl("Genus",all_VCs_full$rank) | grepl("Species",all_VCs_full$rank),]
myplot<-ggplot(all_VCs_plot, aes(x=reorder(rank,100*var_Ad,FUN = median), y=var_Ad*100)) +
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
  #geom_signif(comparisons = list(c("Genus","Species"),c("Clusters","Species")),
  #            map_signif_level = TRUE)
print(myplot)

## Plot the comparisons between the median heritability of all ASVs on the rank with the own rank heritability


Cluster_UMAP_HDBSCAN_ALL <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Cluster_Harmonized/Matrix_Processing/Cluster_UMAP_HDBSCAN_ALL.csv")
clusters=Cluster_UMAP_HDBSCAN_ALL
colnames(clusters)[1]="Genome"
ALL <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Cluster_Harmonized/Matrix_Processing/filtered_tables/ALL.csv")
clusters$Genome=ALL$Genome_ID
clusters$Cluster=gsub("-1","noise",clusters$Cluster)
clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")

filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
filtered_VCs$Genome=taxonomy$Genome[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Cluster=clusters$Cluster[match(filtered_VCs$Genome,clusters$Genome)]
median_var_Ad <- aggregate(var_Ad ~ ., data = filtered_VCs[, c("var_Ad", "Cluster")], FUN = median)
median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
colnames(median_var_Ad)[1]="trait"
if (length(median_var_Ad$trait)>10){
  median_var_Ad=median_var_Ad[1:10,]
}
order=median_var_Ad$trait
filtered_VCs_cluster=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs_cluster=filtered_VCs_cluster[grepl("Cluster",filtered_VCs_cluster$rank),]
filtered_VCs_cluster=filtered_VCs_cluster[filtered_VCs_cluster$trait %in% median_var_Ad$trait,]
filtered_VCs_cluster$LBL[filtered_VCs_cluster$Significance == "Significant (Bonferroni)" ] <- "**"
filtered_VCs_cluster$LBL[filtered_VCs_cluster$Significance == "Significant (FDR)"] <- "*"
filtered_VCs_cluster$LBL[is.na(filtered_VCs_cluster$LBL)] <- ""
filtered_VCs_cluster$heritability=paste(round(filtered_VCs_cluster$var_Ad*100,digit=2),filtered_VCs_cluster$LBL,sep="")
filtered_VCs_cluster$median=median_var_Ad$var_Ad[match(filtered_VCs_cluster$trait,median_var_Ad$trait)]
filtered_VCs_cluster$comp[filtered_VCs_cluster$var_Ad > filtered_VCs_cluster$median] <- "Cluster Heritability > Species Median"
filtered_VCs_cluster$comp[filtered_VCs_cluster$var_Ad < filtered_VCs_cluster$median] <- "Cluster Heritability < Species Median"
filtered_VCs_cluster$comp[is.na(filtered_VCs_cluster$comp)] <- ""
plot_data=filtered_VCs[,c(which(colnames(filtered_VCs)=="var_Ad"),which(colnames(filtered_VCs)=="Cluster"))]
plot_data=plot_data[plot_data$Cluster %in% median_var_Ad$trait,]
# Plot
myplot<-ggplot(plot_data, aes(x = factor(Cluster, levels = order), y = var_Ad*100)) +
  geom_boxplot() +
  geom_label(data = filtered_VCs_cluster,  
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
  scale_fill_manual(values = c("Cluster Heritability > Species Median" = "lightgreen", "Cluster Heritability < Species Median" = "lightcoral"),
                    name = "") +
  #geom_signif(comparisons = list(c("cluster__7","cluster__6"),c("cluster__6","cluster__8"),c("cluster__8","cluster__5")),
  #            map_signif_level = TRUE) +
  labs(x = "", y=paste0("Heritability (%) (Subset: ",subset_id))
print(myplot)

phylotree="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/phylotree_GTDB207.tree"
tree = ape::read.tree(phylotree)

for (subset_id in subset_ids){
  all_VCs_plot=filtered_VCs
  keep_taxa=all_VCs_plot$Genome
  any(is.na(keep_taxa))
  remove_taxa = setdiff(tree$tip.label, keep_taxa)
  any(is.na(remove_taxa))
  pruned_tree = drop.tip(tree, remove_taxa)
  egde_labels <- all_VCs_plot$Significance
  pruned_tree$egde.label <- egde_labels
  pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')
  significance_colors <- ifelse(all_VCs_plot$Cluster == "cluster__6", "red",
                                ifelse(all_VCs_plot$Cluster == "cluster__7", "blue","lightgrey"))
  significance_colors <- ifelse(all_VCs_plot$Genus == "g__Prevotella", "red",
                                ifelse(all_VCs_plot$Genus == "g__Bacteroides", "blue","lightgrey"))
  plot.phylo(pruned_tree, cex = 0.8, type = "unrooted", edge.color = significance_colors, show.tip.label = FALSE, show.node.label = FALSE, main = paste0("Subset: ",subset_id))
  #edge_labels=c("Cluster 6","Cluster 7")
  #significance_colors=c("red","blue")
  #legend("bottom", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
}