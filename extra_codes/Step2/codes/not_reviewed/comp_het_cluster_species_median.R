invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}

subset_ids=c("MI","NY","ALL")
subset_id="ALL"

filtered_noprev <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/filtered_noprev_ALL.txt")
UMAP_HDBSCAN <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/Output_UMAP_HDBSCAN_ALL.csv")
colnames(UMAP_HDBSCAN)[1]="Genome"
UMAP_HDBSCAN$Genome=filtered_noprev$Genome_ID

load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/new_29-01-24/Heritability/heritability.RData")
species=all_VCs_full
load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/09_02_24/Output_Shallow_cluster/Heritability/heritability.RData")
clusters=all_VCs_full

taxonomy=read.table("/users/abaud/data/secondary/indexes/gtdb/genes/gtdb_ids_taxon.txt", sep="\t")
colnames(taxonomy)=c("genome","taxonomy")
taxonomy$taxonomy<-gsub(" ","_",as.character(taxonomy$taxonomy))
taxonomy$phylum<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",2)
taxonomy$class<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",3)
taxonomy$order<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",4)
taxonomy$family<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",5)
taxonomy$genus<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",6)
taxonomy$species<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",7)

anova_total=data.frame()
pdf('species_heritability_by_cluster.pdf',bg='white')
for (subset_id in subset_ids){
  ## Filter Species by study and assign higher taxonomic ranks
  filtered_VCs=species[grepl(subset_id,species$subset_id),]
  filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
  filtered_VCs$Genome=taxonomy$genome[match(filtered_VCs$trait,taxonomy$species)]
  filtered_VCs$Cluster=UMAP_HDBSCAN$Cluster[match(filtered_VCs$Genome,UMAP_HDBSCAN$Genome)]
  filtered_VCs$Cluster=paste0("cluster__",filtered_VCs$Cluster)
  filtered_VCs=filtered_VCs[!grepl("cluster__-1",filtered_VCs$Cluster) & !grepl("cluster__NA",filtered_VCs$Cluster),]
  ## ANOVA test on clusters
  filtered_VCs$norm_heritability = invrank(filtered_VCs$var_Ad)
  anova.herit <- summary(aov(norm_heritability ~ Cluster, data = filtered_VCs))
  anova.herit=anova.herit[[1]]
  anova.herit$R2=anova.herit$`Sum Sq`/sum(anova.herit$`Sum Sq`)
  anova.herit$subset_id=subset_id
  anova_total=rbind(anova_total,anova.herit)
  ## Plot the comparisons between the median heritability of all species on the rank with the own rank heritability
    median_var_Ad <- aggregate(var_Ad ~ ., data = filtered_VCs[, c("var_Ad", "Cluster")], FUN = median)
    colnames(median_var_Ad)[1]="trait"
    median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
    order=median_var_Ad$trait
    filtered_clusters=clusters[grepl(subset_id,clusters$subset_id),]
    filtered_clusters=filtered_clusters[filtered_clusters$trait %in% median_var_Ad$trait,]
    filtered_clusters$LBL[filtered_clusters$Significance == "Significant (Bonferroni)" ] <- "**"
    filtered_clusters$LBL[filtered_clusters$Significance == "Significant (FDR)"] <- "*"
    filtered_clusters$LBL[is.na(filtered_clusters$LBL)] <- ""
    filtered_clusters$heritability=paste(round(filtered_clusters$var_Ad*100,digit=2),filtered_clusters$LBL,sep="")
    filtered_clusters$median=median_var_Ad$var_Ad[match(filtered_clusters$trait,median_var_Ad$trait)]
    filtered_clusters$comp[filtered_clusters$var_Ad > filtered_clusters$median] <- "Cluster Heritability > Species Median"
    filtered_clusters$comp[filtered_clusters$var_Ad < filtered_clusters$median] <- "Cluster Heritability < Species Median"
    filtered_clusters$comp[is.na(filtered_clusters$comp)] <- ""
    # Plot
    myplot<-ggplot(filtered_VCs, aes(x = factor(Cluster, levels = order), y = var_Ad*100)) +
      geom_boxplot() +
      geom_label(data = filtered_clusters,  
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
      labs(x = "Clusters with higher Species median heritability", y=paste0("Heritability (%) (Subset: ",subset_id,")"))
    print(myplot)
}
dev.off()
save(anova_total,
     file="anova_cluster_effect_on_heritability.RData")