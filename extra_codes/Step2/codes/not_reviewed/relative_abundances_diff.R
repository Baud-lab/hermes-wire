subset_id="ALL"

total_family=c()
for (subset_id in subset_ids){
  family_matrix <- get(paste("filtered_prev_Genus",subset_id,sep="_"))
  total_family=unique(c(total_family,rownames(family_matrix)))
}
total_family=c(total_family,"Rare Taxa")
colours = metafolio::gg_color_hue(n =  length(total_family), c = 100, l = 65)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]


# Make into relative abundance
family_matrix <- apply(get(paste("filtered_prev_Genus",subset_id,sep="_")), 2, function(i) i/sum(i))
colSums(family_matrix)[1:5]
# Define a cutoff for rare taxa
maxabundances <- apply(family_matrix, 1, max)
# Meanwhile, transpose the count table for future wrangling.
family_matrix <- data.frame(t(family_matrix))
rare=which(maxabundances < 0.01)
if (length(rare)>1){
  family_matrix$`Rare Taxa` <- rowSums(family_matrix[,maxabundances < 0.01], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
  family_matrix = family_matrix[,c(maxabundances > 0.01, T) ] #`T` to include the `Rare Taxa`.
} else {
  colnames(family_matrix)[which(colnames(family_matrix)==names(maxabundances)[rare])]="Rare Taxa"
}
motch = match(rownames(family_matrix), metadata[[sample_identifier]])
any(is.na(motch))
this_metadata = metadata[motch,]
covariate="Study"
family_matrix[[covariate]]=this_metadata[[covariate]]
family_matrix$Sample=rownames(family_matrix)
# Wrangle the data to long format for easy plotting
formula <- c("Sample")
formula <- c(formula, covariate)
barlong = family_matrix %>%
  pivot_longer(!formula, names_to = c("family"), values_to = "value") %>%
  mutate(family = str_replace(family, ".*_or_", ""))
# Calculate mean values across covariates
grouped_matrix <- barlong %>%
  group_by(family, !!sym(covariate)) %>%
  summarise(mean_value_per_family = mean(value))
mi_data <- grouped_matrix %>%
  filter(Study == "MI")

ny_data <- grouped_matrix %>%
  filter(Study == "NY")

# Join the dataframes by family
merged_data <- merge(mi_data, ny_data, by = "family", suffixes = c("_MI", "_NY"))

# Calculate the difference in mean_value_per_family
merged_data <- merged_data %>%
  mutate(diff_mean = mean_value_per_family_NY - mean_value_per_family_MI) %>%
  select(family, diff_mean)

# If you want to filter out families with NA differences
merged_data <- merged_data[complete.cases(merged_data), ]
save(merged_data,grouped_matrix,file="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow/Cluster_Analyses/genus_abundance_differences.RData")





pdf('heritability_significant_vs_non_significant.pdf',bg='white')
for (subset_id in subset_ids){
  filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  filtered_VCs <- filtered_VCs %>% 
    mutate( All_Sig = case_when(
      Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Significant",
      TRUE ~ "Non-significant"))
  if (method=="16S"){
    if ("ASV" %in% ranks){
      print("Comparisons between significant and non-significant ASVs")
      #
      filtered_VCs=filtered_VCs[grepl("ASV",filtered_VCs$rank),]
      filtered_VCs$Genome_size=catalogue_stats$Size[match(filtered_VCs$trait,catalogue_stats$ASV)]
      filtered_VCs$Number_of_sequences=catalogue_stats$Number_of_sequences[match(filtered_VCs$trait,catalogue_stats$ASV)]
      filtered_VCs$Mean_gene_size=catalogue_stats$mean[match(filtered_VCs$trait,catalogue_stats$ASV)]
      filtered_VCs$family=taxonomy$family[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Genome)]
      rank="ASV"
    }
  } else {
      if ("Species" %in% ranks){
        print("Comparisons between significant and non-significant Species")
        filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
        filtered_VCs$Genome_size=catalogue_stats$Size[match(filtered_VCs$trait,catalogue_stats$Species)]
        filtered_VCs$Number_of_sequences=catalogue_stats$Number_of_sequences[match(filtered_VCs$trait,catalogue_stats$Species)]
        filtered_VCs$Mean_gene_size=catalogue_stats$mean[match(filtered_VCs$trait,catalogue_stats$Species)]
        filtered_VCs$family=taxonomy$family[match(filtered_VCs$trait,taxonomy$Species)]
        filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
        filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
        filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
        filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
        rank="Species"
      }
    }
  
  # Calculate the number of significant (FDR) in each rank
  stats_family=count(filtered_VCs,All_Sig,family)
  stats_class=count(filtered_VCs,All_Sig,Class)
  stats_order=count(filtered_VCs,All_Sig,Order)
  stats_family=count(filtered_VCs,All_Sig,Family)
  stats_genus=count(filtered_VCs,All_Sig,Genus)
  save(stats_family,
       stats_class,
       stats_order,
       stats_family,
       stats_genus, file="statistics_per_rank.RData")
  
  if (length(filtered_VCs[filtered_VCs$All_Sig=="Significant"])>0 ){
    
    # Add relative abundance and prevalence
    matrix=get(paste("filtered_prev",rank,subset_id,sep="_"))
    ab_vs_prev=data.frame(trait=rownames(matrix),
                          abundance=100*apply(sweep(matrix,2,colSums(matrix),"/"), 1, mean, na.rm=TRUE),
                          prevalence=100*rowSums(matrix>0)/ncol(matrix))
    filtered_VCs$Abundance=ab_vs_prev$abundance[match(filtered_VCs$trait,ab_vs_prev$trait)]
    filtered_VCs$Prevalence=ab_vs_prev$prevalence[match(filtered_VCs$trait,ab_vs_prev$trait)]
    # Make the plots
    ## Comparing all types of Significance
    order=c("Non-significant traits","Significant (FDR)","Significant (Bonferroni)")
    ### Relative Abundance
    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=Abundance)) +
      geom_boxplot() +
      xlab("") +
      ylab("Relative Abundance (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
                  map_signif_level = TRUE)
    print(myplot)
  }
}
dev.off()

























load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow/Cluster_Analyses/filtered_prev.RData")
subset_id="ALL"
out="/users/abaud/fmorillo/cluster_test/"

filtered=filtered_prev_objs[[18]]
matrix=data.frame(Genome_ID=rownames(filtered))
matrix=cbind(matrix,as.data.frame(filtered))
write.table(matrix,file=paste0(subset_id,"_PREV.csv"),sep=",",row.names=F,quote=F)



matrix=data.frame(Genome_ID=rownames(residuals_qned_counts))
matrix=cbind(matrix,as.data.frame(residuals_qned_counts))
write.table(matrix,file=paste0(out,subset_id,"_PREV.csv"),sep=",",row.names=F,quote=F)
