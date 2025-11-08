subset_ids=c("MI","NY","ALL")
ranks="Phylum,Class,Order,Family,Genus,Species"
ranks<-na.omit(c((sapply(strsplit(ranks, ","),"[",1)),
                 (sapply(strsplit(ranks, ","), "[",2)),
                 (sapply(strsplit(ranks, ","), "[",3)),
                 (sapply(strsplit(ranks, ","), "[",4)),
                 (sapply(strsplit(ranks, ","), "[",5)),
                 (sapply(strsplit(ranks, ","), "[",6))))
subset_id="ALL"
rank="Phylum"


if ("Species" %in% ranks){
  # 8. Check the importance of the taxonomic rank for heritability values on the Species level
  print("Checking the importance of the taxonomic rank for heritability values on the Species level")
  anova_total=data.frame()
  pdf('species_heritability_by_rank.pdf',bg='white')
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
    ## Plot the comparisons between the median heritability of all species on the rank with the own rank heritability
    for (test_rank in test_ranks){
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
        labs(x = "Top 10 traits with higher Species median heritability", y=paste0("Heritability (%) (Subset: ",subset_id," / Rank: ",test_rank,")"))
        print(myplot)
    }
  }
  dev.off()
  save(anova_total,
       file="anova_rank_effect_on_heritability.RData")
}



########


median_var_Ad$rank_herit=filtered_VCs_rank$var_Ad[match(median_var_Ad[[rank]],filtered_VCs_rank$trait)]
median_var_Ad_long <- pivot_longer(median_var_Ad, cols = c(var_Ad, rank_herit), names_to = "Variable")
colnames(median_var_Ad_long)[1]="trait"
median_var_Ad_long$Significance=filtered_VCs_rank$Significance[match(median_var_Ad_long$trait,filtered_VCs_rank$trait)]