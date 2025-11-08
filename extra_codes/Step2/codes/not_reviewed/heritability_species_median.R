## Plot the comparisons between the median heritability of all species on the rank with the own rank heritability
test_rank="rank"

median_var_Ad <- aggregate(var_Ad ~ ., data = all_VCs_full[, c("var_Ad", test_rank)], FUN = median)
median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
colnames(median_var_Ad)[1]="trait"
#if (length(median_var_Ad$trait)>10){
#  median_var_Ad=median_var_Ad[1:10,]
#}
order=median_var_Ad$trait
filtered_VCs_rank=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs_rank=filtered_VCs_rank[grepl(test_rank,filtered_VCs_rank$rank),]
filtered_VCs_rank=filtered_VCs_rank[filtered_VCs_rank$trait %in% median_var_Ad$trait,]
#filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (Bonferroni)" ] <- "**"
#filtered_VCs_rank$LBL[filtered_VCs_rank$Significance == "Significant (FDR)"] <- "*"
#filtered_VCs_rank$LBL[is.na(filtered_VCs_rank$LBL)] <- ""
filtered_VCs_rank$heritability=paste(round(filtered_VCs_rank$var_Ad*100,digit=2),filtered_VCs_rank$LBL,sep="")
filtered_VCs_rank$median=median_var_Ad$var_Ad[match(filtered_VCs_rank$trait,median_var_Ad$trait)]
filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad > filtered_VCs_rank$median] <- "Rank Heritability > Species Median"
filtered_VCs_rank$comp[filtered_VCs_rank$var_Ad < filtered_VCs_rank$median] <- "Rank Heritability < Species Median"
filtered_VCs_rank$comp[is.na(filtered_VCs_rank$comp)] <- ""
plot_data=filtered_VCs[,c(which(colnames(filtered_VCs)=="var_Ad"),which(colnames(filtered_VCs)==test_rank))]
colnames(plot_data)[which(colnames(plot_data)==test_rank)]="test_rank"
plot_data=plot_data[plot_data$test_rank %in% median_var_Ad$trait,]
# Plot
myplot<-ggplot(plot_data[plot_data$test_rank=="cluster__3" | plot_data$test_rank=="cluster__10" | plot_data$test_rank=="cluster__11",], aes(x = factor(test_rank, levels = order), y = var_Ad*100)) +
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
  labs(x = "Top 10 traits with higher species median heritability", y=paste0("Heritability (%) (Subset: ",subset_id," / Rank: ",test_rank,")"))
print(myplot)
