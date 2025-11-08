### Enterotypes relative abundances per species
# Make into relative abundance
species_matrix <- apply(matrix, 2, function(i) i/sum(i))
colSums(species_matrix)[1:5]
# Convert species_matrix to a dataframe
species_df <- as.data.frame(species_matrix)
species_df <- tibble::rownames_to_column(species_df, var = "Species")
# Gather the Species matrix to long format
species_long <- species_df %>%
  pivot_longer(cols = -Species, names_to = "sample", values_to = "Abundance")
# Merge with ent_samples to get enterotype information
merged_df <- species_long %>%
  inner_join(ent_samples, by = "sample")
# Calculate mean abundance of each species per enterotype
mean_abundance <- merged_df %>%
  group_by(Enterotype, Species) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  ungroup()
mean_abundance$MeanAbundance=100*mean_abundance$MeanAbundance
###
reshaped_data_species <- mean_abundance %>%
  spread(key = Enterotype, value = MeanAbundance) %>%
  mutate(log_fold_change = log2(`enterotype__2` / `enterotype__1`))
reshaped_data_species$diff=reshaped_data_species$enterotype__1 - reshaped_data_species$enterotype__2
reshaped_data_species$Genus=taxonomy$Genus[match(reshaped_data_species$Species,taxonomy$Species)]
reshaped_data_species$Family=taxonomy$Family[match(reshaped_data_species$Species,taxonomy$Species)]
reshaped_data_species$Class=taxonomy$Class[match(reshaped_data_species$Species,taxonomy$Species)]
reshaped_data_species$Phylum=taxonomy$Phylum[match(reshaped_data_species$Species,taxonomy$Species)]
reshaped_data_species$Heritability=filtered_VCs$Heritability[match(reshaped_data_species$Species,filtered_VCs$trait)]
reshaped_data_species$qvalue=filtered_VCs$qvalue[match(reshaped_data_species$Species,filtered_VCs$trait)]

#plot(reshaped_data_species$Heritability,reshaped_data_species$diff)
plot(reshaped_data_species$Heritability[grepl("s__Prevotella",reshaped_data_species$Species)],reshaped_data_species$log_fold_change[grepl("s__Prevotella",reshaped_data_species$Species)])

###
# Select top 10 genera per enterotype
top10_genera <- mean_abundance %>%
  group_by(Enterotype) %>%
  top_n(20, MeanAbundance) %>%
  ungroup()
top10_genera$group=paste(top10_genera$Enterotype,top10_genera$Species,sep="_")
## Create a new column to indicate the direction for the bidirectional plot
#top10_genera <- top10_genera %>%
#  mutate(Direction = ifelse(Enterotype == unique(Enterotype)[1], -MeanAbundance, MeanAbundance))
# Order genera by decreasing mean abundance within each enterotype
top10_genera$Phylum=taxonomy$Phylum[match(top10_genera$Species,taxonomy$Species)]
top10_genera <- top10_genera %>%
  group_by(Enterotype) %>%
  arrange(Enterotype, desc(Phylum)) %>%
  arrange(Enterotype, desc(MeanAbundance)) %>%
  mutate(Species = factor(Species, levels = unique(Species))) %>%
  ungroup()
# Creating the bidirectional barplot with the adjusted theme
top10_genera$Enterotype=gsub("enterotype__","Enterotype ",top10_genera$Enterotype)
top10_genera$LBL[top10_genera$Species == "g__Prevotella" & top10_genera$Enterotype == "enterotype__2"] <- "*"
top10_genera$LBL[top10_genera$Species == "g__UBA3282" & top10_genera$Enterotype == "enterotype__1"] <- "*"
top10_genera$LBL[is.na(top10_genera$LBL)] <- ""
max=max(top10_genera$MeanAbundance)
myplot = ggplot(top10_genera, aes(x = Species, y = MeanAbundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(y = "Mean Mean relative abundance (%)", x = NULL, 
       title = "Top 20 genera by enterotype") +
  theme_minimal() +
  ylim(0, max(top10_genera$MeanAbundance) + 1) +  # Adjusted to use max from the data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    #legend.position = "bottom",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_grid(Enterotype ~ ., scales = "free_y", switch = "y") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_text(data = top10_genera,
            aes(x = Species, y = MeanAbundance,
                label = format(LBL, nsmall = 0, digits = 1, scientific = FALSE)),
            color = "black", vjust = -0.5, hjust = 0.5, size = 7)
plot(myplot)


reshaped_data_species$Heritability=filtered_VCs$Heritability[match(reshaped_data_species$Species,filtered_VCs$trait)]
genera=unique(reshaped_data_species$Genus)


correlations_total=data.frame()
for(genus in genera){
  print(genus)
  matrix=reshaped_data_species[reshaped_data_species$Genus == genus,]
  if (length(matrix$Species)>3){
    correlations = cor.test(matrix$Heritability,matrix$log_fold_change, method = "pearson")
    new_row=data.frame(trait=genus,correlations=correlations$estimate,pvalue=correlations$p.value,logp=-log10(correlations$p.value))
    correlations_total=rbind(correlations_total,new_row)
  }
}
correlations_total$pvalue=round(correlations_total$pvalue,digit=6)
correlations_total$bonf_pvalue = correlations_total$pvalue*dim(correlations_total)[1]
if (length(correlations_total$pvalue)>=10){
  correlations_total$qvalue = qvalue(correlations_total$pvalue)$qvalue
  correlations_total <- correlations_total %>% 
    mutate( Significance = case_when(
      qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
      qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
      TRUE ~ "Non-significant"))
  correlations_total$qvalue = round(correlations_total$qvalue, digits = 2)
} else {
  correlations_total <- correlations_total %>% 
    mutate( Significance = case_when(
      bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
      TRUE ~ "Non-significant"))
}
correlations_total$Direction <- ifelse(correlations_total$correlations > 0, "Positive", "Negative")
correlations_total = correlations_total[order(correlations_total$correlations, decreasing = F),]
correlations_total$Phylum=taxonomy$Phylum[match(correlations_total$trait,taxonomy$Genus)]

myplot=ggplot(correlations_total, aes(x = correlations, y = logp)) +
  geom_point(aes(color = qvalue < 0.1)) +
  scale_color_manual(values = c("darkgrey", "red"),
                     labels = c("Not Significant", "Significant")) +
  theme_bw() +
  labs(#title = paste0("Genus: ",gsub("g__","",genus_text1)," / Function type: ",gsub("_"," ",func)),
       #labs(title = paste0("Genus: ",gsub("g__","",genus)),
       x = "Pearson Correlation (rho): LFC (Ent1 vs Ent2) and Heritability",
       #y = "-Log10 P-value (Mann-Whitney U Test)",
       y = "-Log10 P-value",
       color = "") +
  theme(legend.position="bottom",
        panel.grid.minor = element_line(color = "grey90"), 
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=13),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5, size = 15))
plot(myplot)

# Add a column for text color based on Phylum
correlations_total$TextColor <- ifelse(correlations_total$Phylum == "p__Bacteroidota", "darkgreen",
                                       ifelse(correlations_total$Phylum == "p__Firmicutes_A", "blue", "black"))

# Add a column for bar fill based on qvalue
correlations_total$BarColor <- ifelse(correlations_total$qvalue < 0.10, "red", "gray")

# Plot with text and bar color customization
myplot <- ggplot(correlations_total, aes(x = reorder(trait, correlations, FUN = median), y = correlations)) +
  geom_bar(stat = "identity", aes(fill = BarColor)) +  # Assign fill color for bars
  scale_fill_identity() +  # Use the colors as defined in the BarColor column
  xlab("") +
  ylab("Pearson Correlation (rho): LFC (Ent1 vs Ent2) and Heritability") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 7, color = correlations_total$TextColor),  # Custom text colors
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_flip()
print(myplot)

###
mean_data_by_genus <- reshaped_data_species %>%
  group_by(Genus) %>%
  summarise(
    mean_enterotype_1 = mean(enterotype__1, na.rm = TRUE),
    mean_enterotype_2 = mean(enterotype__2, na.rm = TRUE)
  )

mean_data_by_genus$log_fold_change = log2(mean_data_by_genus$mean_enterotype_2 / mean_data_by_genus$mean_enterotype_1)
mean_data_by_genus$Phylum=taxonomy$Phylum[match(mean_data_by_genus$Genus,taxonomy$Genus)]

# Add a column for text color based on Phylum
mean_data_by_genus$TextColor <- ifelse(mean_data_by_genus$Phylum == "p__Bacteroidota", "darkgreen",
                                       ifelse(mean_data_by_genus$Phylum == "p__Firmicutes_A", "blue", "black"))

# Plot with text and bar color customization
myplot <- ggplot(mean_data_by_genus, aes(x = reorder(Genus, log_fold_change, FUN = median), y = log_fold_change)) +
  geom_bar(stat = "identity", aes(fill = "gray")) +  # Assign fill color for bars
  scale_fill_identity() +  # Use the colors as defined in the BarColor column
  xlab("") +
  ylab("LFC (Ent1 vs Ent2)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 7, color = mean_data_by_genus$TextColor),  # Custom text colors
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_flip()
print(myplot)


###

# Plot with text and bar color customization
myplot <- ggplot(reshaped_data_species, aes(x = reorder(trait, correlations, FUN = median), y = correlations)) +
  geom_bar(stat = "identity", aes(fill = BarColor)) +  # Assign fill color for bars
  scale_fill_identity() +  # Use the colors as defined in the BarColor column
  xlab("") +
  ylab("Pearson Correlation (rho): LFC (Ent1 vs Ent2) and Heritability") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 7, color = mean_data_by_genus$TextColor),  # Custom text colors
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_flip()
print(myplot)



prevotella=reshaped_data_species[grepl("s__Helicobacter",reshaped_data_species$Species),]
plot(prevotella$Heritability,prevotella$log_fold_change,xlab="Heritability (%)",ylab="LFC (Ent 1 to Ent 2)",main="CAG-485")
abline(lm(prevotella$log_fold_change ~ prevotella$Heritability), col = "red", lwd = 2)

cor_result=cor.test(prevotella$Heritability,prevotella$log_fold_change,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topright", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")




order=c("enterotype__1", "enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=host_reads)) +
  geom_boxplot() +
  xlab("") +
  ylab("Host Reads (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1", "enterotype__2")),
              map_signif_level = TRUE)
print(myplot)


myplot=ggplot(reshaped_data_species, aes(x = log_fold_change, y = Heritability)) +
  geom_point(aes(color = qvalue < 0.1)) +
  scale_color_manual(values = c("darkgrey", "red"),
                     labels = c("Not Significant", "Significant")) +
  theme_bw() +
  labs(#title = paste0("Genus: ",gsub("g__","",genus_text1)," / Function type: ",gsub("_"," ",func)),
    #labs(title = paste0("Genus: ",gsub("g__","",genus)),
    x = "LFC (Ent1 vs Ent2)",
    #y = "-Log10 P-value (Mann-Whitney U Test)",
    y = "Heritability (%)",
    color = "") +
  theme(legend.position="bottom",
        panel.grid.minor = element_line(color = "grey90"), 
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=13),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5, size = 15))
plot(myplot)

ent_taxa=data.frame(Species=reshaped_data_species$Species,
                    Phylum=reshaped_data_species$Phylum,
                    LFC=reshaped_data_species$log_fold_change,
                    Heritability=reshaped_data_species$Heritability,
                    All_Sig=reshaped_data_species$All_Sig)
ent_taxa$Enterotype[ent_taxa$LFC<0]="Enterotype 1"
ent_taxa$Enterotype[ent_taxa$LFC>0]="Enterotype 2"

order=c("Enterotype 1", "Enterotype 2")
myplot<-ggplot(ent_taxa, aes(x=factor(Enterotype, levels = order), y=100*Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Enterotype 1", "Enterotype 2")),
              map_signif_level = TRUE)
print(myplot)

stats=count(ent_taxa,Enterotype,Phylum,All_Sig)
View(stats)
