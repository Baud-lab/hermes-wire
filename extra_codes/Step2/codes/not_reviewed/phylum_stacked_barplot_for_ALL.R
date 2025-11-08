load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/filtered_prev_new.RData")

pdf("/users/abaud/fmorillo/paper_figures/stacked_barplots_ALL.pdf")

# Phylum
phylum_matrix=filtered_prev_objs[[3]]
# Define colors
total_phyla=c()
total_phyla=unique(c(total_phyla,rownames(phylum_matrix)))
#total_phyla=c(total_phyla,"Rare")

# Make into relative abundance
phylum_matrix <- apply(phylum_matrix, 2, function(i) i/sum(i))
colSums(phylum_matrix)[1:5]
# Define a cutoff for rare taxa
maxabundances <- apply(phylum_matrix, 1, max)
# Meanwhile, transpose the count table for future wrangling.
phylum_matrix <- data.frame(t(phylum_matrix))
# For every sample, sum up all rare taxa ( < 1% at their highest in this case)
threshold=0.01
rare=which(maxabundances < threshold)
if (length(rare)>1){
  phylum_matrix$Rare <- rowSums(phylum_matrix[,maxabundances < threshold], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
  phylum_matrix = phylum_matrix[,c(maxabundances > threshold, T) ] #`T` to include the Rare.
} else {
  colnames(phylum_matrix)[which(colnames(phylum_matrix)==names(maxabundances)[rare])]="Rare"
}

# Calculate mean values across covariates
phylum_df <- as.data.frame(phylum_matrix)
mean_values <- colMeans(phylum_df, na.rm = TRUE)  # na.rm = TRUE to handle NA values if any
mean_phylum_df <- data.frame(Phylum = names(mean_values), MeanValue = mean_values)
total_phyla=mean_phylum_df$Phylum
colours = metafolio::gg_color_hue(n =  length(total_phyla), c = 100, l = 65)
#colours <- qualitative_hcl(n = length(total_phyla), h = c(0, 360), c = 100, l = 65)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]

# Plot
myplot <- mean_phylum_df %>%
  ggplot(aes(x="", y = MeanValue, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, width = 1) +
  scale_fill_manual(values = colours,
                    labels = total_phyla,
                    limits = total_phyla) +
  labs(color = "") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  xlab("") +
  ylab("Cumulative mean relative abundance") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10),)
print(myplot)

##########

# Family
family_matrix=filtered_prev_objs[[12]]
# Define colors
total_family=c()
total_family=unique(c(total_family,rownames(family_matrix)))
#total_family=c(total_family,"Rare")

# Make into relative abundance
family_matrix <- apply(family_matrix, 2, function(i) i/sum(i))
colSums(family_matrix)[1:5]
# Define a cutoff for rare taxa
maxabundances <- apply(family_matrix, 1, max)
# Meanwhile, transpose the count table for future wrangling.
family_matrix <- data.frame(t(family_matrix))
# For every sample, sum up all rare taxa ( < 1% at their highest in this case)
threshold=0.2
rare=which(maxabundances < threshold)
if (length(rare)>1){
  family_matrix$Rare <- rowSums(family_matrix[,maxabundances < threshold], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
  family_matrix = family_matrix[,c(maxabundances > threshold, T) ] #`T` to include the Rare.
} else {
  colnames(family_matrix)[which(colnames(family_matrix)==names(maxabundances)[rare])]="Rare"
}

# Calculate mean values across covariates
family_df <- as.data.frame(family_matrix)
mean_values <- colMeans(family_df, na.rm = TRUE)  # na.rm = TRUE to handle NA values if any
mean_family_df <- data.frame(Family = names(mean_values), MeanValue = mean_values)
total_family=mean_family_df$Family
colours = metafolio::gg_color_hue(n =  length(total_family), c = 100, l = 65)
#colours <- qualitative_hcl(n = length(total_family), h = c(0, 360), c = 100, l = 65)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]

# Plot
myplot <- mean_family_df %>%
  ggplot(aes(x="", y = MeanValue, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, width = 1) +
  scale_fill_manual(values = colours,
                    labels = total_family,
                    limits = total_family) +
  labs(color = "") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  xlab("") +
  ylab("Cumulative mean relative abundance") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10),)
print(myplot)


##########

# Genus
genus_matrix=filtered_prev_objs[[15]]
# Define colors
total_genus=c()
total_genus=unique(c(total_genus,rownames(genus_matrix)))
#total_genus=c(total_genus,"Rare")

# Make into relative abundance
genus_matrix <- apply(genus_matrix, 2, function(i) i/sum(i))
colSums(genus_matrix)[1:5]
# Define a cutoff for rare taxa
maxabundances <- apply(genus_matrix, 1, max)
# Meanwhile, transpose the count table for future wrangling.
genus_matrix <- data.frame(t(genus_matrix))
# For every sample, sum up all rare taxa ( < 1% at their highest in this case)
threshold=0.1
rare=which(maxabundances < threshold)
if (length(rare)>1){
  genus_matrix$Rare <- rowSums(genus_matrix[,maxabundances < threshold], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
  genus_matrix = genus_matrix[,c(maxabundances > threshold, T) ] #`T` to include the Rare.
} else {
  colnames(genus_matrix)[which(colnames(genus_matrix)==names(maxabundances)[rare])]="Rare"
}

# Calculate mean values across covariates
genus_df <- as.data.frame(genus_matrix)
mean_values <- colMeans(genus_df, na.rm = TRUE)  # na.rm = TRUE to handle NA values if any
mean_genus_df <- data.frame(Genus = names(mean_values), MeanValue = mean_values)
total_genus=mean_genus_df$Genus
colours = metafolio::gg_color_hue(n =  length(total_genus), c = 100, l = 65)
#colours <- qualitative_hcl(n = length(total_genus), h = c(0, 360), c = 100, l = 65)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]

# Plot
myplot <- mean_genus_df %>%
  ggplot(aes(x="", y = MeanValue, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, width = 1) +
  scale_fill_manual(values = colours,
                    labels = total_genus,
                    limits = total_genus) +
  labs(color = "") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  xlab("") +
  ylab("Cumulative mean relative abundance") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10),)
print(myplot)

dev.off()
