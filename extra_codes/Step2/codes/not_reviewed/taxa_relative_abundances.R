metadata <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")
sample_identifier="host_subject_id"

covariates=c("Sex","Study")
subset_id="ALL"
subset_ids=c("MI","NY","ALL")
ranks=c("Phylum","Class","Order","Family","Genus","Species")

# Calculate relative abundances and prevalence
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Cluster_Analyses/filtered_prev.RData")

permutations <- expand.grid(subset_ids,ranks)
colnames(permutations)=c("subset_ids","ranks")
concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
for (i in 1:length(concatenated)){
  filtered_prev=filtered_prev_objs[[i]]
  assign(paste('filtered_prev',concatenated[i],sep='_'), filtered_prev)
}

ranks=c("Phylum","Family","Genus")
pdf("/users/abaud/fmorillo/paper_figures/taxa_relative_abundances.pdf")
for (rank in ranks){
  # Make into relative abundance
  matrix <- apply(get(paste("filtered_prev",rank,subset_id,sep="_")), 2, function(i) i/sum(i))
  colSums(matrix)[1:5]
  # Define a cutoff for rare taxa
  maxabundances <- apply(matrix, 1, max)
  # Define mean relative abundances
  row_means=rowMeans(matrix)
  row_means <- data.frame(taxa=rownames(matrix), MeanValue = row_means)
  assign(paste("row_means",rank,sep="_"),row_means)
  # Meanwhile, transpose the count table for future wrangling.
  matrix <- data.frame(t(matrix))
  # For every sample, sum up all rare taxa ( < 1% at their highest in this case)
  if (rank=="Phylum"){
    threshold=0
  } else {
    threshold=0.08
  }
  rare=which(maxabundances < threshold)
  if (length(rare)>1){
    matrix$`Rare Taxa` <- rowSums(matrix[,maxabundances < threshold], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
    matrix = matrix[,c(maxabundances > threshold, T) ] #`T` to include the `Rare Taxa`.
  } else {
    colnames(matrix)[which(colnames(matrix)==names(maxabundances)[rare])]=paste0("Rare ",rank)
  }
  # Prepare the data for ggplot by adding in metadata here
  motch = match(rownames(matrix), metadata[[sample_identifier]])
  any(is.na(motch))
  this_metadata = metadata[motch,]
  covariates=c("Sex","Study")
  formula <- c("Sample")
  matrix[[covariates[1]]]=this_metadata[[covariates[1]]]
  matrix[[covariates[2]]]=this_metadata[[covariates[2]]]
  matrix$Sample=rownames(matrix)
  # Wrangle the data to long format for easy plotting
  formula <- c(formula, "Study","Sex")
  barlong <- matrix %>%
    pivot_longer(!formula, names_to = rank, values_to = "value") %>%
    mutate(!!rank := str_replace(!!sym(rank), ".*_or_", ""))
  # Calculate mean values across covariates
  grouped_matrix <- barlong %>%
    group_by(!!sym(rank), Study, Sex) %>%
    summarise(mean_value_per_rank = mean(value))
  grouped_matrix$Sex[grouped_matrix$Sex=="F"]="Female"
  grouped_matrix$Sex[grouped_matrix$Sex=="M"]="Male"
  grouped_matrix[[rank]]=gsub("\\.","-",grouped_matrix[[rank]])
  ###
  # Calculate by covariate
  ## Sex
  grouped_matrix_Sex <- barlong %>%
    group_by(!!sym(rank), Sex) %>%
    summarise(mean_value_per_rank = mean(value))
  grouped_matrix_Sex$Sex[grouped_matrix_Sex$Sex=="F"]="Female"
  grouped_matrix_Sex$Sex[grouped_matrix_Sex$Sex=="M"]="Male"
  grouped_matrix_Sex[[rank]]=gsub("\\.","-",grouped_matrix_Sex[[rank]])
  assign(paste("grouped_matrix_Sex",rank,sep="_"),grouped_matrix_Sex)
  ## Study
  grouped_matrix_Study <- barlong %>%
    group_by(!!sym(rank), Study) %>%
    summarise(mean_value_per_rank = mean(value))
  grouped_matrix_Study[[rank]]=gsub("\\.","-",grouped_matrix_Study[[rank]])
  assign(paste("grouped_matrix_Study",rank,sep="_"),grouped_matrix_Study)
  ###
  levels_rank <- levels(factor(grouped_matrix[[rank]]))
  set.seed(120)
  num_colors <- length(levels_rank)
  all_colors <- colors()
  random_colors <- sample(all_colors, num_colors, replace = FALSE)
  names(random_colors) <- levels_rank
  myplot <- grouped_matrix %>%
    ggplot(aes(x = Sex, y = mean_value_per_rank, fill = !!sym(rank))) +
    geom_bar(stat = "identity", position = "stack", col = "white", linewidth = 2, width = 1) +
    facet_wrap(. ~ Study) +
    labs(color = "") +
    ggtitle(paste0("Subset: ",subset_id)) +
    xlab(paste0("")) +
    ylab("Mean relative abundance") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = random_colors) +
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
  myplot <- grouped_matrix_Study %>%
    ggplot(aes(x = Study, y = mean_value_per_rank, fill = !!sym(rank))) +
    geom_bar(stat = "identity", position = "stack", col = "white", linewidth = 2, width = 1) +
    labs(color = "") +
    ggtitle(paste0("Subset: ",subset_id)) +
    xlab(paste0("")) +
    ylab("Mean relative abundance") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = random_colors) +
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
}
dev.off()

#####

filtered_df <- grouped_matrix_Study_Genus %>%
  filter(Study %in% c("MI", "NY"))
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = Study, values_from = mean_value_per_rank, names_prefix = "mean_value_per_rank_")
diff_dataset_Genus <- pivoted_df %>%
  mutate(difference_MI_NY = mean_value_per_rank_MI - mean_value_per_rank_NY) %>%
  select(Genus, difference_MI_NY)

filtered_df <- grouped_matrix_Sex_Phylum %>%
  filter(Sex %in% c("Male", "Female"))
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = Sex, values_from = mean_value_per_rank, names_prefix = "mean_value_per_rank_")
diff_dataset_Phylum <- pivoted_df %>%
  mutate(difference_sex = mean_value_per_rank_Male - mean_value_per_rank_Female) %>%
  select(Phylum, difference_sex)
