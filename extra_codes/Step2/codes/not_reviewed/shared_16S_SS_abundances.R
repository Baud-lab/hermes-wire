# 1. Read gzipped or flat files
read_gzipped <- function(input_file) {
  input <- input_file
  if (file_ext(input_file)=="gz") {
    input = gzfile(input_file)  
  }
  return(input)
}

read_tab_file <- function(input_file, header=FALSE) {
  file_content <- read.delim(read_gzipped(input_file), header=header, check.names=FALSE)		
  return(file_content)
}

read_csv_file <- function(input_file) {
  file_content <- read.csv(read_gzipped(input_file), stringsAsFactors = FALSE)	
  return(file_content)
}

# 4. Data transformation (CLR)
calculate_gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}
do_clr_default = function(data){
  transform_col = function(x){
    gMean = calculate_gm_mean(x)
    transformed_col = log(x / gMean)
    return(transformed_col)
  }
  transformed_data = apply(data, MAR = 2, FUN = transform_col)
  return(transformed_data)
}

ranks=c("Phylum","Class","Order","Family","Genus","Species")
subset_ids=c("ALL")
subset_id="ALL"
module_comparisons="YES"
methods=c("16S","Shallow")
harmonization=c("sample_harmonization","sample_and_taxa_harmonization")
sample_identifier= "host_subject_id"

for (harm in harmonization){
  # Load counts
  if (harm=="sample_harmonization") {
    load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Sample_Harm/Harmonization/harmonized_samples.RData")
    if (module_comparisons=="YES"){
      permutations <- expand.grid(subset_ids,ranks,methods)
      colnames(permutations)=c("subset_ids","ranks","methods")
      concatenated <- paste(permutations$ranks, permutations$subset_ids, permutations$methods,sep = "_")
    } else {
      permutations <- expand.grid(subset_ids,ranks)
      colnames(permutations)=c("subset_ids","ranks")
      concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
    }
    for (i in 1:length(concatenated)){
      filtered_prev=filtered_prev_objs[[i]]
      assign(paste('filtered_prev',concatenated[i],harm,sep='_'), filtered_prev)
    }
  } else {
    load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Sample_Harm/Harmonization/harmonized_samples_and_taxa.RData")
    if (module_comparisons=="YES"){
      permutations <- expand.grid(subset_ids,ranks,methods)
      colnames(permutations)=c("subset_ids","ranks","methods")
      concatenated <- paste(permutations$ranks, permutations$subset_ids, permutations$methods,sep = "_")
    } else {
      permutations <- expand.grid(subset_ids,ranks)
      colnames(permutations)=c("subset_ids","ranks")
      concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
    }
    for (i in 1:length(concatenated)){
      filtered_prev=filtered_prev_objs[[i]]
      assign(paste('filtered_prev',concatenated[i],harm,sep='_'), filtered_prev)
    }
    
    ## Residuals
    #load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Matrix_Processing/residuals_qned_counts.RData")
  }
}

pdf("/users/abaud/fmorillo/paper_figures/16S_SS_scatter_plots.pdf")
for (rank in ranks){
  #ab_vs_prev_total=data.frame()
  rank="Genus" #### Remove
  for (method in methods) {
    matrix=get(paste('filtered_prev',rank,subset_id,method,harmonization[2],sep="_"))
    # Matrix relative abundance %
    matrix.relative <- sweep(matrix,2,colSums(matrix),"/")
    matrix.relative=matrix.relative*100
    assign(paste("matrix.relative",rank,method,sep="_"),matrix.relative)
    # Mean relative abundances and prevalence
    all_harm=get(paste('filtered_prev',rank,subset_id,method,harmonization[2],sep="_"))
    ab_vs_prev=data.frame(trait=rownames(all_harm),
                          abundance=100*apply(sweep(all_harm,2,colSums(all_harm),"/"), 1, mean, na.rm=TRUE),
                          prevalence=100*rowSums(all_harm>0)/ncol(all_harm))
    assign(paste("ab_vs_prev",rank,method,sep="_"),ab_vs_prev)
    # CLR
    offset = 0.00001
    matrix = matrix + offset
    clr_counts = do_clr_default(matrix)
    assign(paste('filtered_clr_counts',method,sep='_'), clr_counts)
  }
  ## Scatter plot - Relative abundances
  plot(matrix.relative_16S,matrix.relative_Shallow, xlab="Relative abundance - 16S (%)", ylab="Relative abundance - Shallow Shotgun (%)", cex.lab = 1.3, main=paste0("Rank: ",rank))
  cor_result=cor.test(matrix.relative_16S,matrix.relative_Shallow,method="spearman")
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
  legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
  
  ## Scatter plot - CLR
  plot(filtered_clr_counts_16S,filtered_clr_counts_Shallow, xlab="Relative abundance - 16S (CLR)", ylab="Relative abundance - Shallow Shotgun (CLR)", cex.lab = 1.3, main=paste0("Rank: ",rank))
  cor_result=cor.test(filtered_clr_counts_16S,filtered_clr_counts_Shallow,method="spearman")
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
  legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
  
  ## Scatter plot - Mean relative abundances
  plot(ab_vs_prev_16S$abundance,ab_vs_prev_Shallow$abundance, xlab="Mean relative abundance - 16S (%)", ylab="Mean relative abundance - Shallow Shotgun (%)", cex.lab = 1.3, main=paste0("Rank: ",rank))
  cor_result=cor.test(ab_vs_prev_16S$abundance,ab_vs_prev_Shallow$abundance,method="spearman")
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
  legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
  
  ## Scatter plot - Prevalence
  plot(ab_vs_prev_16S$prevalence,ab_vs_prev_Shallow$prevalence, xlab="Prevalence - 16S (%)", ylab="Prevalence - Shallow Shotgun (%)", cex.lab = 1.3, main=paste0("Rank: ",rank))
  cor_result=cor.test(ab_vs_prev_16S$prevalence,ab_vs_prev_Shallow$prevalence,method="spearman")
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
  legend("bottomright", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
}
dev.off()

######
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Heritability/heritability.RData")
all_shallow=all_VCs_full
all_shallow=all_shallow[grepl("ALL",all_shallow$subset_id),]
all_shallow$method="Shallow"
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Heritability/heritability.RData")
all_16S=all_VCs_full
all_16S=all_16S[grepl("ALL",all_16S$subset_id),]
all_16S$method="16S"

all_16S=all_16S[intersect(rownames(all_16S),rownames(all_shallow)),]
all_shallow=all_shallow[intersect(rownames(all_16S),rownames(all_shallow)),]


for (rank in ranks){
  ab_vs_prev_16S <- get(paste("ab_vs_prev",rank,"16S",sep="_"))
  ab_vs_prev_Shallow <- get(paste("ab_vs_prev",rank,"Shallow",sep="_"))
  all_16S=all_16S[grepl(rank,all_16S$rank),]
  all_shallow=all_shallow[grepl(rank,all_shallow$rank),]
  matrix_16S=get(paste("matrix.relative",rank,"16S",sep="_"))
  matrix_Shallow=get(paste("matrix.relative",rank,"Shallow",sep="_"))
  correlations_abundance_total=data.frame()
  pdf("/users/abaud/fmorillo/paper_figures/genus_abundance_16S_SS.pdf")
  for (trait in rownames(matrix_16S)){
    print(trait)  
    cor_result=cor.test(matrix_16S[rownames(matrix_16S)==trait,],matrix_Shallow[rownames(matrix_Shallow)==trait,],method="pearson")
    r_value <- cor_result$estimate
    p_value <- cor_result$p.value
    correlations_abundance=data.frame(trait=trait,r_value_abundance=r_value,p_value_abundance=p_value)
    correlations_abundance_total=rbind(correlations_abundance_total,correlations_abundance)
    
    plot(matrix_16S[rownames(matrix_16S)==trait,],
         matrix_Shallow[rownames(matrix_Shallow)==trait,],
         main=paste("Correlations of Relative Abundances (%) - ",rank,": ",gsub("g__","",trait)),
         xlab="16S",
         ylab="Shallow")
    abline(lm(matrix_Shallow[rownames(matrix_Shallow)==trait,] ~ matrix_16S[rownames(matrix_16S)==trait,]), col="red", lwd=2)
    if (p_value < 0.01){
      P="P < 0.01"
    } else {
      if (p_value > 0.01 && p_value < 0.05){
        P="P < 0.05"
      } else {
        P="P > 0.05"
      }
    }
    legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
    
  }
  dev.off()
  
  correlations_abundance_total$heritability_16S=round(all_16S$var_Ad[match(correlations_abundance_total$trait,all_16S$trait)]*100,digit=2)
  correlations_abundance_total$heritability_Shallow=round(all_shallow$var_Ad[match(correlations_abundance_total$trait,all_shallow$trait)]*100,digit=2)
  correlations_abundance_total$diff_heritability=abs(correlations_abundance_total$heritability_16S-correlations_abundance_total$heritability_Shallow)
  correlations_abundance_total$prevalence_16S=ab_vs_prev_16S$prevalence[match(correlations_abundance_total$trait,ab_vs_prev_16S$trait)]
  correlations_abundance_total$prevalence_Shallow=ab_vs_prev_Shallow$prevalence[match(correlations_abundance_total$trait,ab_vs_prev_Shallow$trait)]
  correlations_abundance_total$diff_prevalence=abs(correlations_abundance_total$prevalence_16S-correlations_abundance_total$prevalence_Shallow)
  correlations_abundance_total$avg_abundance_16S=ab_vs_prev_16S$abundance[match(correlations_abundance_total$trait,ab_vs_prev_16S$trait)]
  correlations_abundance_total$avg_abundance_Shallow=ab_vs_prev_Shallow$abundance[match(correlations_abundance_total$trait,ab_vs_prev_Shallow$trait)]
  correlations_abundance_total$diff_abundance=abs(correlations_abundance_total$avg_abundance_16S-correlations_abundance_total$avg_abundance_Shallow)
  assign(paste("correlations_abundance",rank,sep="_"),correlations_abundance_total)
  
  #cor_herit_abundance_cor=cor.test(correlations_abundance_total$r_value,correlations_abundance_total$diff_heritability,method="spearman")
  plot(correlations_abundance_total$r_value,correlations_abundance_total$diff_heritability)
  #cor_herit_abundance_dif=cor.test(correlations_abundance_total$diff_heritability, correlations_abundance_total$diff_abundance,method="spearman")
  #cor_herit_prevalence=cor.test(correlations_abundance_total$diff_prevalence, correlations_abundance_total$diff_heritability,method="spearman")
  #cor_abundance_prevalence=cor.test(correlations_abundance_total$diff_prevalence, correlations_abundance_total$diff_abundance,method="spearman")
}

save(correlations_abundance_Genus,file="/users/abaud/fmorillo/paper_figures/correlations_abundance_Genus.RData")

## Filter common "trait" with "prevalence" > 95% in both dataframes
#filtered_16S <- subset(get(paste("ab_vs_prev",rank,"16S",sep="_")), prevalence > 90)
#filtered_Shallow <- subset(ab_vs_prev_Genus_Shallow, prevalence > 90)
#
## Find common "trait" between both filtered dataframes
#common_traits <- intersect(filtered_Genus_16S$trait, filtered_Genus_Shallow$trait)
#common_traits <-ab_vs_prev_Genus_16S$trait



######




pdf("/users/abaud/fmorillo/paper_figures/16S_SS_abundance_prevalence.pdf")
for (rank in ranks[-length(ranks)]){
  all_harm=get(paste('filtered_prev',rank,subset_id,methods[1],harmonization[2],sep="_"))
  shared=rownames(all_harm)
  ab_vs_prev_shared=data.frame()
  for (method in methods){
    sample_harm=get(paste('filtered_prev',rank,subset_id,method,harmonization[1],sep="_"))
    ab_vs_prev=data.frame(trait=rownames(sample_harm),
                                  abundance=100*apply(sweep(sample_harm,2,colSums(sample_harm),"/"), 1, mean, na.rm=TRUE),
                                  prevalence=100*rowSums(sample_harm>0)/ncol(sample_harm))
    ab_vs_prev <- ab_vs_prev %>%
      mutate(shared = ifelse(trait %in% shared, "Shared", "Not shared"))
    if (rank == "Genus" || rank == "Species"){
      max_ab=20
    } else {
      if (rank == "Family" || rank == "Order"){
        max_ab=50
      } else {
        max_ab=100
      }
    }
    plot(x=1,
         type = "n",
         ylab=paste0("Mean Abundance (%) (Method: ",method," / Rank: ",rank,")"),
         xlab=paste0("Prevalence (%) (Method: ",method," / Rank: ",rank,")"),
         ylim=c(0,max_ab),
         xlim=c(50,100))
    points(x = ab_vs_prev$prevalence[ab_vs_prev$shared == "Not shared"],
           y = ab_vs_prev$abundance[ab_vs_prev$shared == "Not shared"],
           pch = 19,
           col = "grey")
    points(x = ab_vs_prev$prevalence[ab_vs_prev$shared == "Shared"],
           y = ab_vs_prev$abundance[ab_vs_prev$shared == "Shared"],
           pch = 19,
           col = "blue")
    edge_labels=c("Shared","Not shared")
    significance_colors=c("blue","grey")
    legend("top", legend = edge_labels, col = significance_colors, pch = 15, horiz = TRUE, xpd = TRUE, inset = c(0, -0.1))
    order=c("Shared","Not shared")
    myplot<-ggplot(ab_vs_prev, aes(x=factor(shared, levels = order), y=abundance)) +
      geom_boxplot(
        #outlier.shape = NA
                   ) +
      xlab("") +
      ylab("Relative Abundance (%)") +
      ggtitle(paste0("Method: ",method," / Rank:",rank)) +
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
      geom_signif(comparisons = list(c("Shared","Not shared")),
                  map_signif_level = TRUE)
    print(myplot)
    myplot<-ggplot(ab_vs_prev, aes(x=factor(shared, levels = order), y=prevalence)) +
      geom_boxplot(
        #outlier.shape = NA
      ) +
      xlab("") +
      ylab("Prevalence (%)") +
      ggtitle(paste0("Method: ",method," / Rank:",rank)) +
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
      geom_signif(comparisons = list(c("Shared","Not shared")),
                  map_signif_level = TRUE)
    print(myplot)
    # Compare shared
    ab_vs_prev$method=method
    ab_vs_prev=ab_vs_prev[ab_vs_prev$trait %in% shared,]
    ab_vs_prev_shared=rbind(ab_vs_prev_shared,ab_vs_prev)
    print(paste0(method,": ",rank))
  }
  order=c("16S","Shallow")
  myplot<-ggplot(ab_vs_prev_shared, aes(x=factor(method, levels = order), y=abundance)) +
    geom_boxplot(
      #outlier.shape = NA
    ) +
    xlab("") +
    ylab("Relative Abundance (%)") +
    ggtitle(paste0("Shared taxa: ",rank)) +
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
    geom_signif(comparisons = list(c("16S","Shallow")),
                map_signif_level = TRUE)
  print(myplot)
  
myplot<-ggplot(ab_vs_prev_shared, aes(x=factor(method, levels = order), y=prevalence)) +
    geom_boxplot(
      #outlier.shape = NA
    ) +
    xlab("") +
    ylab("Prevalence (%)") +
    ggtitle(paste0("Shared taxa: ",rank)) +
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
    geom_signif(comparisons = list(c("16S","Shallow")),
                map_signif_level = TRUE)
  print(myplot)
}
dev.off()


pdf("/users/abaud/fmorillo/paper_figures/16S_SS_proportion_of_zeroes.pdf")
for (harm in harmonization){
  if (harm=="sample_harmonization") {
    title="Only samples harmonized"
  } else {
    title="Samples and taxa harmonized"
  }
  all_ranks=data.frame()
  for (rank in ranks){
    proportion_zeroes_all=data.frame()
    for (method in methods){
      sample_harm=get(paste('filtered_prev',rank,subset_id,method,harm,sep="_"))
      proportion_zeroes <- apply(sample_harm, 2, function(column) {
        sum(column == 0) / length(column)
      })
      proportion_zeroes=data.frame(sample=colnames(sample_harm),proportion=100*proportion_zeroes)
      proportion_zeroes$method=method
      proportion_zeroes$rank=rank
      proportion_zeroes_all=rbind(proportion_zeroes_all,proportion_zeroes)
      all_ranks=rbind(all_ranks,proportion_zeroes)
    }
    order=c("16S","Shallow")
    p_value <- wilcox.test(proportion_zeroes_all$proportion[proportion_zeroes_all$method == "16S"], 
                           proportion_zeroes_all$proportion[proportion_zeroes_all$method == "Shallow"])$p.value
    myplot<-ggplot(proportion_zeroes_all, aes(x=factor(method, levels = order), y=proportion)) +
      geom_boxplot(
        #outlier.shape = NA
      ) +
      xlab("") +
      ylab("Proportion of zeroes (%)") +
      ggtitle(paste0(title,": ",rank)) +
      coord_cartesian(ylim = c(0, 100)) +
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
      geom_signif(comparisons = list(c("16S","Shallow")), 
                  annotations = sprintf("p = %.2e", p_value),
                  test = "wilcox.test",
                  y_position = 95, map_signif_level = FALSE)
    print(myplot)
    myplot = ggplot(proportion_zeroes_all, aes(x = proportion, color = method, fill = method)) +
      geom_density(alpha = 0.5) +
      labs(fill = "Sequencing Method", color = "Sequencing Method") +
      labs(title = paste0(title,": ",rank),
           x = "Proportion of zeroes (%)",
           y = "Density") +
      scale_fill_manual(values = c("blue", "red")) +
      scale_color_manual(values = c("blue", "red")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            legend.position = "bottom")
    print(myplot)
    print(paste0(title,": ",rank))
  }
  order=c("Phylum","Class","Order","Family","Genus","Species")
  myplot<-ggplot(all_ranks, aes(x=factor(rank, levels = order), y=proportion)) +
    geom_boxplot(
      #outlier.shape = NA
    ) +
    facet_grid(. ~ method)+
    xlab("") +
    ylab("Proportion of zeroes (%)") +
    ggtitle(paste0(title)) +
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
    scale_x_discrete(guide = guide_axis(angle = 60))
  print(myplot)
}
dev.off()


########

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Matrix_Processing/filtered_prev.RData")
filtered_prev_Species=filtered_prev_objs[[6]]
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Matrix_Processing/filtered_prev.RData")
filtered_prev_EC4=filtered_prev_objs[[4]]

methods=c("Species","EC4")

pdf("/users/abaud/fmorillo/paper_figures/tax_func_proportion_of_zeroes.pdf")

proportion_zeroes_all=data.frame()
for (method in methods){
  sample_harm=get(paste('filtered_prev',method,sep="_"))
  proportion_zeroes <- apply(sample_harm, 2, function(column) {
    sum(column == 0) / length(column)
  })
  proportion_zeroes=data.frame(sample=colnames(sample_harm),proportion=100*proportion_zeroes)
  proportion_zeroes$method=method
  proportion_zeroes_all=rbind(proportion_zeroes_all,proportion_zeroes)
}
order=methods
p_value <- wilcox.test(proportion_zeroes_all$proportion[proportion_zeroes_all$method == "Species"], 
                       proportion_zeroes_all$proportion[proportion_zeroes_all$method == "EC4"])$p.value
myplot<-ggplot(proportion_zeroes_all, aes(x=factor(method, levels = order), y=proportion)) +
  geom_boxplot(
    #outlier.shape = NA
  ) +
  xlab("") +
  ylab("Proportion of zeroes (%)") +
  #ggtitle(paste0(title,": ",rank)) +
  coord_cartesian(ylim = c(0, 50)) +
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
  geom_signif(comparisons = list(c("Species","EC4")), 
              annotations = sprintf("p = %.2e", p_value),
              test = "wilcox.test",
              y_position = 48, map_signif_level = FALSE)
print(myplot)


myplot = ggplot(proportion_zeroes_all, aes(x = proportion, color = method, fill = method)) +
  geom_density(alpha = 0.5) +
  labs(fill = "Sequencing Method", color = "Sequencing Method") +
  labs(x = "Proportion of zeroes (%)",
       y = "Density") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "bottom")
print(myplot)

dev.off()
