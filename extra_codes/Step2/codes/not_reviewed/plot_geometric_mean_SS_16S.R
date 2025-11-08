library(ggplot2)
library(tidyverse)
library(compositions)
library(coda.base)

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

#######

pdf("/users/abaud/fmorillo/paper_figures/geomean_SS_16S.pdf")
for (rank in ranks){
  # Add an offset to deal with zeroes
  for (method in methods) {
    matrix=get(paste('filtered_prev',rank,subset_id,method,harmonization[1],sep="_"))
    offset = 0.00001
    matrix = matrix + offset
    #matrix=as.matrix(cmultRepl(matrix, method = 'SQ'))
    centre <- coda.base::center(t(matrix))
    assign(paste('geomean',method,sep='_'), centre)
  }
  
  # Convert to data frames for easy manipulation
  geomean_16S_df <- data.frame(Trait = names(geomean_16S), geomean_16S = geomean_16S)
  #geomean_16S_df=geomean_16S_df[geomean_16S_df$Trait!="g__",]
  geomean_Shallow_df <- data.frame(Trait = names(geomean_Shallow), geomean_Shallow = geomean_Shallow)

  # Join data frames
  GM_mean <- full_join(geomean_16S_df, geomean_Shallow_df, by = 'Trait')
  log2_GM_mean <- log2(GM_mean[,2:3])
  rownames(log2_GM_mean)=GM_mean$Trait
  
  # Handle NA values
  GM_mean[is.na(GM_mean)] <- 0
  log2_GM_mean[is.na(log2_GM_mean)] <- 0
  
  # Perform inner join and calculate Spearman correlation
  GM_mean_both <- inner_join(geomean_16S_df, geomean_Shallow_df, by = 'Trait')
  log2_GM_both <- log2(GM_mean_both[,2:3])
  #cor(GM_mean_both$geomean_16S, GM_mean_both$geomean_Shallow, method = 'spearman')
  ro=round(cor(log2_GM_both$geomean_16S, log2_GM_both$geomean_Shallow, method = 'spearman'),digits = 3)
  #ro=round(cor( 
  #  log10_GM_mean2_mean_counts_16S[log10_GM_mean2_mean_counts_Shallow != 0 & log10_GM_mean2_mean_counts_16S != 0],
  #  log10_GM_mean2_mean_counts_Shallow[log10_GM_mean2_mean_counts_Shallow != 0 & log10_GM_mean2_mean_counts_16S != 0],
  #  method = 'spearman'
  #),digits = 3)
  # 0.6167
  
  ## Plotting
  # Prepare log-transformed data for plotting
  log10_GM_mean2_mean_counts_16S <- -log10(2^(log2_GM_mean$geomean_16S))
  log10_GM_mean2_mean_counts_Shallow <- -log10(2^(log2_GM_mean$geomean_Shallow))
  
  ### Adjust values for better visualization in plots
  #log10_GM_mean2_mean_counts_16S[log10_GM_mean2_mean_counts_16S<(-4.5)] <- -4.61
  #log10_GM_mean2_mean_counts_Shallow[log10_GM_mean2_mean_counts_Shallow<(-8.5)] <- -8.5
  
  # Set up plot layout
  #par(mfrow = c(3,1))
  
  # Plot
  col_plot <- rep(rgb(0.4,0.4,0.8,0.6), length(log2_GM_mean$geomean_16S))
  col_plot[log2_GM_mean$geomean_16S == 0] <- rgb(0, 0.4, 0, 0.6)
  col_plot[log2_GM_mean$geomean_Shallow == 0] <- rgb(0.8, 0, 0, 0.6)
  plot(log10_GM_mean2_mean_counts_16S, log10_GM_mean2_mean_counts_Shallow, axes = FALSE, xlab = '-log10 (16S)', ylab = '-log10 (Shallow Shotgun)', col = col_plot,
       pch = 16, 
       #cex = 1.5,
       main = paste0(rank,' abundance of 16S and Shotgun'),
       #xlim = c(-4.5,0),
       #ylim = c(-8.5,0),
       #cex.main = 2.2,
       #cex.lab = 1.9
  )
  abline(lm(log10_GM_mean2_mean_counts_Shallow[log10_GM_mean2_mean_counts_Shallow != 0 & log10_GM_mean2_mean_counts_16S != 0] ~ 
              log10_GM_mean2_mean_counts_16S[log10_GM_mean2_mean_counts_Shallow != 0 & log10_GM_mean2_mean_counts_16S != 0]), 
         col = 'slateblue3', lwd = 2)
  box()
  axis(1, #at = seq(-5, 0, 1), #labels = c(expression(10**-5),  expression(10**-4), expression(10**-3), expression(10**-2), expression(10**-1), expression(10**0)), 
       cex.axis = 1.2)
  axis(2, #at = seq(-8, 0, 1), labels = c(expression(10**-8),'', expression(10**-6),'', expression(10**-4),'', expression(10**-2),'', expression(10**-0)), las = 1, 
       cex.axis = 1.2)
  legend("topright", legend = c(paste0('Common ',rank), 'Not detected by Shotgun', 'Not detected by 16S'), 
         col = c(rgb(0.4,0.4,0.8,0.7), rgb(0.8, 0, 0, 0.7), rgb(0, 0.4, 0, 0.7)),
         cex = 0.9, pch = 16)
  # Get current plot limits using par("usr")
  plot_limits <- par("usr")
  x_pos <- plot_limits[2] * 0.85  # 95% of the x-axis range (right side)
  y_pos <- plot_limits[4] * 0.30  # 95% of the y-axis range (top)
  
  # Display the rho value in the top-right of the plot
  text(x_pos, y_pos, bquote(rho == .(ro)), col = 'slateblue3', cex = 1.2)
  
  print(rank)
}
dev.off()
