library(ggvenn)

ranks=c("Phylum","Class","Order","Family","Genus","Species")
subset_ids=c("ALL")
module_comparisons="YES"
methods=c("16S","Shallow")
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
  assign(paste('filtered_prev',concatenated[i],sep='_'), filtered_prev)
}

pdf("/users/abaud/fmorillo/paper_figures/venn_16S_shallow.pdf")
for (rank in ranks){#[-length(ranks)]){
  print(rank)
  matrix_SS=get(paste('filtered_prev',rank,"ALL_Shallow",sep="_"))
  matrix_16S=get(paste('filtered_prev',rank,"ALL_16S",sep="_"))
  taxa_SS=rownames(matrix_SS)
  taxa_16S=rownames(matrix_16S)
  taxa_list <- list(`Shallow Shotgun` = taxa_SS, `16S` = taxa_16S)
  # Generate and save the Venn diagram
  myplot=ggvenn(taxa_list, 
         fill_color = c("lightcoral", "lightblue"), 
         stroke_size = 0.5, 
         set_name_size = 0,
         text_size = 5,
         text_color = "black",
         show_percentage = TRUE) +
    labs(title = NULL) +
    annotate("text", x = 0, y = 1.6, label = paste0("Venn Diagram: ",rank), size = 6, fontface = "bold", hjust = 0.5) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(t = 1, r = 1, b = 3, l = 1, "cm")  # Adjust margins to provide space
    ) +
    annotate("text", x = -0.7, y = 1.2, label = "Shallow Shotgun", size = 6, hjust = 0.5) + 
    annotate("text", x = 0.7, y = 1.2, label = "16S", size = 6, hjust = 0.5)
  print(myplot)
}
dev.off()
