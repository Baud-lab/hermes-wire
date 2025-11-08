subset_id="ALL"
# "ALL" "MI" "NY"
method="Shallow"

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

metadata <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")

residuals="~/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Matrix_Processing/residuals_qned_counts.RData"
load(residuals)
residuals_ALL=residuals_qned_counts_objs[[3]]
colnames(residuals_ALL)=sapply(strsplit(colnames(residuals_ALL), "_"),"[",1)
samples=unique(c(colnames(residuals_ALL)))

print("Load enterotypes")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Enterotypes/enterotypes.RData")
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
ent_samples=enterotypes_samples_objs[[number]]
#ent_names=unique(ent_samples$Enterotype)

# GRM
print("Load GRM")
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/grm.RData")
grm=grm[intersect(rownames(grm), samples),intersect(rownames(grm), samples)]

this_metadata=metadata[metadata$RFID %in% rownames(grm),]
rownames(this_metadata)=this_metadata$RFID
ent_samples=ent_samples[ent_samples$sample %in% this_metadata$host_subject_id,]
this_metadata$Enterotype=ent_samples$Enterotype[match(this_metadata$host_subject_id,ent_samples$sample)]

grm_dist <- as.dist(1 - grm)
grm.pcoa <- ape::pcoa(grm_dist)
grm.pc <- data.frame(grm.pcoa$vectors)
grm.pc <- grm.pc [rownames(this_metadata),]

permanova.grm <- adonis2(grm.pc ~ Enterotype, data = this_metadata, permutations = 999)
p_value <- round(permanova.grm[1,5],digits=3)

pdf("/users/abaud/fmorillo/paper_figures/enterotype_population_structure.pdf")
myplot<-ggplot(grm.pc, aes(x=Axis.1, y=Axis.2, colour = as.factor(this_metadata$Enterotype))) +
  geom_point() + 
  ggtitle("Principal Coordinates of GRM by Enterotype") +
  scale_colour_manual(values = c("enterotype__1" = "#1f77b4", "enterotype__2" = "#ff7f0e")) +
  xlab("PCo1") +
  ylab("PCo2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.line = element_line(),
        plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=13),
        legend.title = element_blank(),
        legend.position="bottom") +
  annotate("label", 
           x = max(grm.pc$Axis.1), 
           y = max(grm.pc$Axis.2), 
           label = paste0("P = ",p_value),
           hjust = 1, vjust = 1,
           size = 5,
           label.size = 0.5,
           label.r = unit(0.2, "lines"),
           label.padding = unit(0.25, "lines"),
           fill = "white",
           color = "black",
           label.colour = "black")
print(myplot)
dev.off()
