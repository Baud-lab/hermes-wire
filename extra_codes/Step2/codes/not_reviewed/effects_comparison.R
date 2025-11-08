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

# Read taxonomy
print("Reading taxonomy info")
input_taxonomy="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt"
if (input_taxonomy != "") {
  taxonomy<-read_tab_file(input_taxonomy, FALSE)
  if (grepl("GTDB", input_taxonomy)) {
    colnames(taxonomy)=c("Genome","Taxonomy")
    taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
    taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
  } else {
    if (grepl("Greengenes2", input_taxonomy)){
      colnames(taxonomy) <- c("Genome","Taxonomy","Confidence","GTDB")
      taxonomy$Taxonomy<-gsub("\\s(?!.*\\s)", "_", as.character(taxonomy$Taxonomy), perl = TRUE)
      taxonomy$Taxonomy<-gsub(" ", "", as.character(taxonomy$Taxonomy), perl = TRUE)
    } else {
      stop("The taxonomy provided should be GTDB for shotgun data or Greengenes2 for 16S")
    }
  }
  taxonomy$Phylum<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",2)
  taxonomy$Class<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",3)
  taxonomy$Order<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",4)
  taxonomy$Family<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",5)
  taxonomy$Genus<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",6)
  taxonomy$Species<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",7)
  taxonomy$ASV=paste0("asv__",rownames(taxonomy))
} else {
  stop("No taxonomy file provided for the main dataset")
}

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/old_Outputs/old_maternal/Output_Maternal/Heritability/heritability.RData")
subset_id="ALL"
filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
#filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]

filtered_VCs$Host_Aggregate_Genetics=100*filtered_VCs$var_Ad
filtered_VCs$Cohousing=100*filtered_VCs$var_C
filtered_VCs$Maternal_Effects=100*filtered_VCs$var_M
filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]

new_df <- filtered_VCs %>%
  pivot_longer(cols = c("Host_Aggregate_Genetics", "Cohousing", "Maternal_Effects"), 
               names_to = "Effect", 
               values_to = "Value")

max_y <- max(new_df[new_df$Phylum=="p__Bacteroidota",]$Value)
buffer <- 2 # Add a buffer to avoid overlapping
y_positions <- max_y + c(1, 1, 3) * buffer

myplot<-ggplot(new_df[new_df$Phylum=="p__Bacteroidota",], aes(x=reorder(Effect,Value,FUN = median), y=Value)) +
  geom_boxplot() +
  xlab("") +
  ylab("Effect size (%)") +
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
  geom_signif(comparisons = list(c("Cohousing","Host_Aggregate_Genetics"),c("Cohousing","Maternal_Effects"),c("Host_Aggregate_Genetics","Maternal_Effects")),
              map_signif_level = TRUE, y_position = y_positions)
print(myplot)

#Genus

max_y <- max(new_df[new_df$Genus=="g__CAG-485",]$Value)
buffer <- 2 # Add a buffer to avoid overlapping
y_positions <- max_y + c(1, 1, 3) * buffer

myplot<-ggplot(new_df[new_df$Genus=="g__CAG-485",], aes(x=reorder(Effect,Value,FUN = median), y=Value)) +
  geom_boxplot() +
  xlab("") +
  ylab("Effect size (%)") +
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
  geom_signif(comparisons = list(c("Cohousing","Maternal_Effects"),c("Host_Aggregate_Genetics","Maternal_Effects"),c("Cohousing","Host_Aggregate_Genetics")),
              map_signif_level = TRUE, y_position = y_positions)
print(myplot)
