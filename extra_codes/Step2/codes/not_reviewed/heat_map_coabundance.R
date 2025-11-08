# Load correlation matrix
SparCC_Output_ALL <- read.csv("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/SparCC_Output_ALL.csv")
rownames(SparCC_Output_ALL)=SparCC_Output_ALL$Genome_ID
SparCC_Output_ALL=SparCC_Output_ALL[,-1]
matrix=as.matrix(SparCC_Output_ALL)


# Read taxonomy
input_taxonomy="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt"
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
print("Reading taxonomy info")
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
} else {
  stop("No taxonomy file provided for the main dataset")
}


# Collapse into species
colnames(matrix)=taxonomy$Species[match(colnames(matrix), taxonomy$Genome)]
rownames(matrix)=taxonomy$Species[match(rownames(matrix), taxonomy$Genome)]

#Load heritability
load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/new_29-01-24/Heritability/heritability.RData")
all_VCs_filt=all_VCs_full[grepl("ALL",all_VCs_full$subset_id),]
herit=100*all_VCs_filt$var_Ad
names(herit)=all_VCs_filt$trait

# Filter by genus
species="s__Bacteroides"
filt_matrix=matrix[grepl(species,rownames(matrix)),grepl(species,colnames(matrix))]
filt_herit=herit[names(herit) %in% rownames(filt_matrix)]


# Heat map
heatmap(filt_matrix, 
        symm = TRUE)

ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Correlation Heatmap with Extra Value") +
  theme_minimal()
