####

library(readxl)
library(rhdf5, lib.loc = "/nfs/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(biomformat)
library(xfun)
library(ggplot2)

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

## Define samples
#print("Defining samples based on the latest SS analysis")
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2//Output_Tax_Mat_Net/Cluster_Analyses/residuals_qned_counts_new.RData")
#samples=colnames(residuals_qned_counts_objs[[1]])
#samples=sapply(strsplit(samples, "_"), "[",1)

# Load stats table
print("Load stats table")
#stats <- read.csv("/nfs/users/abaud/fmorillo/P50/microbiome/output/Woltka/stats.tsv")
stats <- read.delim("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/statistics.txt")
stats$prop_gtdb=100*stats$Mapped_Reads/stats$Non_Host_Reads
stats$Sample=sapply(strsplit(stats$Sample, "_"), "[",1)
#stats=stats[stats$Sample %in% samples,]

## Load non-host reads from qiita
print("Load non-host reads from qiita")
qiita_after_filters <- read_excel("/nfs/users/abaud/fmorillo/P50/microbiome/output/Woltka/qiita_after_filters.xlsx")
qiita_after_filters=qiita_after_filters[!grepl("BLANK",qiita_after_filters$filename),]
qiita_after_filters$type=sapply(strsplit(qiita_after_filters$filename, "_"), "[",4)
qiita_after_filters$filename=sapply(strsplit(qiita_after_filters$filename, "_"), "[",1)
qiita_after_filters=qiita_after_filters[!grepl("R2",qiita_after_filters$type),]
for(i in 1:length(qiita_after_filters$filename)){
  if(!nchar(qiita_after_filters$filename[i])>=10){qiita_after_filters$filename[i] <- paste0("000", qiita_after_filters$filename[i])}
}
qiita_after_filters=qiita_after_filters[qiita_after_filters$filename %in% stats$Sample,]
qiita_after_filters <- qiita_after_filters[match(stats$Sample, qiita_after_filters$filename), ]
stats$available_reads_qiita=qiita_after_filters$reads
#stats$available_reads_qiita_2x=2*qiita_after_filters$reads

# Load genome counts table
print("Load Woltka genome counts table")
input_matrix="/users/abaud/fmorillo/P50/microbiome/output/Woltka/192967_none.biom"
counts_table<-as.matrix(biom_data(read_biom(input_matrix)))
#### Customized!!!!
counts_table=counts_table[,!grepl("BLANK",colnames(counts_table))]
ids <- sapply(strsplit(colnames(counts_table), "\\."), function(x) tail(x, 1))
if(!all(nchar(ids)>=10)){
  for(i in 1:length(ids)){
    if(!nchar(ids[i])>=10){ids[i] <- paste0("000", ids[i])}
  }
}
colnames(counts_table)<-ids
counts_table=counts_table[,-c(which(colnames(counts_table)=="00077E9A45"), which(colnames(counts_table)=="00077E7BC3"))]

woltka_depth=data.frame(sample=colnames(counts_table),reads=colSums(counts_table))
woltka_depth <- woltka_depth[match(stats$Sample, woltka_depth$sample), ]
stats$aligned_to_wol=woltka_depth$reads/2
stats$prop_wol_qiita=100*stats$aligned_to_wol/stats$available_reads_qiita

#load("/users/abaud/fmorillo/P50/microbiome/output/Woltka/stats.RData")

gtdb_stats=data.frame(value=stats$prop_gtdb,Method="GTDB with Kaiju")
wol_stats=data.frame(value=stats$prop_wol_qiita,Method="WoL with Woltka")


# Filter out low abundance (<0.01%)

#GTDB
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Matrix_Processing/filtered_traits_matrix_noprev.RData")
#matrix=as.matrix(matrix.abundance_filtered)
matrix=filtered_trait_objs[[1]]
species_gtdb=rownames(matrix)
colnames(matrix)<-sapply(strsplit(colnames(matrix), "_"), "[",1)
rownames(stats)<-stats$Sample
stats_filt=stats[rownames(stats) %in% colnames(matrix),]
stats_filt=stats[colnames(matrix),]
stats_filt$Mapped_Reads_After_Abundance_Filtering_gtdb=colSums(matrix)/2
stats_filt$Clear_Mapping_Ratio_gtdb=100*(stats_filt$Mapped_Reads_After_Abundance_Filtering_gtdb)/stats_filt$Non_Host_Reads

#Woltka
print("Filter WoL species based on abundance")
print("Load samples to be removed")
exclude_file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/remove_ids_Shallow.txt"
if (exclude_file != "") {
  print("Reading files with ids to be removed")
  exclude_ids<-read_csv_file(exclude_file)
  colnames(exclude_ids)[1]="Files"
}
exclude_ids<-c(exclude_ids[,1])
exclude_ids<-sapply(strsplit(exclude_ids, "_"), "[",1)
print(paste0("Counts table dimensions: ",dim(counts_table)," (Before removing bad samples and samples with initial low depths)"))
matrix.sample_filtered <- counts_table[,!colnames(counts_table) %in% exclude_ids]
print(paste0("Counts table dimensions: ",dim(matrix.sample_filtered)," (After removing bad samples and samples with initial low depths)"))
# Identify low abundant traits
print("Identify and remove low abundant traits")
## Re-scale the table in relative values, make sure the sum of each row is 1
matrix.relative <- sweep(matrix.sample_filtered,2,colSums(matrix.sample_filtered),"/")
## Calculate the mean relative abundance of each trait
spe.meanAbun <- apply(matrix.relative, 1, mean, na.rm=TRUE)
## Filter based on the mean relative abundance
min_abundance=0.0001
spe.keep <- rownames(matrix.relative)[spe.meanAbun > min_abundance]
species_wol=spe.keep
print(paste0("Counts table dimensions: ",dim(matrix.sample_filtered)," (Before filtering out taxa)"))
matrix.abundance_filtered <- matrix.sample_filtered[spe.keep,]
print(paste0("Counts table dimensions: ",dim(matrix.abundance_filtered)," (After filtering out taxa)"))


counts_table_filt=matrix.abundance_filtered[,colnames(matrix.abundance_filtered) %in% colnames(matrix)]
counts_table_filt=matrix.abundance_filtered[,colnames(matrix)]
stats_filt$Mapped_Reads_After_Abundance_Filtering_wol=colSums(counts_table_filt)/2
stats_filt$Clear_Mapping_Ratio_wol=100*stats_filt$Mapped_Reads_After_Abundance_Filtering_wol/stats_filt$available_reads_qiita

print("Mean mapping ratio GTDB (Clear DB)")
mean_clear_mr_gtdb=mean(stats_filt$Clear_Mapping_Ratio_gtdb)
print(mean_clear_mr_gtdb)
#54.55
print("Mean mapping ratio WoL (Clear DB)")
mean_clear_mr_wol=mean(stats_filt$Clear_Mapping_Ratio_wol)
print(mean_clear_mr_wol)
#5.34

print("Saving stats")
save(stats_filt,
     gtdb_stats,
     wol_stats,
     file="/users/abaud/fmorillo/paper_figures/funnel_gtdb_vs_wol.RData"
)

## Common species
metadata <- read.delim("/nfs/users/abaud/fmorillo/P50/microbiome/output/Woltka/metadata.tsv")
rownames(counts_table_filt)=metadata$assembly_accession[match(rownames(counts_table_filt),metadata$X.genome)]
species_wol <- unique(rownames(counts_table_filt))
print("Number of species - WoL")
print(length(species_wol))
print("Saving species wol")
save(species_wol,
     file="/users/abaud/fmorillo/paper_figures/species_wol.RData")
#234
species_gtdb <- unique(rownames(matrix))
print("Number of species - GTDB")
print(length(species_gtdb))
#1181
print("Saving species gtdb")
save(species_gtdb,
     file="/users/abaud/fmorillo/paper_figures/species_gtdb.RData")
common=intersect(species_gtdb,species_wol)
print("Number of shared species - WoL and GTDB")
print(length(common))
#40

### Plot
db=rbind(gtdb_stats,wol_stats)
db=na.omit(db)

# Plot the histogram
pdf("~/paper_figures/gtdb_vs_wol.pdf")
ggplot(db, aes(x = value, fill = Method)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "",
       x = "Mapping ratio (%)",
       y = "Frequency",
       fill = "Profiling Method") +
  scale_fill_manual(values = c("darkblue", "lightblue")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),        # Remove grid lines
        axis.line = element_line(color = "grey"),
        legend.position = "bottom",
        legend.title.align = 0.5)
dev.off()

#########################################################################################################################################################

# Load proportion aligned from HumanN3:
humann=read.delim("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_alignments/profiles/humann_mapped_reads.tsv")
ids <- sapply(strsplit(humann$sample, "\\."), function(x) tail(x, 1))
if(!all(nchar(ids)>=10)){
  for(i in 1:length(ids)){
    if(!nchar(ids[i])>=10){ids[i] <- paste0("000", ids[i])}
  }
}
humann$sample<-ids
stats=stats[stats$sample %in% humann$sample,]
humann <- humann[match(stats$sample, humann$sample), ]
stats$prop_humann=humann$prop_mapped


boxplot(stats$prop_gtdb, stats$prop_wol_qiita, stats$prop_humann, names=c("GTDB_207","WoL_2","Humann_3"), xlab="Database", ylab="Mapping ratio (%)")
#summary(stats$prop_gtdb)
#summary(stats$prop_wol_qiita)

###



stats_no_out=stats[stats$prop_wol_qiita <= 100,]
boxplot(stats_no_out$prop_gtdb, stats_no_out$prop_wol_qiita, names=c("GTDB","WoL"), xlab="Database", ylab="Mapping ratio (%)")

boxplot(stats$available_reads, stats$available_reads_qiita)


#######

#Prromenade
EC4 <- read.delim("~/P50/microbiome/output/microbiome_profiles/functional_profile/prromenade/step1/Counts_Functional/EC4.out")
EC4_matrix=as.matrix(EC4[,-1])
colnames(EC4_matrix)=gsub("^X","",colnames(EC4_matrix))
colnames(EC4_matrix)=sapply(strsplit(colnames(EC4_matrix), "_"), "[",1)
EC4_matrix=EC4_matrix[,colnames(EC4_matrix) %in% stats$sample]
EC4_matrix=EC4_matrix[,match(stats$sample, colnames(EC4_matrix))]
EC4_depth=colSums(EC4_matrix)
stats$Prromenade_EC4_mapped=EC4_depth
stats$prop_Prromenade_EC4=100*stats$Prromenade_EC4_mapped/stats$available_reads
boxplot(stats$prop_Prromenade_EC4)
summary(stats$prop_Prromenade_EC4)

#WoL

#
#
############# Mice paper
## Load non-host reads from qiita
#qiita_after_filters <- read_excel("/users/abaud/fmorillo/P50/microbiome/output/Woltka/mice_paper/zengler_mice_raw_reads.xlsx")
#qiita_after_filters=qiita_after_filters[!grepl("etoh",qiita_after_filters$filename),]
#qiita_after_filters=qiita_after_filters[!grepl("ctrl",qiita_after_filters$filename),]
#qiita_after_filters$type=sapply(strsplit(qiita_after_filters$filename, "_"), "[",3)
#qiita_after_filters$filename=sapply(strsplit(qiita_after_filters$filename, "_"), "[",1)
#qiita_after_filters=qiita_after_filters[!grepl("R2",qiita_after_filters$type),]
##for(i in 1:length(qiita_after_filters$filename)){
##  if(!nchar(qiita_after_filters$filename[i])>=10){qiita_after_filters$filename[i] <- paste0("000", qiita_after_filters$filename[i])}
##}
#stats=data.frame(sample=qiita_after_filters$filename,available_reads_qiita_2x=2*qiita_after_filters$reads)
#
#
## Load genome counts table
#input_matrix="/users/abaud/fmorillo/P50/microbiome/output/Woltka/mice_paper/zengler_mice_genomes.biom"
#counts_table<-as.matrix(biom_data(read_biom(input_matrix)))
#counts_table=counts_table[,!grepl("etoh",colnames(counts_table))]
#counts_table=counts_table[,!grepl("ctrl",colnames(counts_table))]
#ids=gsub("13052.","",colnames(counts_table))
#colnames(counts_table)=ids
#
#woltka_depth=data.frame(sample=colnames(counts_table),reads=colSums(counts_table))
#woltka_depth <- woltka_depth[match(stats$sample, woltka_depth$sample), ]
#stats$aligned_to_wol=woltka_depth$reads
#stats$prop_wol_qiita=100*stats$aligned_to_wol/stats$available_reads_qiita_2x
#
#boxplot(stats$prop_wol_qiita,ylab="Mapping ratio (%)")
#boxplot(stats$available_reads_qiita_2x,stats$aligned_to_wol,names=c("Raw reads","Mapped to WoL2"),ylab="Number of reads (millions)")
