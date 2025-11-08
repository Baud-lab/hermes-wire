#Functions
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

#cobalamin=c("K02188",
#            "K02224",
#            "K05936",
#            "K00595",
#            "K06042",
#            "K01771",
#            "K02231",
#            "K02233",
#            "K02227",
#            "K02226",
#            "K00768")
#cobalamin=data.frame(ko=cobalamin,module="Cobalamin Biosynthesis")

subset_id="ALL"
genera=c(
  #"g__Paramuribaculum",
  "g__Bacteroides",
  "g__Prevotella"
  #"g__CAG-485",
  #"g__CAG-873",
  #"g__Duncaniella",
  #"g__Acetatifactor"
  #"g__Phocaeicola",
  #"all"
  #"g__Helicobacter_D"
  #"g__UBA4372",
  #"g__Paraprevotella"
)

#genera=c(
#  #"g__Paramuribaculum",
#  "Bacteroides_low",
#  "Prevotella_low",
#  "Bacteroides_high",
#  "Prevotella_high"
#)

functions=c(
 # "eggNOG_OGs",
  "EC"
  #"KEGG_ko",
  #"KEGG_Module" 
  #"PFAMs",
  #"KEGG_Pathway",
  #"CAZy"
  #"COG_category",
  #"COG_pathway"
)

# Load functional annotations
load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_eggNOG_OGs.RData")
load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_matrices_v4.RData")

# Load the filtered functions
#load("/nfs/users/abaud/fmorillo/paper_figures/enrichment_t.test_final.RData")

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

## Load Heritability
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Heritability/heritability.RData")

filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
filtered_VCs$Heritability=as.numeric(filtered_VCs$var_Ad)*100
filtered_VCs <- filtered_VCs %>% 
  mutate( All_Sig = case_when(
    Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Heritable",
    TRUE ~ "Non-heritable"))
filtered_VCs$Genome=taxonomy$Genome[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]

#bac_mean=mean(filtered_VCs$Heritability[filtered_VCs$Genus=="g__Bacteroides"])
#prev_mean=mean(filtered_VCs$Heritability[filtered_VCs$Genus=="g__Prevotella"])
#
filtered_VCs$Genus[filtered_VCs$Genus=="g__Bacteroides" & filtered_VCs$All_Sig=="Heritable"]="Heritable Bacteroides"
filtered_VCs$Genus[filtered_VCs$Genus=="g__Bacteroides" & filtered_VCs$All_Sig!="Heritable"]="Non-heritable Bacteroides"
filtered_VCs$Genus[filtered_VCs$Genus=="g__Prevotella" & filtered_VCs$All_Sig=="Heritable"]="Heritable Prevotella"
filtered_VCs$Genus[filtered_VCs$Genus=="g__Prevotella" & filtered_VCs$All_Sig!="Heritable"]="Non-heritable Prevotella"

# Filter matrices based on the differential functions and selected genera
mean_matrix_total=NULL
func_matrix_total=NULL

ec_codes <- c(
  "EC4__3.1.1.11", "EC4__1.1.1.58", "EC4__5.5.1.4", "EC4__2.4.1.11",
  "EC4__4.1.1.19", "EC4__4.1.99.1", "EC4__3.5.1.1", "EC4__4.2.1.49",
  "EC4__2.1.2.5", "EC4__5.4.99.2", "EC4__5.1.99.1", "EC4__4.1.1.9",
  "EC4__2.5.1.55", "EC4__2.3.1.129", "EC4__3.5.1.108", "EC4__2.4.1.83",
  "EC4__1.1.1.290", "EC4__2.5.1.9", "EC4__2.4.2.19", "EC4__4.1.3.36",
  "EC4__1.3.5.1", "EC4__1.6.5.11", "EC4__1.1.1.37", "EC4__1.2.7.3"
  #"EC4__2.3.1.8", "EC4__1.10.3.14", "EC4__4.2.1.10", "EC4__1.15.1.2",
  #"EC4__2.1.1.192", "EC4__1.1.1.262", "EC4__2.7.8.13", "EC4__2.1.2.2", "EC4__6.3.2.9"
)
EC4_functions <- read_excel("paper_figures/EC4_functions.xlsx")
EC4_functions=EC4_functions[EC4_functions$EC4 %in% ec_codes,]
annotation=data.frame(#Origin = as.factor(EC4_functions$Origin),
                      Bioprocess = as.factor(EC4_functions$Bioprocess),
                      row.names=EC4_functions$EC4)
genera=c(
  "Heritable Bacteroides",
  "Non-heritable Bacteroides",
  "Heritable Prevotella",
  "Non-heritable Prevotella"
)
for (genus in genera){
  #genus=gsub("g__","",genus)
  #genus=gsub("-","",genus)
  print(genus)
  this_filtered_VCs=filtered_VCs[filtered_VCs$Genus==genus,]
  #this_filtered_VCs = this_filtered_VCs[order(this_filtered_VCs$Heritability, decreasing = T),]
  for(func in functions){
    # Filter functions and genera
    func_matrix=get(paste("func_matrix",func,sep="_"))
    #func_matrix=func_matrix[grepl("M00153",rownames(func_matrix)) | grepl("M00076",rownames(func_matrix)) | grepl("M00077",rownames(func_matrix)),]
    ### Fix it ###
    x=data.frame(genome=colnames(func_matrix)[grepl("^GCA_", colnames(func_matrix)) & !grepl("\\.1$", colnames(func_matrix))])
    x$fixed=paste0(x$genome,".1")
    y=data.frame(genome=colnames(func_matrix)[!colnames(func_matrix) %in% x$genome])
    y$fixed=y$genome
    all_genomes=rbind(x,y)
    colnames(func_matrix)=all_genomes$fixed[match(colnames(func_matrix),all_genomes$genome)]
    ###
    colnames(func_matrix)=taxonomy$Species[match(colnames(func_matrix),taxonomy$Genome)]
    func_matrix=func_matrix[,colnames(func_matrix) %in% this_filtered_VCs$trait,drop=F]
    func_matrix_total=cbind(func_matrix_total,func_matrix)
    #filtered_func=get(paste("all_cogs",genus,func,sep="_"))
    #filtered_func=filtered_func$COG[filtered_func$Significance!="Non-significant"]
    #if (length(filtered_func)>1){
   #clean_codes=c("M00153",
   #                "M00076",
   #                "M00077",
   #                "M00400",
   #                "M00403",
   #                "M00334",
   #                "M00607",
   #              "M00212"
   #                #"M00434","M00080","M00131"
   #                )
    # Original EC codes
    clean_codes <- sub("^EC\\s+", "EC4__", ec_codes)
    rownames(func_matrix)=paste0("EC4__",rownames(func_matrix))
    func_matrix=func_matrix[rownames(func_matrix) %in% clean_codes,]
    #setdiff(this_filtered_VCs$trait,colnames(func_matrix))
    #this_filtered_VCs=this_filtered_VCs[this_filtered_VCs$trait %in% colnames(func_matrix),]
    row_mean_matrix <- matrix(rowMeans(func_matrix), ncol = 1)
    colnames(row_mean_matrix)=genus
    rownames(row_mean_matrix)=rownames(func_matrix)
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00080"]="Molybdenum cofactor biosynthesis"
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00434"]="Formaldehyde assimilation"
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00131"]="Inositol phosphate metabolism"
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00076"]="Dermatan sulfate degradation"
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00077"]="Chondroitin sulfate degradation"
    #rownames(row_mean_matrix)[rownames(row_mean_matrix)=="M00153"]="Cytochrome bd ubiquinol oxidase"
    assign(paste("mean_matrix",genus,func,sep="_"),row_mean_matrix)
    mean_matrix_total=cbind(mean_matrix_total,row_mean_matrix)
    #if (func=="KEGG_ko"){
    #  rownames(func_matrix)=gsub("ko.","",rownames(func_matrix))
    #  annot_cobalamin=data.frame(ko=rownames(func_matrix))
    #  annot_cobalamin$`Cobalamin biosyntesis`=cobalamin$module[match(annot_cobalamin$ko,cobalamin$ko)]
    #  rownames(annot_cobalamin)=annot_cobalamin$ko
    #  annot_cobalamin=annot_cobalamin[,2,drop=F]
    #}
    #annot_heritability=data.frame(Heritability = this_filtered_VCs$Heritability)
    #colnames(func_matrix)=this_filtered_VCs$trait
    #rownames(annot_heritability) <- colnames(func_matrix)
  }
}

    
# Make the plot
#attr(func_matrix, "dimnames")[[2]] <- NULL
# Filter out non-mapped functions
#row_sums <- rowSums(mean_matrix_total)
#non_zero_columns <- row_sums != 0
#matrix_filt <- mean_matrix_total[non_zero_columns,]
#
#cazyme_classification <- read_excel("NextFlow/Microbiome_Profiling/step2/input/cazyme_classification.xlsx")
#cazyme_classification=cazyme_classification[cazyme_classification$CAZY %in% rownames(matrix_filt),]
#annotation=data.frame(Target = as.factor(cazyme_classification$Target))
#rownames(annotation)=rownames(matrix_filt)

total_bioprocesses=unique(c(annotation$Bioprocess))
colours <- metafolio::gg_color_hue(
  n = length(total_bioprocesses), 
  hue_min = 0,    # start at red
  hue_max = 330,  # nearly full circle
  c = 100,        # full chroma for vivid colors
  l = 60          # moderate luminance
)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]
names(colours) <- total_bioprocesses

# Define a named list as required by pheatmap
ann_colors <- list(
  Bioprocess = colours
  #Origin = setNames(c("coral", "lightgreen", "purple"), c("Correlation", "Enrichment - High heritability", "Enrichment - Low cohousing effects"))
)

#pdf("~/paper_figures/copies_of_ecs_bac_prev.pdf")
myplot=pheatmap(
  t(mean_matrix_total),
  cluster_rows = F,
  cluster_cols = T,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  show_rownames = T,
  legend = T,
  fontsize_row = 10,       # font size for row names
  fontsize_col = 10,       # font size for column names
  #cellwidth = 10,          # width of each cell
  #cellheight = 10,         # height of each cell
  legend_labels = "Mean copies",
  main = "Mean number of copies across mapped genomes",
  fontsize_main = 11,
)
print(myplot)
#dev.off()

write.table(mean_matrix_total,file="~/paper_figures/tables/ec4_copies_bac_prev.csv",sep=",",col.names=T,row.names=T,quote = F)

df=as.data.frame(mean_matrix_total)
df$diff_bac=df$`Heritable Bacteroides`-df$`Non-heritable Bacteroides`
df$diffprev=df$`Heritable Prevotella`-df$`Non-heritable Prevotella`
#df$diff=df$Bacteroides-df$Prevotella
df$func=rownames(df)
#df$target=cazyme_classification$Target[match(rownames(df),cazyme_classification$CAZY)]
View(df)

#grid.text(paste0(gsub("_"," ",func)," Counts"), x = 0.9, y = 0.5, gp = gpar(fontsize = 10))
# Save the filtered matrix
assign(paste("func_matrix",genus,func,sep="_"),func_matrix)
#}

#############################

### Genome homogeneity

library(tidyverse)

# Load your data: assume each row is a genome with genus, species, and EC columns
df <- ifelse(func_matrix_total > 0, 1, 0)

# Filter for Bacteroides and Prevotella genomes
mat_bac <- t(df[,grepl("Bacteroides",colnames(df))])
mat_pre <- t(df[,grepl("Prevotella",colnames(df))])

############## T-test

library(vegan)
# Jaccard distance
dist_bac <- vegdist(mat_bac, method = "jaccard")
dist_pre <- vegdist(mat_pre, method = "jaccard")

# Convert to similarity
sim_bac <- 1 - as.matrix(dist_bac)
sim_pre <- 1 - as.matrix(dist_pre)

mean_sim_bac <- mean(sim_bac[upper.tri(sim_bac)])
mean_sim_pre <- mean(sim_pre[upper.tri(sim_pre)])

dists_bac <- as.vector(dist_bac)
dists_pre <- as.vector(dist_pre)
t_tes_jaccard=t.test(dists_bac, dists_pre)

# Suppose dists_pre and dists_bac are your distance vectors
qqnorm(dists_pre, main = "Q‑Q Plot: Prevotella distances")
qqline(dists_pre, col = "red")

qqnorm(dists_bac, main = "Q‑Q Plot: Bacteroides distances")
qqline(dists_bac, col = "red")

shap_pre <- shapiro.test(dists_pre)
shap_bac <- shapiro.test(dists_bac)

print(shap_pre)
print(shap_bac)


########################### Jaccard test

library(jaccard)

# Filter out invariant
mat_bac <- mat_bac[, apply(mat_bac, 2, function(col) length(unique(col)) > 1)]
mat_pre <- mat_pre[, apply(mat_pre, 2, function(col) length(unique(col)) > 1)]

# Helper function to compute pairwise Jaccard tests with chosen method
compute_pairwise_jaccard <- function(mat, method = "bootstrap"){
  n <- nrow(mat)
  pvals <- matrix(NA, n, n, dimnames = list(rownames(mat), rownames(mat)))
  stats <- matrix(NA, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for(i in 1:(n - 1)){
    for(j in (i+1):n){
      out <- jaccard.test(mat[i,], mat[j,], method = method, verbose = FALSE)
      pvals[i,j] <- out$pvalue
      stats[i,j] <- out$statistics
      pvals[j,i] <- out$pvalue
      stats[j,i] <- out$statistics
    }
  }
  list(pvalues = pvals, statistics = stats)
}

res_bac_mca <- compute_pairwise_jaccard(mat_bac, method = "bootstrap")
res_pre_mca <- compute_pairwise_jaccard(mat_pre, method = "bootstrap")

pvec_bac_mca <- na.omit(as.vector(res_bac_mca$pvalues[upper.tri(res_bac_mca$pvalues)]))
pvec_pre_mca <- na.omit(as.vector(res_pre_mca$pvalues[upper.tri(res_pre_mca$pvalues)]))

qvec_bac_mca <- p.adjust(pvec_bac_mca, method = "BH")  # BH-adjusted p-values
qvec_pre_mca <- p.adjust(pvec_pre_mca, method = "BH")

mean(qvec_bac_mca < 0.05)
mean(qvec_pre_mca < 0.05)