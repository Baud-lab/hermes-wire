#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(biomformat))
suppressMessages(library(MASS))
suppressMessages(library(rhdf5))
suppressMessages(library(vegan))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(metafolio))
suppressMessages(library(philr))
suppressMessages(library(argparse))
suppressMessages(library(arrow))
suppressMessages(library(xfun))
suppressMessages(library(iNEXT.3D))
suppressMessages(library(MASS))
suppressMessages(library(doParallel))
suppressMessages(library(ggsignif))

# -----------------------------------------------------------------------------
# Functions

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

# 2. Assign covariates to subsets
assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames <- str_split(x[2], pattern=",", simplify = TRUE) 
  assign(paste('covariates',subset_id,sep='_'), covnames, envir = .GlobalEnv)
}


# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Modules
parser$add_argument("-mod_ana", "--module_analysis", type="character", help="Inclusion of module: Analysis",default="Pooled")
parser$add_argument("-mod_herit", "--module_heritability", type="character", help="Inclusion of module: Heritability Analysis",default="NO")
parser$add_argument("-mod_gwas", "--module_gwas", type="character", help="Inclusion of module: GWAS Analysis",default="NO")
## Matrix processing
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-exclude", "--exclude_ids", type="character", help="File with list of ids to be removed (optional)")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-lib", "--input_libsize", type="character", help="File with library depth info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the sample",default="host_subject_id")
parser$add_argument("-mid", "--min_initial_depth", type="integer", help="Minimum number of reads accepted", default=0)
parser$add_argument("-ma", "--min_abundance", type="double", help="Minimum mean relative abundance accepted", default=0.0001)
parser$add_argument("-mfd", "--min_final_depth", type="integer", help="Minimum number of reads accepted after filtering noise", default=30000)
parser$add_argument("-mat", "--sample_material", type="character", help="Define the type of material used", default="Cecal")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-mp", "--min_prevalence", type="double", help="Minimum prevalence accepted", default=0.5)
## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Modules
module_analysis<-args$module_analysis
module_heritability<-args$module_heritability
if (module_heritability != "YES" & module_heritability != "NO") {
  stop("Option for module heritability should be either 'YES' or 'NO' in capital letters", "\n")
}
module_gwas<-args$module_gwas
if (module_gwas != "YES" & module_gwas != "NO") {
  stop("Option for module gwas should be either 'YES' or 'NO' in capital letters", "\n")
}
## Matrix processing
input_matrix<-args$input_matrix
exclude_file<-args$exclude_ids
input_meta<-args$input_meta
input_libsize<-args$input_libsize
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file", "\n")
}
min_initial_depth<-args$min_initial_depth
min_abundance<-args$min_abundance
min_final_depth<-args$min_final_depth
min_prevalence<-args$min_prevalence
Material=args$sample_material
if (Material == "" & (module_heritability != "NO" || module_gwas != "NO")) {
  stop("Please inform the sample type under investigation if you want to run any Host Genetic Effects Analyses. It should correspond to the tissue type where the metadata was collected from (i.e. Cecal, Fecal, etc.)", "\n")
}
main_method=args$method

## Taxonomic
input_taxonomy<-args$input_taxonomy
ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}

# -----------------------------------------------------------------------------


### Step 1 ####
print("Step 1: Read input files and arguments and do the custom removal of bad samples")

# Read counts table and remove extension
if (input_matrix != "") {
  print("Reading counts table")
  if (tail(strsplit(input_matrix, "\\.")[[1]], 1)=="biom"){
    counts_table=as.matrix(biom_data(read_biom(input_matrix)))
    colnames(counts_table) <- str_sub(colnames(counts_table), - 10, - 1)
  } else {
    counts_table<-read_tab_file(input_matrix, TRUE)
    colnames(counts_table)<-tools::file_path_sans_ext(colnames(counts_table))
    row.names(counts_table)<-counts_table$Genome_id
    counts_table$Genome_id<-NULL 
  }
} else {
  stop("No counts table provided")
}

# Read files with ids to be removed
if (exclude_file != "") {
  print("Reading files with ids to be removed")
  exclude_ids<-read_csv_file(exclude_file)
  colnames(exclude_ids)[1]="Files"
}

# Extract metadata info
if (input_meta != "") {
  print("Extracting metadata info")
  metadata<-read_tab_file(input_meta, TRUE)
  if (ncol(metadata)==1){
    metadata<-read_csv_file(input_meta)
  }
  if (length(which(colnames(metadata)==sample_identifier)) == 0) {
    stop("The metadata file has no column with the same name of the sample identifier informed")
  }
} else {
  stop("No metadata file provided")
}

# Read library depth info
print("Reading library depth info")
if (main_method=="16S"){
  depths<-data.frame(Sample=metadata[[sample_identifier]],Raw_Reads=metadata$Depth_16S)
} else {
  if (input_libsize != "") {
    depths<-read_tab_file(input_libsize, TRUE)
  } else {
    stop("No library depth file provided")
  }
}

# Read covariates info
if (list_of_covariates != "") {
  print("Reading covariates and covariate types info")
  covs_data<-read_tab_file(list_of_covariates, TRUE)
  if (ncol(covs_data)==1){
    covs_data<-read_csv_file(list_of_covariates)
  }
  subset_ids <- covs_data$subset_id
  if (nrow(covs_data)==0) {
    stop("No subset was provided. If you only worked with one cohort/Study, please include subset 'ALL' on the input file covariates.txt and assign known covariates to it")
  }
  for (i in 1:nrow(covs_data)) {
    if (covs_data$covariates[i] == "") {
      stop(paste0("Subset ",covs_data$subset_id[i]," has no covariate assigned. Please inform at least one known covariate for this subset on the input file covariates.txt and remember that the name of the covariate should correspond with the name of one of the columns of your matadata file. If more than one covariate are assigned for this subset, remember to list them separated by comma with no spaces.")) 
    }
  }
  for (i in 1:nrow(covs_data)) {
    if (nrow(covs_data) > 1 & covs_data$subset_id[i] == "ALL" & !grepl("Study",covs_data$covariates[i])) {
      covs_data$covariates[i]=paste(covs_data$covariates[i],"Study",sep=",")
      print(paste0("Warning: Subset ",covs_data$subset_id[i]," had no covariate 'Study' assigned to it on the input file covariates.txt. Given that other subsets were listed, this covariate was automaticaly included.")) 
    }
  }
  apply(covs_data, 1,  assignCovariate) ## please check: I changed hard coded 1 for subset_id
  covnames <- unique(c(str_split(covs_data$covariates[covs_data$subset_id=="ALL"], pattern=",", simplify = TRUE)))
  total_covs=c()
  for (subset_id in subset_ids[-which(subset_ids=="ALL")]){
    covariate=get(paste("covariates",subset_id,sep="_"))
    total_covs=unique(c(total_covs,covariate))
  }
  missing_covs <- total_covs[!(total_covs %in% covnames)]
  if (length(missing_covs) == 0) {
    print("All elements in the subsets are present in the subset 'ALL'.")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
  if (list_of_covariate_types != "") {
    covs_types<-read_tab_file(list_of_covariate_types, TRUE)
    if (ncol(covs_types)==1){
      covs_types<-read_csv_file(list_of_covariate_types)
    }
    if (nrow(covs_types)==0) {
      stop("No covariate was provided to the covariate types file. Please, include at least one known covariate on the input file covariates.classification.txt")
    }
    for (i in 1:nrow(covs_types)) {
      if (covs_types$classification[i] == "") {
        stop(paste0("Subset ",covs_types$covariate[i]," has no classification assigned. Please, inform on the input file covariates.classification.txt if the covariate is 'categorical' or 'continuous'.")) 
      }
    }
    for (i in 1:length(covnames)) {
      if (length(which(covs_types$covariate==covnames[i])) == 0) {
        stop(paste0("No classification was assigned to covariate: ",covnames[i],". Please, inform on the input file covariates.classification.txt if this covariate is 'categorical' or 'continuous'"))
      }
    }
    for (subset_id in subset_ids) {
      covariates=get(paste("covariates",subset_id,sep="_"))
      covs=c()
      for (i in 1:length(covariates)){
        cat=covs_types$classification[match(covariates[i],covs_types$covariate)]
        if (cat=="categorical"){
          covs=c(covs,covariates[i])
        }
      }
      assign(paste('covariates',subset_id,sep='_'), covs)
    }
  }  else {
    stop("No covariate types file provided")
  }
} else {
  stop("No covariates file provided")
}

# Read taxonomy
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
      taxonomy$ASV=paste0("asv__",rownames(taxonomy))
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

# Read taxonomic ranks
print("Reading taxonomic ranks")
ranks<-na.omit(c((sapply(strsplit(ranks, ","),"[",1)),
                 (sapply(strsplit(ranks, ","), "[",2)),
                 (sapply(strsplit(ranks, ","), "[",3)),
                 (sapply(strsplit(ranks, ","), "[",4)),
                 (sapply(strsplit(ranks, ","), "[",5)),
                 (sapply(strsplit(ranks, ","), "[",6)),
                 (sapply(strsplit(ranks, ","), "[",7))))

# Remove controls to be used only on Beta Diversity Analyses
print("Isolating controls")
# Subset metadata info
motch = match(colnames(counts_table),metadata[[sample_identifier]]) 
metadata = metadata[motch,]
# Isolate controls
metadata_controls=metadata[grepl("Positive Control",metadata$Material) | grepl("Negative Control",metadata$Material),]
metadata=metadata[!metadata[[sample_identifier]] %in% metadata_controls[[sample_identifier]],]
matrix.sample_controls=counts_table[,match(metadata_controls[[sample_identifier]],colnames(counts_table))]
counts_table=counts_table[,!colnames(counts_table) %in% colnames(matrix.sample_controls)]
save(matrix.sample_controls,
     metadata_controls,
     file = 'controls.RData') 

# Filter by Material (This is important to keep only one sample per host for host genetic analyses)
print("Filtering by material type")
if (module_heritability != "NO" || module_gwas != "NO"){
  metadata = metadata[grepl(Material,metadata$Material),]
  counts_table = counts_table[,match(metadata[[sample_identifier]],colnames(counts_table))]
}

# Identify low covered samples
print("Identifying low covered sample")
low_depth_ids<-depths[depths$"Raw_Reads"<min_initial_depth, ]$Sample

# Remove bad samples from a custom list and/or with low depth
print("Removing bad samples from a custom list and/or with low depth")
if (exclude_file != "") {
  exclude_ids<-c(exclude_ids[,1], low_depth_ids)
} else {
  exclude_ids<-c(low_depth_ids)
}
print(paste0("Counts table dimensions: ",dim(counts_table)," (Before removing bad samples and samples with initial low depths)"))
matrix.sample_filtered <- counts_table[,!colnames(counts_table) %in% exclude_ids]
print(paste0("Counts table dimensions: ",dim(matrix.sample_filtered)," (After removing bad samples and samples with initial low depths)"))

# Subset metadata info
print("Subseting metadata info")
motch = match(colnames(matrix.sample_filtered),metadata[[sample_identifier]]) 
if (any(is.na(motch))) {
  missing_samples <- colnames(matrix.sample_filtered)[!(colnames(matrix.sample_filtered) %in% metadata[[sample_identifier]])]
  stop(paste("These samples were not found in metadata:", paste(missing_samples, collapse = ", ")))
}
metadata = metadata[motch,]

### End of Step 1 ###

### Step 2: Taxa filtering based on relative abundance (Important to remove noise) ###
print("Step 2: Filter low abundance taxa")

if (module_analysis=="Pooled"){
  studies=subset_ids[-which(subset_ids=="ALL")]
  subset_ids="ALL"
}

filtered_trait_objs <- list()
for (subset_id in subset_ids) {
  if(subset_id != "ALL") {
    matrix = matrix.sample_filtered[,which(metadata[,"Study"] == subset_id)]  
  } else {
    samp.keep=metadata[[sample_identifier]][metadata$Study %in% studies]
    matrix = matrix.sample_filtered[,samp.keep]
  }
  print(paste0(subset_id,": ",dim(matrix)," (Before filtering out taxa)"))
  # Identify low abundant traits
  ## Re-scale the table in relative values, make sure the sum of each row is 1
  matrix.relative <- sweep(matrix,2,colSums(matrix),"/")
  ## Calculate the mean relative abundance of each trait
  spe.meanAbun <- apply(matrix.relative, 1, mean, na.rm=TRUE)
  ## Filter based on the mean relative abundance
  spe.keep <- rownames(matrix.relative)[spe.meanAbun > min_abundance]
  matrix.abundance_filtered <- matrix[spe.keep,]
  print(paste0(subset_id,": ",dim(matrix.abundance_filtered)," (After filtering out taxa)"))
  if (subset_id == "ALL") {
    pdf('depth_funnel.pdf')
    if (main_method !="16S"){
      # Make the depth funnel plot
      print("Making the depth funnel plot")
      depth <- colSums(matrix.abundance_filtered)
      depths=depths[depths$Sample %in% colnames(matrix.abundance_filtered),]
      depths=depths[match(colnames(matrix.abundance_filtered), depths$Sample),]
      depths$Mapped_Reads_After_Abundance_Filtering=depth/2
      depths$Clear_Mapping_Ratio=100*depths$Mapped_Reads_After_Abundance_Filtering/depths$Non_Host_Reads
      depths$Raw_Reads=depths$Raw_Reads/1000
      depths$Reads_After_Trimming=depths$Reads_After_Trimming/1000
      depths$Non_Host_Reads=depths$Non_Host_Reads/1000
      depths$Mapped_Reads=depths$Mapped_Reads/1000
      depths$Mapped_Reads_After_Abundance_Filtering=depths$Mapped_Reads_After_Abundance_Filtering/1000
      depths_long <- tidyr::pivot_longer(depths, cols = c(Raw_Reads, Reads_After_Trimming, Non_Host_Reads, Mapped_Reads, Mapped_Reads_After_Abundance_Filtering),
                                         names_to = "Funnel Step", values_to = "Depth (Thousands of Pairs of Reads)")
      order1 <- c("Raw_Reads", "Reads_After_Trimming", "Non_Host_Reads", "Mapped_Reads", "Mapped_Reads_After_Abundance_Filtering")
      order2 <- c("Initial_Rate","Surviving_Rate_Trimming", "Surviving_Rate_Alignment", "Raw_Mapping_Ratio", "Clear_Mapping_Ratio")
      labels <- c("Raw Reads", "Reads After Trimming", "Non-Host Reads", "Mapped Reads", "Mapped Reads After Abundance Filtering")
      mean_values <- depths %>%
        summarise(
          Initial_Rate=100,
          Surviving_Rate_Trimming = mean(Surviving_Rate_Trimming),
          Surviving_Rate_Alignment = mean(Surviving_Rate_Alignment),
          Raw_Mapping_Ratio = mean(Raw_Mapping_Ratio),
          Clear_Mapping_Ratio = mean(Clear_Mapping_Ratio)
        ) %>%
        pivot_longer(cols = everything(), names_to = "Surviving Rates", values_to = "Mean Surviving Rate")
      mean_values$`Mean Surviving Rate`=round(mean_values$`Mean Surviving Rate`,digits = 2)
      scale=mean(depths$Raw_Reads)/mean(mean_values$`Mean Surviving Rate`)
      myplot<-ggplot(depths_long, aes(x = factor(`Funnel Step`, levels = order1), y = `Depth (Thousands of Pairs of Reads)`)) +
        geom_boxplot() +
        geom_label(data = mean_values,  
                   aes(x = as.numeric(factor(`Surviving Rates`, levels = order2)),
                       y = `Mean Surviving Rate` *scale, 
                       label = paste0(`Mean Surviving Rate`, "%"), fill = factor(`Surviving Rates`, levels = order2)),
                   vjust = -0.5, size = 4) +
        scale_x_discrete(labels = labels) +
        scale_fill_manual(values = c("Initial_Rate" = "lightcoral","Surviving_Rate_Trimming" = "beige", "Surviving_Rate_Alignment" = "lightgreen", "Raw_Mapping_Ratio" = "lightblue", "Clear_Mapping_Ratio" = "lightpink"),
                          labels = c("Starting Point", "Surviving: After Trimming", "Survinging: Non-Host", "Raw Mapping Ratio", "Clean Mapping Ratio"),
                          name = "Mean Filtering Rates") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 13),
              legend.text=element_text(size=11),
              legend.title=element_text(size=13),) +
        guides(fill = guide_legend(override.aes = list(color = NA)), 
               color = "none", 
               shape = "none") +
        labs(x = NULL)
      print(myplot)
    }
    dev.off()
    # Make asymptotic richness plots
    print("Making asymptotic richness plots")
    pdf('asymptotic_richness_vs_depth.pdf')
    matrix=round(matrix.abundance_filtered, digits = 0)
    depth <- colSums(matrix)
    names(depth) <- colnames(matrix)
    # Calculate Richness in parallel
    cl <- makeCluster(detectCores())
    print(cl)
    registerDoParallel(cl)
    alpha_table <- foreach(i = 1:ncol(matrix), .combine = rbind) %dopar% {
      col_data <- matrix[,i]
      names(col_data)=rownames(matrix)
      iNEXT.3D::AO3D(data=col_data,
                     diversity = "TD",
                     q = 0,
                     datatype = "abundance",
                     nboot = 50,
                     conf = 0.95,
                     method = "Asymptotic")
    }
    alpha_table$Assemblage=colnames(matrix)
    alpha_table$depth=depth
    stopImplicitCluster()
    # Plot asymptotic curve
    plot(alpha_table$depth,alpha_table$qTD,
         xlab="Depth (Reads per sample)",
         ylab="Asymptotic Richness (Expected number of species per sample)",
         ylim=c(0,max(alpha_table$qTD)+5))
    title(main = paste0("Method: ",main_method), line = 2.5)
    abline(v = min_final_depth, col = "red")
    # Zoom in
    filtered_alpha=alpha_table[alpha_table$depth < min_final_depth + 10000, ]
    if (length(rownames(filtered_alpha))>0){
      plot(filtered_alpha$depth,filtered_alpha$qTD,
           xlab="Depth (Reads per sample)",
           ylab="Asymptotic Richness (Expected number of species per sample)",
           ylim=c(0,max(filtered_alpha$qTD)+5),
           xlim=c(0,min_final_depth + 10000))
      title(main = paste0("Method: ",main_method," (Zoom in)"), line = 2.5)
      abline(v = min_final_depth, col = "red")
    }
    dev.off()
  }
  # Define objects to store the data
  filtered_trait_obj <- paste('filtered_trait',subset_id,sep='_')
  assign(filtered_trait_obj, matrix.abundance_filtered)
  filtered_trait_objs<- append(filtered_trait_objs, list(filtered_trait_obj=matrix.abundance_filtered))
}

### End of Step 2 ###

### Step 3: Collapsing tables based on higher taxonomic levels ###
print("Step 3: Collapse tables based on higher taxonomic levels")
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered_trait',subset_id,sep='_'))
    matrix = t(matrix)
    colnames(matrix)=taxonomy[[rank]][match(colnames(matrix), taxonomy$Genome)]
    collapsed_matrix=t(matrix)
    collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
    assign(paste('collapsed',rank,subset_id,sep='_'), collapsed_matrix)
    print(paste0(rank,"/", subset_id,": ",dim(get(paste('collapsed',rank,subset_id,sep='_')))))
  } 
}

### End of Step 3 ###

### Step 4: Filter samples based on the final depth ###
print("Step 4: Filter samples based on the final depth")
filtered_trait_objs <- list()
for (subset_id in subset_ids){
  #load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/keep.RData")
  #samp.keep <- get(paste("keep",subset_id,sep="_"))
  for (rank in ranks){
    matrix.abundance_filtered=get(paste('collapsed',rank,subset_id,sep="_"))
    print(paste0(rank,"/", subset_id,": ",dim(matrix.abundance_filtered)," (Before filtering out low depth samples)"))
    depth <- colSums(matrix.abundance_filtered)
    names(depth) <- colnames(matrix.abundance_filtered)
    samp.keep <- colnames(matrix.abundance_filtered)[depth > min_final_depth]
    matrix.abundance_filtered <- matrix.abundance_filtered[,samp.keep]
    # Define objects to store the data
    assign(paste("filtered.sample",rank,subset_id,sep="_"),matrix.abundance_filtered)
    print(paste0(rank,"/", subset_id,": ",dim(matrix.abundance_filtered)," (After filtering out low depth samples)"))
  }
  if (main_method=="16S") {
    matrix=get(paste('filtered.sample_ASV',subset_id,sep="_"))
    matrix = t(matrix)
    colnames(matrix)=taxonomy$Genome[match(colnames(matrix), taxonomy$ASV)]
  } else {
    matrix=get(paste('filtered.sample_Species',subset_id,sep="_"))
    matrix = t(matrix)
    colnames(matrix)=taxonomy$Genome[match(colnames(matrix), taxonomy$Species)]
  }
  collapsed_matrix=t(matrix)
  collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
  print(paste0(subset_id,": ",dim(collapsed_matrix)," (New 'filtered traits table' after filtering out low depth samples)"))
  # Define objects to store the data
  filtered_trait_obj <- paste('filtered_trait',subset_id,sep='_')
  assign(filtered_trait_obj, collapsed_matrix)
  filtered_trait_objs<- append(filtered_trait_objs, list(filtered_trait_obj=collapsed_matrix))
  
}

# Save filtered tables in a list of objects for posterior analyses
save(filtered_trait_objs, 
     file = 'filtered_traits_matrix_noprev.RData')

# Write abundance filtered tables as text files for cluster analyses
for (subset_id in subset_ids){
  filtered=get(paste("filtered_trait",subset_id,sep="_"))
  matrix=data.frame(Genome_ID=rownames(filtered))
  matrix=cbind(matrix,as.data.frame(filtered))
  write.table(matrix,file=paste0(subset_id,".csv"),sep=",",row.names=F,quote=F)
}

### End of Step 4 ###

### Step 5: Taxa filtering based on prevalence (Important for Differential Abundance Analyses) ###
print("Step 5: Taxa filtering based on prevalence")

filtered_prev_objs <- list()
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered.sample',rank,subset_id,sep='_'))
    spe.nonZero  <- rowSums(matrix>0)/ncol(matrix) # Prevalence
    spe.keep <- rownames(matrix)[spe.nonZero > min_prevalence]
    matrix.filtered <- matrix[spe.keep,]
    # Define objects to store the data
    filtered_prev_obj <- paste('filtered_prev',rank,subset_id,sep='_')
    assign(filtered_prev_obj, matrix.filtered)
    filtered_prev_objs<- append(filtered_prev_objs, list(filtered_prev_obj=matrix.filtered))
    print(paste0(rank,"/", subset_id,": ",dim(matrix.filtered)))
  } 
}

# Save results for posterior analyses

save(filtered_prev_objs,
     file="filtered_prev.RData")

### End of Step 5 ###

### Step 6: Relative abundance plots ###
print("Step 6: Phylum relative abundance plots")
pdf('phylum_stacked_barplots.pdf')
if ("Phylum" %in% ranks) {
  # Define colors
  total_phyla=c()
  for (subset_id in subset_ids){
    phylum_matrix <- get(paste("filtered_prev_Phylum",subset_id,sep="_"))
    total_phyla=unique(c(total_phyla,rownames(phylum_matrix)))
  }
  total_phyla=c(total_phyla,"Rare Taxa")
  colours = metafolio::gg_color_hue(n =  length(total_phyla), c = 100, l = 65)
  shuffled_indices <- sample(length(colours))
  colours <- colours[shuffled_indices]
  # Make the plot
  for (subset_id in subset_ids){
    # Make into relative abundance
    phylum_matrix <- apply(get(paste("filtered_prev_Phylum",subset_id,sep="_")), 2, function(i) i/sum(i))
    colSums(phylum_matrix)[1:5]
    # Define a cutoff for rare taxa
    maxabundances <- apply(phylum_matrix, 1, max)
    # Meanwhile, transpose the count table for future wrangling.
    phylum_matrix <- data.frame(t(phylum_matrix))
    # For every sample, sum up all rare taxa ( < 1% at their highest in this case)
    rare=which(maxabundances < 0.01)
    if (length(rare)>1){
      phylum_matrix$`Rare Taxa` <- rowSums(phylum_matrix[,maxabundances < 0.01], na.rm = TRUE) #Remove the individual rare taxa now that they're summed up
      phylum_matrix = phylum_matrix[,c(maxabundances > 0.01, T) ] #`T` to include the `Rare Taxa`.
    } else {
      colnames(phylum_matrix)[which(colnames(phylum_matrix)==names(maxabundances)[rare])]="Rare Taxa"
    }
    # No covariates
    colnames(phylum_matrix)=gsub("p__","",colnames(phylum_matrix))
    phylum_matrix$Sample=rownames(phylum_matrix)
    # Wrangle the data to long format for easy plotting
    barlong = phylum_matrix %>%
      pivot_longer(!Sample, names_to = c("Phylum"), values_to = "value") %>%
      mutate(Phylum = str_replace(Phylum, ".*_or_", ""))
    # Calculate mean values across covariates
    grouped_matrix <- barlong %>%
      group_by(Phylum) %>%
      summarise(mean_value_per_Phylum = mean(value, na.rm = TRUE), .groups = "drop")
    myplot <- grouped_matrix %>%
      ggplot(aes(x = factor(1), y = mean_value_per_Phylum, fill = Phylum)) +
      geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
      scale_fill_manual(values = colours) +
      labs(color = "") +
      ggtitle(paste0("Subset: ",subset_id)) +
      xlab("") +
      ylab("Mean relative abundance") +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.text=element_text(size=10),)
    print(myplot)
    # Prepare the data for ggplot by adding in metadata here
    motch = match(rownames(phylum_matrix), metadata[[sample_identifier]])
    any(is.na(motch))
    this_metadata = metadata[motch,]
    covariates=get(paste("covariates",subset_id,sep="_"))
    formula <- c("Sample")
    for (covariate in covariates){
      phylum_matrix[[covariate]]=this_metadata[[covariate]]
      phylum_matrix$Sample=rownames(phylum_matrix)
      # Wrangle the data to long format for easy plotting
      formula <- c(formula, covariate)
      barlong = phylum_matrix %>%
        pivot_longer(!formula, names_to = c("Phylum"), values_to = "value") %>%
        mutate(Phylum = str_replace(Phylum, ".*_or_", ""))
      # Calculate mean values across covariates
      grouped_matrix <- barlong %>%
        group_by(Phylum, !!sym(covariate)) %>%
        summarise(mean_value_per_Phylum = mean(value))
      # With covariates
      myplot <- grouped_matrix %>%
        ggplot(aes(x = !!sym(covariate), y = mean_value_per_Phylum, fill = Phylum)) +
        geom_bar(stat = "identity", position = "stack", col = "NA", linewidth = 0, width = 0.8) +
        scale_fill_manual(values = colours) +
        labs(color = "") +
        ggtitle(paste0("Subset: ",subset_id)) +
        xlab(paste0("Covariate: ",covariate)) +
        ylab("Mean relative abundance") +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text=element_text(size=10),)
      print(myplot)
    }
  }
} else {
  print("WARNING: As 'Phylum' was not selected as a rank no Phylum relative abundance plots will be generated")
}
dev.off()

### End of Step 6 ###
