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

# 3. Clean ids
clean_ids<-function(oldids) {
  raw_ids<-str_split(oldids, pattern="_", simplify = TRUE)
  ids<-raw_ids[,1]
  if(!all(nchar(ids)>=10)){
    for(i in 1:length(ids)){
      if(!nchar(ids[i])>=10){ids[i] <- paste0("000", ids[i])}
    }
  }
  return(ids)
}

# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Matrix processing
parser$add_argument("-mod_ana", "--module_analysis", type="character", help="Inclusion of module: Analysis",default="Pooled")
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-input_harm", "--input_for_harmonization", type="character", help="File with input matrix")
parser$add_argument("-exclude", "--exclude_ids", type="character", help="File with list of ids to be removed (optional)")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-lib", "--input_libsize", type="character", help="File with library depth info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-mets", "--methods", type="character", help="Define the methods to be compared", default="Shallow,16S")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the sample",default="host_subject_id")
parser$add_argument("-mid", "--min_initial_depth", type="integer", help="Minimum number of raw reads accepted for the main dataset", default=0)
parser$add_argument("-ma", "--min_abundance", type="double", help="Minimum median relative abundance accepted", default=0.0001)
parser$add_argument("-mfd", "--min_final_depth", type="integer", help="Minimum number of reads accepted after filtering noise for the main dataset", default=30000)
parser$add_argument("-mid_harm", "--min_initial_depth_for_harmonization", type="integer", help="Minimum number of raw reads accepted for the secondary dataset", default=0)
parser$add_argument("-mfd_harm", "--min_final_depth_for_harmonization", type="integer", help="Minimum number of reads accepted after filtering noise for the secondary dataset", default=6000)
parser$add_argument("-mp", "--min_prevalence", type="double", help="Minimum prevalence accepted", default=0.5)
parser$add_argument("-mat", "--sample_material", type="character", help="Define the type of material used", default="Cecal")
parser$add_argument("-hyb", "--hybrid", type="character", help="Define if hybrid matrices will be built with both datasets", default="NO")
parser$add_argument("-prop", "--proportion", type="double", help="Porportion of the main dataset that the hybrid matrix should have", default=0.0)

## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-taxo_harm", "--input_taxonomy_for_harmonization", type="character", help="File with taxonomy info")
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

## Matrix processing
module_analysis<-args$module_analysis
input_matrix<-args$input_matrix
input_matrix_harm<-args$input_for_harmonization
exclude_file<-args$exclude_ids
compare_ids<-args$compare_ids
input_meta<-args$input_meta
input_libsize<-args$input_libsize
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
min_initial_depth<-args$min_initial_depth
min_abundance<-args$min_abundance
min_final_depth<-args$min_final_depth
min_initial_depth_harm<-args$min_initial_depth_for_harmonization
min_final_depth_harm<-args$min_final_depth_for_harmonization
min_prevalence<-args$min_prevalence
methods<-args$methods
main_method<-args$method
Material=args$sample_material
if (Material == "") {
  stop("Please inform the sample type under investigation if you want to run any Host Genetic Effects Analyses. It should correspond to the tissue type where the metadata was collected from (i.e. Cecal, Fecal, etc.)", "\n")
}
hybrid=args$hybrid
proportion=args$proportion

## Taxonomic
input_taxonomy<-args$input_taxonomy
input_taxonomy_harm<-args$input_taxonomy_for_harmonization
ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus or Species" ,"\n")
}

### Step 1 ####
print("Step 1: Read input files and arguments and do the custom removal of bad samples")

# Read methods
print("Reading methods")
methods<-na.omit(c((sapply(strsplit(methods, ","),"[",1)),
                   (sapply(strsplit(methods, ","), "[",2))))
if (methods[1]!=main_method){
  stop("ERROR: The 1st method pointed on the 'methods' field should be the main method to be used for profiling.")
}

# Reading input files
for (method in methods){
  if (method == methods[1]){
    if (input_matrix != "") {
      print("Reading main counts table")
      print(method)
      if (tail(strsplit(input_matrix, "\\.")[[1]], 1)=="biom"){
        counts_table=as.matrix(biom_data(read_biom(input_matrix)))
        colnames(counts_table) <- str_sub(colnames(counts_table), - 10, - 1)
      } else {
        counts_table<-read_tab_file(input_matrix, TRUE)
        colnames(counts_table)<-tools::file_path_sans_ext(colnames(counts_table))
        row.names(counts_table)<-counts_table$Genome_id
        counts_table$Genome_id<-NULL
      }
      assign(paste("counts_table",method,sep="_"),counts_table)
    } else {
      stop("No main counts table provided")
    }
  } else {
    if (input_matrix_harm != ""){
      print("Reading secondary counts table")
      if (tail(strsplit(input_matrix_harm, "\\.")[[1]], 1)=="biom"){
        counts_table=as.matrix(biom_data(read_biom(input_matrix_harm)))
        colnames(counts_table) <- str_sub(colnames(counts_table), - 10, - 1)
      } else {
        counts_table<-read_tab_file(input_matrix_harm, TRUE)
        colnames(counts_table)<-tools::file_path_sans_ext(colnames(counts_table))
        row.names(counts_table)<-counts_table$Genome_id
        counts_table$Genome_id<-NULL
      }
      assign(paste("counts_table",method,sep="_"),counts_table)
    } else {
      stop("No secondary counts table provided")
    }
  }
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
for (method in methods){
  if (method=="16S"){
    depths<-data.frame(Sample=metadata[[sample_identifier]],Raw_Reads=metadata$Depth_16S)
  } else {
    if (input_libsize != "") {
      depths<-read_tab_file(input_libsize, TRUE)
    } else {
      stop("No library depth file provided")
    }
  }
  assign(paste("depths",method,sep="_"),depths)
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

if (module_analysis=="Pooled"){
  studies=subset_ids[-which(subset_ids=="ALL")]
  subset_ids="ALL"
}

for (method in methods){
  if (method == methods[1]) {
    print("Reading taxonomy info for the main dataset")
    if (input_taxonomy != "") {
      taxonomy<-read_tab_file(input_taxonomy, FALSE)
      if (grepl("GTDB", input_taxonomy)) {
        colnames(taxonomy)=c("Genome","Taxonomy")
        taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
        taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
      } else {
        if (grepl("Greengenes2", input_taxonomy)){
          colnames(taxonomy) <- c("Genome","Taxonomy","Confidence")
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
      assign(paste("taxonomy",method,sep="_"),taxonomy)
    } else {
      stop("No taxonomy file provided for the main dataset")
    }
  } else {
    print("Reading taxonomy info for the secondary dataset")
    if (input_taxonomy_harm != "") {
      taxonomy <- read.delim(input_taxonomy_harm,header=TRUE)
      if (grepl("GTDB", input_taxonomy_harm)) {
        colnames(taxonomy)=c("Genome","Taxonomy")
        taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
        taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
      } else {
        if (grepl("Greengenes2", input_taxonomy_harm)){
          colnames(taxonomy) <- c("Genome","Taxonomy","Confidence")
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
      assign(paste("taxonomy",method,sep="_"),taxonomy)
    } else {
      stop("No taxonomy file provided for the secondary dataset")
    }
  }
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

# Perform custom clean of samples
print("Performing custom cleaning")

for (method in methods){
  counts_table=get(paste("counts_table",method,sep="_"))
  colnames(counts_table)<-clean_ids(colnames(counts_table))
  print(method)
  assign(paste("counts_table",method,sep="_"),counts_table)
  depths=get(paste("depths",method,sep="_"))
  depths$Sample<-clean_ids(depths$Sample)
  assign(paste("depths",method,sep="_"),depths)
}
if (exclude_file != "") {
  exclude_ids$Files<-clean_ids(exclude_ids$Files)
}

# Remove controls to be used only on Beta Diversity Analyses
print("Removing controls to be used only on Beta Diversity Analyses")
matrix.sample_controls_objs=c()
metadata_controls_objs=c()
for (method in methods){
  counts_table=get(paste("counts_table",method,sep="_"))
  motch = match(colnames(counts_table),metadata[[sample_identifier]]) 
  #if (any(is.na(motch))) {
  #  missing_samples <- colnames(counts_table)[!(colnames(counts_table) %in% metadata[[sample_identifier]])]
  #  stop(paste("These samples were not found in metadata:", paste(missing_samples, collapse = ", ")))
  #}
  this_metadata = metadata[motch,]
  metadata_controls=this_metadata[grepl("Positive Control",this_metadata$Material) | grepl("Negative Control",this_metadata$Material),]
  this_metadata=this_metadata[!this_metadata[[sample_identifier]] %in% metadata_controls[[sample_identifier]],]
  matrix.sample_controls=counts_table[,match(metadata_controls[[sample_identifier]],colnames(counts_table))]
  counts_table=counts_table[,!colnames(counts_table) %in% colnames(matrix.sample_controls)]
  # Define objects to store the data
  matrix.sample_controls_obj <- paste('matrix.sample_controls',method,sep='_')
  assign(matrix.sample_controls_obj, matrix.sample_controls)
  matrix.sample_controls_objs<- append(matrix.sample_controls_objs, list(fmatrix.sample_controls_obj=matrix.sample_controls))
  metadata_controls_obj <- paste('metadata_controls',method,sep='_')
  assign(metadata_controls_obj, metadata_controls)
  metadata_controls_objs<- append(metadata_controls_objs, list(fmetadata_controls_obj=metadata_controls))
  assign(paste("counts_table",method,sep="_"),counts_table)
  print(method)
}

# Save controls
save(matrix.sample_controls_objs,
     metadata_controls_objs,
     file = 'controls.RData')

# Filter by Material (This is important to keep only one sample per host for host genetic analyses)
print("Filtering by Material (This is important to keep only one sample per host for host genetic analyses)")
metadata = metadata[grepl(Material,metadata$Material),]
for (method in methods){
  counts_table=get(paste("counts_table",method,sep="_"))
  counts_table = counts_table[,match(metadata[[sample_identifier]],colnames(counts_table))]
  assign(paste("counts_table",method,sep="_"),counts_table)
  print(paste0(method," : ",dim(counts_table)," (After filtering by material)"))
}


# Identify low covered samples
for (method in methods){
  depths=get(paste("depths",method,sep="_"))
  if(method==methods[1]){
    low_depth_ids<-depths[depths$"Raw_Reads"<min_initial_depth, ]$Sample
  } else{
    low_depth_ids<-depths[depths$"Raw_Reads"<min_initial_depth_harm, ]$Sample 
  }
  assign(paste("low_depth_ids",method,sep="_"),low_depth_ids)
}


# Remove bad samples from a custom list and/or with low depth
if (exclude_file != "") {
  exclude_ids<-c(exclude_ids[,1], get(paste("low_depth_ids",methods[1],sep="_")),get(paste("low_depth_ids",methods[2],sep="_")))
} else {
  exclude_ids<-c(get(paste("low_depth_ids",methods[1],sep="_")),get(paste("low_depth_ids",methods[2],sep="_")))
}
for (method in methods){
  counts_table=get(paste("counts_table",method,sep="_"))
  matrix.sample_filtered <- counts_table[,!colnames(counts_table) %in% exclude_ids]
  print(paste0(method," : ",dim(matrix.sample_filtered)," (After removing bad samples and samples with initial low depths)"))
  assign(paste("matrix.sample_filtered",method,sep="_"),matrix.sample_filtered)
}

# 1st Sample harmonization
print("1st Sample harmonization")
matrix.sample_filtered=get(paste("matrix.sample_filtered",methods[1],sep="_"))
matrix.sample_filtered_harm=get(paste("matrix.sample_filtered",methods[2],sep="_"))
print(paste0("Main dataset: ",dim(matrix.sample_filtered)," (Before sample harmonization)"))
matrix.sample_filtered=matrix.sample_filtered[,colnames(matrix.sample_filtered) %in% colnames(matrix.sample_filtered_harm)]
print(paste0("Main dataset: ",dim(matrix.sample_filtered)," (After sample harmonization)"))
print(paste0("Secondary dataset: ",dim(matrix.sample_filtered_harm)," (Before sample harmonization)"))
matrix.sample_filtered_harm=matrix.sample_filtered_harm[,colnames(matrix.sample_filtered_harm) %in% colnames(matrix.sample_filtered)]
print(paste0("Secondary dataset: ",dim(matrix.sample_filtered_harm)," (After sample harmonization)"))
assign(paste("matrix.sample_filtered",methods[1],sep="_"),matrix.sample_filtered)
assign(paste("matrix.sample_filtered",methods[2],sep="_"),matrix.sample_filtered_harm)

# Subset metadata info

matrix=get(paste("matrix.sample_filtered",methods[1],sep="_"))
motch = match(colnames(matrix),metadata[[sample_identifier]])
if (any(is.na(motch))) {
  print("WARNING SOME SAMPLE IS MISSING IN THE METADATA")
}
metadata = metadata[motch,]

### End of Step 1 ###

### Step 2: Taxa filtering based on relative abundance (Important to remove noise) ###
print("Step 2: Filter low abundance taxa per study")

filtered_trait_objs <- list()
for (method in methods) {
  matrix_full=get(paste("matrix.sample_filtered",method,sep="_"))
  for (subset_id in subset_ids) {
    if(subset_id != "ALL") {
      matrix = matrix_full[,which(metadata[,"Study"] == subset_id)]  
    } else {
      samp.keep=metadata[[sample_identifier]][metadata$Study %in% studies]
      matrix = matrix_full[,samp.keep]
      }
    # Identify low abundant traits
    # Re-scale the table in relative values, make sure the sum of each row is 1
    print(paste0(method," / ",subset_id,": ",dim(matrix)," (Before filtering out low abundance taxa)"))
    matrix.relative <- sweep(matrix,2,colSums(matrix),"/")
    # Calculate the mean relative abundance of each trait
    spe.meanAbun <- apply(matrix.relative, 1, mean, na.rm=TRUE)
    # Filter based on the mean relative abundance
    spe.keep <- rownames(matrix.relative)[spe.meanAbun > min_abundance]
    matrix.abundance_filtered <- as.matrix(matrix[spe.keep,])
    print(paste0(method," / ",subset_id,": ",dim(matrix.abundance_filtered)," (After filtering out low abundance taxa)"))
    if (subset_id == "ALL" & method !="16S") {
      # Make the depth funnel plot
      print("Making the depth funnel plot")
      depth <- colSums(matrix.abundance_filtered)
      depths=get(paste("depths",method,sep="_"))
      depths=depths[depths$Sample %in% colnames(matrix.abundance_filtered),]
      depths=depths[match(colnames(matrix.abundance_filtered), depths$Sample),]
      depths$Mapped_Reads_After_Abundance_Filtering=depth
      depths$Clear_Mapping_Ratio=100*depths$Mapped_Reads_After_Abundance_Filtering/depths$Non_Host_Reads
      save(depths,file = 'depth_funnel_focal_method.RData')
      depths$Raw_Reads=depths$Raw_Reads/1000
      depths$Reads_After_Trimming=depths$Reads_After_Trimming/1000
      depths$Non_Host_Reads=depths$Non_Host_Reads/1000
      depths$Mapped_Reads=depths$Mapped_Reads/1000
      depths$Mapped_Reads_After_Abundance_Filtering=depths$Mapped_Reads_After_Abundance_Filtering/1000
      pdf('depth_funnel_focal_method.pdf')
      depths_long <- tidyr::pivot_longer(depths, cols = c(Raw_Reads, Reads_After_Trimming, Non_Host_Reads, Mapped_Reads, Mapped_Reads_After_Abundance_Filtering),
                                         names_to = "Funnel Step", values_to = "Depth (Thousands of Reads)")
      order1 <- c("Raw_Reads", "Reads_After_Trimming", "Non_Host_Reads", "Mapped_Reads", "Mapped_Reads_After_Abundance_Filtering")
      order2 <- c("Initial_Rate","Surviving_Rate_Trimming", "Surviving_Rate_Alignment", "Raw_Mapping_Ratio", "Clear_Mapping_Ratio")
      labels <- c("Raw Reads", "Reads After Trimming", "Non-Host Reads", "Mapped Reads", "Mapped Reads After Abundance Filtering")
      median_values <- depths %>%
        summarise(
          Initial_Rate=100,
          Surviving_Rate_Trimming = mean(Surviving_Rate_Trimming),
          Surviving_Rate_Alignment = mean(Surviving_Rate_Alignment),
          Raw_Mapping_Ratio = mean(Raw_Mapping_Ratio),
          Clear_Mapping_Ratio = mean(Clear_Mapping_Ratio)
        ) %>%
        pivot_longer(cols = everything(), names_to = "Surviving Rates", values_to = "Median Surviving Rate")
      median_values$`Median Surviving Rate`=round(median_values$`Median Surviving Rate`,digits = 2)
      scale=median(depths$Raw_Reads)/median(median_values$`Median Surviving Rate`)
      myplot<-ggplot(depths_long, aes(x = factor(`Funnel Step`, levels = order1), y = `Depth (Thousands of Reads)`)) +
        geom_boxplot() +
        geom_label(data = median_values,  
                   aes(x = as.numeric(factor(`Surviving Rates`, levels = order2)),
                       y = `Median Surviving Rate` *scale, 
                       label = paste0(`Median Surviving Rate`, "%"), fill = factor(`Surviving Rates`, levels = order2)),
                   vjust = -0.5, size = 3.5) +
        scale_x_discrete(labels = labels) +
        scale_fill_manual(values = c("Initial_Rate" = "lightcoral","Surviving_Rate_Trimming" = "beige", "Surviving_Rate_Alignment" = "lightgreen", "Raw_Mapping_Ratio" = "lightblue", "Clear_Mapping_Ratio" = "lightpink"),
                          labels = c("Starting Point", "Surviving: After Trimming", "Survinging: Non-Host", "Raw Mapping Ratio", "Clean Mapping Ratio"),
                          name = "Mean Filtering Rates") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        guides(fill = guide_legend(override.aes = list(color = NA)), 
               color = "none", 
               shape = "none") +
        labs(x = NULL)
      print(myplot)
      dev.off()
    }
    # Define objects to store the data
    filtered_trait_obj <- paste('filtered_trait',subset_id,method,sep='_')
    assign(filtered_trait_obj, matrix.abundance_filtered)
    filtered_trait_objs<- append(filtered_trait_objs, list(filtered_trait_obj=matrix.abundance_filtered))
  }
}

# Salve abundance filtered tables
save(filtered_trait_objs,
     file='filtered_traits_matrix_noprev.RData')


# Make asymptotic richness plots
print("Making asymptotic richness plots")
pdf('asymptotic_richness_vs_depth.pdf')
for (method in methods) {
  matrix=get(paste("filtered_trait_ALL",method,sep="_"))
  matrix=round(matrix, digits = 0)
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
  if (method==methods[1]){
    save(alpha_table,
         file="alpha_table_main.RData")
    this_min_final_depth=min_final_depth
  } else {
    save(alpha_table,
         file="alpha_table_secondary.RData")
    this_min_final_depth=min_final_depth_harm
  }
  if (method=="16S"){
    text="ASVs"
  } else {text="OGUs"}
  plot(alpha_table$depth,alpha_table$qTD,xlab="Depth (Pairs of reads per sample)",ylab=paste0("Asymptotic Richness (Expected number of ",text," per sample)"))
  title(main = paste0("Method: ",method), line = 2.5)
  abline(v = this_min_final_depth, col = "red")
  # Zoom in
  filtered_alpha=alpha_table[alpha_table$depth < this_min_final_depth + 10000, ]
  plot(filtered_alpha$depth,filtered_alpha$qTD,xlab="Depth (Pairs of reads per sample)",ylab=paste0("Asymptotic Richness (Expected number of ",text," per sample)"))
  title(main = paste0("Method: ",method," (Zoom in)"), line = 2.5)
  abline(v = this_min_final_depth, col = "red")
  print(method)
}
dev.off()


### End of Step 2 ###

### Step 3: Collapse into taxonomic ranks ###
print("Step 3: Collapse into taxonomic ranks")

for (method in methods){
  for (rank in ranks){
    for (subset_id in subset_ids) {
      matrix = get(paste('filtered_trait',subset_id,method,sep='_'))
      taxonomy=get(paste('taxonomy',method,sep='_'))
      matrix = t(matrix)
      colnames(matrix)=taxonomy[[rank]][match(colnames(matrix), taxonomy$Genome)]
      collapsed_matrix=t(matrix)
      collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
      assign(paste('collapsed',rank,subset_id,method,sep='_'), collapsed_matrix)
      print(paste0(rank,"/", subset_id,"/", method,": ",dim(get(paste('collapsed',rank,subset_id,method,sep='_')))))
    } 
  }
}

### End of Step 3 ###

### Step 4: Filter samples based on the final depth ###
print("Step 4: Filter samples based on the final depth")
for (method in methods){
  for (subset_id in subset_ids){
    for (rank in ranks){
      matrix.abundance_filtered=get(paste('collapsed',rank,subset_id,method,sep="_"))
      print(paste0(rank,"/", subset_id,"/",method,": ",dim(matrix.abundance_filtered)," (Before filtering out low depth samples)"))
      if (method==methods[1]){
        this_min_final_depth=min_final_depth
      } else {
        this_min_final_depth=min_final_depth_harm
      }
      depth <- colSums(matrix.abundance_filtered)
      names(depth) <- colnames(matrix.abundance_filtered)
      samp.keep <- colnames(matrix.abundance_filtered)[depth > this_min_final_depth]
      matrix.abundance_filtered <- matrix.abundance_filtered[,samp.keep]
      # Define objects to store the data
      assign(paste("filtered.sample",rank,subset_id,method,sep="_"),matrix.abundance_filtered)
      print(paste0(rank,"/", subset_id,"/",method,": ",dim(matrix.abundance_filtered)," (After filtering out low depth samples)"))
    }
  }
}

### End of Step 4 ###

### Step 5: 2nd Sample harmonization ###
print("Step 5: 2nd Sample harmonization")
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste("filtered.sample",rank,subset_id,methods[1],sep='_'))
    matrix_harm = get(paste("filtered.sample",rank,subset_id,methods[2],sep='_'))
    print(paste0(rank,"/", subset_id,"/","Main matrix: ",dim(matrix)," (Before sample harmonization)"))
    print(paste0(rank,"/", subset_id,"/","Secondary matrix: ",dim(matrix_harm)," (Before sample harmonization)"))
    matrix=matrix[,colnames(matrix) %in% colnames(matrix_harm)]
    matrix_harm=matrix_harm[,colnames(matrix_harm) %in% colnames(matrix)]
    print(paste0(rank,"/", subset_id,"/","Main matrix: ",dim(matrix)," (After sample harmonization)"))
    print(paste0(rank,"/", subset_id,"/","Secondary matrix: ",dim(matrix_harm)," (After sample harmonization)"))
    assign(paste('harmonized.sample',rank,subset_id,methods[1],sep='_'), matrix)
    assign(paste('harmonized.sample',rank,subset_id,methods[2],sep='_'), matrix_harm)
  } 
}

### End of Step 5 ###

# Step 6: Filter by prevalence
print("Step 6: Taxa filtering based on prevalence")

filtered_prev_objs <- list()
for (method in methods){
  for (rank in ranks){
    for (subset_id in subset_ids) {
      matrix = get(paste('harmonized.sample',rank,subset_id,method,sep='_'))
      print(paste0(rank,"/", subset_id,"/",method,": ",dim(matrix)," (Before filtering out low prevalence taxa)"))
      # Calculate the prevalence of each trait
      spe.nonZero  <- rowSums(matrix>0)/ncol(matrix) # Prevalence
      # Get the list of traits that are present in more than 50% samples
      spe.keep <- rownames(matrix)[spe.nonZero > min_prevalence]
      length(spe.keep) # Check how many traits are in the list.
      #Remove low prevalent traits
      matrix.filtered <- matrix[spe.keep,]
      print(paste0(rank,"/", subset_id,"/",method,": ",dim(matrix.filtered)," (After filtering out low prevalence taxa)"))
      # Define objects to store the data
      filtered_prev_obj <- paste('filtered_prev',rank,subset_id,method,sep='_')
      assign(filtered_prev_obj, matrix.filtered)
      filtered_prev_objs<- append(filtered_prev_objs, list(filtered_prev_obj=matrix.filtered))
      print(paste0(method,"/",rank,"/", subset_id,": ",dim(matrix.filtered)))
    } 
  }
}

save(filtered_prev_objs,
     file="harmonized_samples.RData")


### End of Step 6 ###


# Step 7: Taxa harmonization
print("Step 7: Taxa harmonization")
harmonized_taxa_objs <- list()
harmonized_taxa_harm_objs <- list()
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered_prev',rank,subset_id,methods[1],sep='_'))
    matrix_harm = get(paste('filtered_prev',rank,subset_id,methods[2],sep='_'))
    print(paste0(rank,"/", subset_id,"/","Main matrix: ",dim(matrix)," (Before taxa harmonization)"))
    print(paste0(rank,"/", subset_id,"/","Secondary matrix: ",dim(matrix_harm)," (Before taxa harmonization)"))
    matrix_harm=matrix_harm[rownames(matrix_harm) %in% rownames(matrix),]
    matrix=matrix[rownames(matrix) %in% rownames(matrix_harm),]
    print(paste0(rank,"/", subset_id,"/","Main matrix: ",dim(matrix)," (After taxa harmonization)"))
    print(paste0(rank,"/", subset_id,"/","Secondary matrix: ",dim(matrix_harm)," (After taxa harmonization)"))
    # Define objects to store the data
    harmonized_taxa_obj <- paste('harmonized_taxa',rank,subset_id,methods[1],sep='_')
    assign(harmonized_taxa_obj, matrix)
    harmonized_taxa_harm_obj <- paste('harmonized_taxa',rank,subset_id,methods[2],sep='_')
    assign(harmonized_taxa_harm_obj, matrix_harm)
    harmonized_taxa_objs<- append(harmonized_taxa_objs, list(harmonized_taxa_obj=matrix))
    harmonized_taxa_harm_objs<- append(harmonized_taxa_harm_objs, list(harmonized_taxa_harm_obj=matrix_harm))
  } 
}

if (hybrid=="YES"){
  hybrid_matrix_objs <- list()
  hybrid_matrix_harm_objs <- list()
  for (rank in ranks){
    for (subset_id in subset_ids) {
      matrix = get(paste('harmonized_taxa',rank,subset_id,methods[1],sep='_'))
      matrix_harm = get(paste('harmonized_taxa',rank,subset_id,methods[2],sep='_'))
      hybrid_matrix <- (proportion*matrix) + ((1-proportion)*matrix_harm)
      hybrid_matrix_harm <- (proportion*matrix_harm) + ((1-proportion)*matrix)
      # Define objects to store the data
      hybrid_matrix_obj <- paste('hybrid_matrix',rank,subset_id,methods[1],sep='_')
      assign(hybrid_matrix_obj, hybrid_matrix)
      hybrid_matrix_harm_obj <- paste('hybrid_matrix',rank,subset_id,methods[2],sep='_')
      assign(hybrid_matrix_harm_obj, hybrid_matrix_harm)
      hybrid_matrix_objs<- append(hybrid_matrix_objs, list(hybrid_matrix_obj=hybrid_matrix))
      hybrid_matrix_harm_objs<- append(hybrid_matrix_harm_objs, list(hybrid_matrix_harm_obj=hybrid_matrix_harm))
    }
  }
}

print("Saving harmonized samples and taxa")
filtered_prev_objs <- list()
for (method in methods){
  for (rank in ranks) {
    for (subset_id in subset_ids){
      if (hybrid=="YES"){
        matrix = get(paste('hybrid_matrix',rank,subset_id,method,sep='_'))
      } else {
        matrix = get(paste('harmonized_taxa',rank,subset_id,method,sep='_'))
      }
      # Define objects to store the data
      filtered_prev_obj <- paste('filtered_prev',rank,subset_id,methods[1],sep='_')
      assign(filtered_prev_obj, matrix)
      filtered_prev_objs<- append(filtered_prev_objs, list(filtered_prev_obj=matrix))
    }
  }
}
save(filtered_prev_objs,
     file="harmonized_samples_and_taxa.RData")

### End of Step 7 ###

#Step 8: Relative abundance plots
print("Step 8: Phylum relative abundance plots")

if ("Phylum" %in% ranks) {
  # Define colors
  harmonizations=c("Only Samples Harmonized","Samples and Taxa Harmonized")
  pdf('phylum_stacked_barplots.pdf',bg='white')
  for (harmonization in harmonizations){
    total_phyla=c()
    for (method in methods){
      for (subset_id in subset_ids){
        if (harmonization==harmonizations[1]){
          phylum_matrix = get(paste('filtered_prev_Phylum',subset_id,method,sep='_'))
        } else {
          phylum_matrix = get(paste('harmonized_taxa_Phylum',subset_id,method,sep='_'))
        }
        total_phyla=unique(c(total_phyla,rownames(phylum_matrix)))
      }
    }
      total_phyla=c(total_phyla,"Rare Taxa")
      colours = metafolio::gg_color_hue(n =  length(total_phyla), c = 100, l = 65)
      shuffled_indices <- sample(length(colours))
      colours <- colours[shuffled_indices]
      # Make the plot
      for (method in methods) {
        for (subset_id in subset_ids){
          # Make into relative abundance
          if (harmonization==harmonizations[1]){
            phylum_matrix <- apply(get(paste('filtered_prev_Phylum',subset_id,method,sep='_')), 2, function(i) i/sum(i))
          } else {
            phylum_matrix <- apply(get(paste('harmonized_taxa_Phylum',subset_id,method,sep='_')), 2, function(i) i/sum(i))
          }
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
            myplot <- grouped_matrix %>%
              ggplot(aes(x = !!sym(covariate), y = mean_value_per_Phylum, fill = Phylum)) +
              geom_bar(stat = "identity", position = "stack", col = "white", linewidth = 2, width = 1) +
              scale_fill_manual(values = colours,
                                labels = total_phyla,
                                limits = total_phyla) +
              labs(color = "") +
              ggtitle(paste0(harmonization,"/",method,"/",subset_id)) +
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
                    axis.title = element_text(size = 15),
                    legend.text=element_text(size=10),)
            print(myplot)
          }
        }
      }
  }
  dev.off()
} else {
  print("WARNING: As 'Phylum' was not selected as a rank no Phylum relative abundance plots will be generated")
}