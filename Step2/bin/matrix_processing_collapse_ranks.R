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
## Matrix processing
parser$add_argument("-input", "--input_matrix", type="character", help="File with input matrix")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the host",default="host_subject_id")
parser$add_argument("-mfd", "--min_final_depth", type="integer", help="Minimum number of reads accepted after filtering noise", default=30000)
parser$add_argument("-mp", "--min_prevalence", type="double", help="Minimum prevalence accepted", default=0.5)
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Matrix processing
input_matrix<-args$input_matrix
input_meta<-args$input_meta
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
min_final_depth<-args$min_final_depth
min_prevalence<-args$min_prevalence
method=args$method

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
print("Step 1: Read input files and arguments")

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

# Read input matrices
if (input_matrix != "") {
  print("Reading filtered matrices per subset")
  load(input_matrix)
  for (i in 1:length(filtered_trait_objs)){
    filtered_trait=filtered_trait_objs[[i]]
    assign(paste('filtered_trait',subset_ids[i],sep='_'), filtered_trait)
  }
} else {
  stop("No matrices provided")
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

### End of Step 1 ###

### Step 2: Collapsing tables based on higher taxonomic levels ###
print("Step 2: Collapse tables based on higher taxonomic levels")
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

### Step 3: Filter samples based on the final depth ###
print("Step 3: Filter samples based on the final depth")
for (subset_id in subset_ids){
  for (rank in ranks){
    matrix.abundance_filtered=get(paste('collapsed',rank,subset_id,sep="_"))
    print(paste0(rank,"/", subset_id,": ",dim(matrix.abundance_filtered)," (Before filtering out low depth samples)"))
    depth <- colSums(matrix.abundance_filtered)
    names(depth) <- colnames(matrix.abundance_filtered)
    samp.keep <- colnames(matrix.abundance_filtered)[depth >= min_final_depth]
    matrix.abundance_filtered <- matrix.abundance_filtered[,samp.keep]
    # Define objects to store the data
    assign(paste("filtered.sample",rank,subset_id,sep="_"),matrix.abundance_filtered)
    print(paste0(rank,"/", subset_id,": ",dim(matrix.abundance_filtered)," (After filtering out low depth samples)"))
  }
}

### End of Step 3 ###

### Step 4: Taxa filtering based on prevalence (Important for Differential Abundance Analyses) ###
print("Step 4: Taxa filtering based on prevalence")

filtered_prev_objs <- list()
for (rank in ranks){
  for (subset_id in subset_ids) {
    matrix = get(paste('filtered.sample',rank,subset_id,sep='_'))
    # Calculate the prevalence of each trait
    spe.nonZero  <- rowSums(matrix>0)/ncol(matrix) # Prevalence
    # Get the list of traits that are present in more than 50% samples
    spe.keep <- rownames(matrix)[spe.nonZero > min_prevalence]
    length(spe.keep) # Check how many traits are in the list.
    #Remove low prevalent traits
    matrix.filtered <- matrix[spe.keep,]
    if (method=="16S"){
      remove=c("s__","g__","f__","c__","o__","p__")
      matrix.filtered=matrix.filtered[! rownames(matrix.filtered) %in% remove,]
    }
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

### End of Step 4 ###

### Step 5: Relative abundance plots ###
print("Step 5: Phylum relative abundance plots")
pdf('phylum_stacked_barplots.pdf',bg='white')
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
              axis.title = element_text(size = 15),
              legend.text=element_text(size=10),)
      print(myplot)
    }
  }
} else {
  print("WARNING: As 'Phylum' was not selected as a rank no Phylum relative abundance plots will be generated")
}
dev.off()

### End of Step 4 ###