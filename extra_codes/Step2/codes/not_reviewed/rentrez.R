# Load required libraries
library(rentrez, lib.loc = "/nfs/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(XML, lib.loc = "/nfs/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(xfun, lib.loc = "/nfs/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(dplyr)
library(ggplot2)
library(tidyverse)

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

#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_Maternal/Cluster_Analyses/residuals_qned_counts_new.RData")
#species_gtdb=rownames(residuals_qned_counts_objs[[18]])
#species_gtdb=taxonomy$Genome[match(species_gtdb,taxonomy$Species)]
load("/users/abaud/fmorillo/paper_figures/species_gtdb.RData")
species_gtdb=c("Genome",species_gtdb)
load("/users/abaud/fmorillo/paper_figures/species_wol.RData")
species_wol=c("Genome",species_wol)


write.table(
  species_gtdb,            # Replace with your data frame or matrix
  "/users/abaud/fmorillo/paper_figures/species_gtdb.csv",         # Output file name
  sep = ",",            # Comma separator for CSV
  col.names = FALSE,    # Exclude column headers
  row.names = FALSE,    # Exclude row names
  quote = FALSE         # Disable quoting of values
)

write.table(
  species_wol,            # Replace with your data frame or matrix
  "/users/abaud/fmorillo/paper_figures/species_wol.csv",         # Output file name
  sep = ",",            # Comma separator for CSV
  col.names = FALSE,    # Exclude column headers
  row.names = FALSE,    # Exclude row names
  quote = FALSE         # Disable quoting of values
)


################################################################################################

pdf("~/paper_figures/gtdb_vs_wol_isolation_sources.pdf")

wol_classified <- read.csv("/nfs/users/abaud/fmorillo/paper_figures/wol_classified.csv")
wol_classified$Database="WoL"
wol_undefined=wol_classified$Type[wol_classified$Type=="Undefined"]

gtdb_classified <- read.csv("/nfs/users/abaud/fmorillo/paper_figures/gtdb_classified.csv")
gtdb_classified$Database="GTDB"
gtdb_undefined=gtdb_classified$Type[gtdb_classified$Type=="Undefined"]

concat=rbind(wol_classified,gtdb_classified)
stats_undefined=count(concat[concat$Type=="Undefined",],Database)
print(stats_undefined)

save(gtdb_undefined,
     wol_undefined,
     stats_undefined,
     file="/nfs/users/abaud/fmorillo/paper_figures/undefined_wol_gtdb.RData"
     )

concat=concat[concat$Type!="Undefined",]
unique(concat$Type)

concat$Host=NA
concat$Host[grepl("Rat ",concat$Type)]="Rat"
concat$Host[grepl("Mouse ",concat$Type)]="Mouse"
concat$Host[grepl("Other rodents ",concat$Type)]="Other rodents"
concat$Host[grepl("Human ",concat$Type)]="Human"
concat$Host[grepl("Other mammals ",concat$Type)]="Other mammals"
concat$Host[grepl("Other animals ",concat$Type)]="Other animals"
concat$Host[grepl("Undefined host ",concat$Type)]="Undefined host"
concat$Host[is.na(concat$Host)]="Other sources"
unique(concat$Host)

concat$Tissues=NA
concat$Tissues[grepl("gut",concat$Type) & !grepl("non-gut ",concat$Type)]="Gut"
concat$Tissues[grepl("non-gut",concat$Type)]="Non-gut"
concat$Tissues[grepl("Undefined tissue",concat$Type)]="Undefined tissue"
concat$Tissues[is.na(concat$Tissues)]="Non-host associated"
unique(concat$Tissues)

hosts=c(
  "Rat",
  "Mouse",
  "Other rodents",
  "Human",
  "Other mammals",
  "Other animals",
  "Undefined host",
  "Other sources"
)

statistics_host <- concat %>%
  count(Database, Host) %>%                     # your original step
  complete(Database, Host, fill = list(n = 0))  # <-- add the missing rows

max_value=max(statistics_host$n)
myplot <- ggplot(statistics_host, aes(x = factor(Host, levels = hosts), y = n, fill = Database)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, max_value+10)) +
  labs(x = "Host", y = "Number of mapped species") +
scale_fill_manual(values = c("GTDB" = "darkgreen", "WoL" = "lightgreen"),
                  labels = c("GTDB", "WoL"),
                  name = "Database") +
scale_x_discrete(guide = guide_axis(angle = 60)) +
geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5)
print(myplot)

tissues=c(
  "Gut",
  "Non-gut",
  "Undefined tissue",
  "Non-host associated"
)

statistics_tissues <- concat %>%
  count(Database, Tissues) %>%                     # your original step
  complete(Database, Tissues, fill = list(n = 0))  # <-- add the missing rows

max_value=max(statistics_tissues$n)
myplot <- ggplot(statistics_tissues, aes(x = factor(Tissues, levels = tissues), y = n, fill = Database)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, max_value+10)) +
  labs(x = "Tissue", y = "Number of mapped species") +
  scale_fill_manual(values = c("GTDB" = "darkblue", "WoL" = "lightblue"),
                    labels = c("GTDB", "WoL"),
                    name = "Database") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5)
print(myplot)


types=c(
  "Rat gut",
  "Rat non-gut",
  "Rat (Undefined tissue)",
  "Mouse gut",
  "Mouse non-gut",
  "Mouse (Undefined tissue)",
  "Other rodents gut",
  "Other rodents non-gut",
  "Other rodents (Undefined tissue)",
  "Human gut",
  "Human non-gut",
  "Human (Undefined tissue)",
  "Other mammals gut",
  "Other mammals non-gut",
  "Other mammals (Undefined tissue)",
  "Other animals gut",
  "Other animals non-gut",
  "Other animals (Undefined tissue)",
  "Gut (Undefined host)",
  "Non-gut (Undefined host)",
  "Food",
  "Environment"
  #"Undefined"
)

statistics_types <- concat %>%
  count(Database, Type) %>%                     # your original step
  complete(Database, Type, fill = list(n = 0))  # <-- add the missing rows

max_value=max(statistics_types$n)
myplot <- ggplot(statistics_types, aes(x = factor(Type, levels = types), y = n, fill = Database)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, max_value+10)) +
  labs(x = "Source", y = "Number of mapped species") +
  scale_fill_manual(values = c("GTDB" = "darkred", "WoL" = "darkorange"),
                    labels = c("GTDB", "WoL"),
                    name = "Database") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5)
print(myplot)

dev.off()

#max_value=max(statistics_plot$n)
#myplot <- ggplot(statistics_plot, aes(x = factor(Database, levels = c("WoL","GTDB")), y = n, fill = factor(Host, levels = order))) +
#  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
#  theme_bw() +
#  theme(panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.line.x = element_line(),
#        axis.line.y = element_line(),
#        axis.text = element_text(size = 11),
#        axis.title = element_text(size = 13),
#        legend.text = element_text(size = 11),
#        legend.title = element_text(size = 13)) +
#  coord_cartesian(ylim = c(0, max_value+10)) +
#  labs(x = "", y = "Number of mapped species")
#  #scale_fill_manual(#values = c("GTDB" = "darkgreen", "WoL" = "lightgreen"),
#  #                  labels = c("GTDB", "WoL"),
#  #                  name = "Isolation source")
#  #scale_x_discrete(guide = guide_axis(angle = 60)) +
#  #geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5)
#print(myplot)




## Function to get the BioSample accession number from a genome accession
#get_biosample_accession <- function(genome_accession) {
#  # Search for the genome accession in the assembly database
#  search_result <- entrez_search(db = "assembly", term = genome_accession)
#  
#  if (length(search_result$ids) == 0) {
#    return(NA)  # Return NA if no result is found
#  }
#  
#  # Fetch the summary of the assembly to extract the BioSample accession
#  summary <- entrez_summary(db = "assembly", id = search_result$ids[1])
#  
#  # Extract the BioSample accession (if available)
#  biosample_accession <- summary$biosampleaccn
#  
#  return(biosample_accession)
#}
#
## Function to get the BioSample isolation source using the BioSample accession
#get_biosample_isolation_source <- function(biosample_accession) {
#  # Fetch the full XML record for the BioSample with error handling
#  biosample_xml <- tryCatch({
#    entrez_fetch(db = "biosample", id = biosample_accession, rettype = "xml")
#  }, error = function(e) {
#    cat("Error in fetching BioSample for accession:", biosample_accession, "\n")
#    return(NA)  # Return NA if fetch fails
#  })
#  
#  if (is.na(biosample_xml)) {
#    return("Not Available")
#  }
#  
#  # Parse the XML content
#  xml_parsed <- xmlParse(biosample_xml)
#  
#  # Use XPath to find the isolation source attribute
#  isolation_source <- xpathSApply(xml_parsed, "//Attribute[@attribute_name='isolation source' or @attribute_name='isolation-source' or @attribute_name='isolation_source' or @attribute_name='tissue' or @attribute_name='environment']", xmlValue)
#  
#  # Return the isolation source or "Not Available" if not found
#  return(ifelse(length(isolation_source) > 0, isolation_source, "Not Available"))
#}
#
## List of genome accession numbers
##genome_accessions <- species_gtdb  # Replace with your list of accessions
#genome_accessions <- species_wol  # Replace with your list of accessions
#
## Initialize a data frame to store results
#results <- data.frame(GenomeAccession = genome_accessions, BioSampleAccession = character(length(genome_accessions)), IsolationSource = character(length(genome_accessions)), stringsAsFactors = FALSE)
#
## Loop through each genome accession and query for the BioSample isolation source
#for (i in 1:length(genome_accessions)) {
#  genome_acc <- genome_accessions[i]
#  
#  # Get the BioSample accession with retry
#  biosample_acc <- tryCatch({
#    get_biosample_accession(genome_acc)
#  }, error = function(e) {
#    cat("Error in fetching BioSample accession for genome accession:", genome_acc, "\n")
#    return(NA)
#  })
#  
#  # Get the isolation source if BioSample accession is found
#  if (!is.na(biosample_acc)) {
#    isolation_source <- get_biosample_isolation_source(biosample_acc)
#  } else {
#    isolation_source <- "Not Found"
#  }
#  
#  # Store the results
#  results$BioSampleAccession[i] <- ifelse(is.na(biosample_acc), "Not Found", biosample_acc)
#  results$IsolationSource[i] <- isolation_source
#  
#  # Print progress and save intermediate results every 100 queries
#  cat("Genome accession:", i, "of", length(genome_accessions), "/ Isolation Source:", results$IsolationSource[i], "\n")
#  if (i %% 100 == 0) {
#    write.csv(results, "intermediate_biosample_isolation_sources.csv", row.names = FALSE)
#  }
#  
#  # Pause to avoid overloading the server
#  Sys.sleep(0.5)
#}
#
#
## Save the final results to a CSV file
#write.csv(results, "biosample_isolation_sources_wol.csv", row.names = FALSE)
#
#save(results,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Matrix_Processing/biosamples_wol.RData")
#
#
########
#
#sourcers=unique(results$IsolationSource)
#write.csv(sourcers, "sources_wol.csv", row.names = FALSE)
#
##############################
#
##sources <- read.csv("/nfs/users/abaud/fmorillo/sources.csv")
#sources <- read.csv("~/NextFlow/Microbiome_Profiling/code/extra_codes/sources_wol.csv", header=FALSE)
#colnames(sources)=c("isolation_source","type")
##load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Matrix_Processing/biosamples.RData")
#load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Matrix_Processing/biosamples_wol.RData")
#results$Type=sources$type[match(results$IsolationSource,sources$isolation_source)]
#stats=count(results,Type)
#pie(stats$n,stats$Type)
#slices=stats$n
#labels=stats$Type
#percentages <- round(slices / sum(slices) * 100)
#labels <- paste(labels, percentages, "%", sep=" ")
#
## Create pie chart with colors
#pie(slices, labels = labels,col = rainbow(length(slices)))
#
#save(results,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Matrix_Processing/biosamples.RData")
#
##########
#
## Create the new column based on the presence of any of the keywords
#results$Environment <- ifelse(sapply(results$IsolationSource, function(x) any(grepl(paste(keywords_plant, collapse = "|"), x, ignore.case = TRUE))), "Plant",
#                              ifelse(sapply(results$IsolationSource, function(x) any(grepl(paste(keywords_animal_gut, collapse = "|"), x, ignore.case = TRUE))), "Animal Gut", NA))
#
######
#
#load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Matrix_Processing/biosamples.RData")
#results$Genus=taxonomy$Genus[match(results$GenomeAccession,taxonomy$Genome)]
