library(factoextra)

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


# 4. Data transformation (CLR)
calculate_gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}
do_clr_default = function(data){
  transform_col = function(x){
    gMean = calculate_gm_mean(x)
    transformed_col = log(x / gMean)
    return(transformed_col)
  }
  transformed_data = apply(data, MAR = 2, FUN = transform_col)
  return(transformed_data)
}

subset_id="ALL"
rank="Species"
method="Residuals"

if (rank=="Species"){
  index=18
} else {
  if (rank=="Genus") {
    index=15
  } else {
    if (rank=="Phylum"){
      index=3
    } else {
      if (rank=="Family") {
        index=12
      }
    }
  }
}


if (method=="Residuals"){
  #Residuals
  load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/residuals_qned_counts_new.RData")
  residuals=residuals_qned_counts_objs[[index]]
  matrix=t(residuals)
} else {
  #Prevalence
  load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/filtered_prev_new.RData")
  filtered_prev=filtered_prev_objs[[index]]
  offset = 0.00001
  filtered_prev = filtered_prev + offset
  clr_counts = do_clr_default(filtered_prev)
  matrix=t(clr_counts)
}

#optimal_enterotypes <-2
#ent_rank= "Genus"
ent_kmax= 5
ent_nstart= 25
# Determine the optimal number of enterotypes using the silhouette analysis
nbclust_plot <- fviz_nbclust(matrix, kmeans, method = "silhouette", k.max = ent_kmax)
#plot(nbclust_plot)
silhouette_results <- fviz_nbclust(matrix, kmeans, method = "silhouette", k.max = ent_kmax)
silhouette=silhouette_results$data
silhouette$abs_diff_y <- c(NA, abs(diff(silhouette$y)))
optimal_enterotypes <- which.max(silhouette$abs_diff_y)

## Define enterotypes
set.seed(123)
kmeans_results <- kmeans(matrix, centers = optimal_enterotypes, nstart = ent_nstart)

## Plot enterotypes
fviz_cluster_plot <- fviz_cluster(kmeans_results, data = matrix,
                                  label = NULL,
                                  geom = "point",
                                  palette = c("blue","orange"),
                                  pointsize = 3,
                                  ellipse = TRUE,
                                  show.clust.cent = TRUE,
                                  ellipse.type = "convex",
                                  ggtheme = theme_minimal()) +
  ggtitle(paste0("Cluster Plot with Ellipses: ",subset_id)) +
  theme(plot.title = element_text(hjust = 0.5))
plot(fviz_cluster_plot)
#
#
## Centroids of the clusters
#kmeans_centers <- kmeans_results$centers
#print(kmeans_centers)
#
## Calculate differences between cluster centroids
#centroid_differences <- abs(kmeans_centers[1, ] - kmeans_centers[2, ])
#print(centroid_differences)
#
## Identify the top features driving the separation
#main_features <- sort(centroid_differences, decreasing = TRUE)
#print(main_features)[1:5]


# Perform PCA
pca_results <- prcomp(matrix, scale. = TRUE)

## Visualize the PCA with clusters
#fviz_pca_ind(pca_results, geom.ind = "point", col.ind = as.factor(kmeans_results$cluster),palette = c("blue","orange")) +
#  ggtitle(paste0("PCA Plot: ", subset_id)) +
#  theme(plot.title = element_text(hjust = 0.5))

# Contributions of variables to the first principal component
pca_contributions <- abs(pca_results$rotation[, 1])  # Contributions to PC1
main_features_pca <- sort(pca_contributions, decreasing = TRUE)
print(main_features_pca)[1:5]

# Get the contributions of variables to the first principal component (PC1)
pc1_contributions <- pca_results$rotation[, 1]
pc1_drivers=as.data.frame(pc1_contributions)
pc1_drivers$Species=rownames(pc1_drivers)
pc1_drivers$Phylum=taxonomy$Phylum[match(pc1_drivers$Species,taxonomy$Species)]
pc1_drivers$Class=taxonomy$Class[match(pc1_drivers$Species,taxonomy$Species)]
pc1_drivers$Order=taxonomy$Order[match(pc1_drivers$Species,taxonomy$Species)]
pc1_drivers$Family=taxonomy$Family[match(pc1_drivers$Species,taxonomy$Species)]
pc1_drivers$Genus=taxonomy$Genus[match(pc1_drivers$Species,taxonomy$Species)]

# Identify the most relevant positive and negative features
most_relevant_positive_feature <- names(which.max(pc1_contributions))
most_relevant_negative_feature <- names(which.min(pc1_contributions))

print(paste("Most relevant positive feature:", most_relevant_positive_feature))
print(paste("Most relevant negative feature:", most_relevant_negative_feature))

# Create a biplot, showing only the most relevant positive and negative features
fviz_pca_biplot(pca_results, 
                geom.ind = "point",  # Points for individuals (samples)
                col.ind = as.factor(kmeans_results$cluster),  # Color by clusters
                palette = c("blue","orange"),  # Color palette
                #addEllipses = TRUE,  # Add ellipses for clusters
                label = "var",  # Label variables (features)
                select.var = list(name = c(most_relevant_positive_feature, most_relevant_negative_feature)),  # Show only the most relevant features
                col.var = "black",  # Set vector arrows to black
                labelsize = 5,  # Adjust label size if needed
                col.circle = "black",  # Set vector labels to black
                repel = TRUE,  # Avoid text overlap
                ggtheme = theme_minimal())

########

# Read centrality values
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Cage/Matrix_Processing/networks.RData")
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
centrality=centralities_objs[[number]]
hubs=hubs_objs[[number]]
hubs$hub="Hubs"




load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Heritability/heritability.RData")
filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
filtered_VCs <- filtered_VCs %>% 
  mutate( All_Sig = case_when(
    Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Heritable",
    TRUE ~ "Non-heritable"))
filtered_VCs$Heritability=100*filtered_VCs$var_Ad
filtered_VCs$pc1_contributions=abs(pc1_contributions[match(filtered_VCs$trait,names(pc1_contributions))])
filtered_VCs$Degree=as.numeric(centrality$Degree[match(filtered_VCs$trait,centrality$Species)])
filtered_VCs$Betweeness=as.numeric(centrality$Betweeness[match(filtered_VCs$trait,centrality$Species)])
filtered_VCs$Closeness=as.numeric(centrality$Closeness[match(filtered_VCs$trait,centrality$Species)])


plot(filtered_VCs$Heritability,filtered_VCs$pc1_contributions)
cor.test(filtered_VCs$Heritability,filtered_VCs$pc1_contributions,method="spearman")
plot(filtered_VCs$Degree,filtered_VCs$pc1_contributions)
cor.test(filtered_VCs$Degree,filtered_VCs$pc1_contributions,method="spearman")
plot(filtered_VCs$Closeness,filtered_VCs$pc1_contributions)
cor.test(filtered_VCs$Closeness,filtered_VCs$pc1_contributions,method="spearman")
