subset_id="ALL"
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


# Load clusters
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Beta_TD_PD/Cluster_Analyses/filtered_prev_clusters.RData")
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
clusters_ALL=filtered_prev_cluster_objs[[number]]

# Load ECs
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input_func/EC_ALL_RES.RData")
residuals_ECs=residuals_ECs

clusters_ALL=clusters_ALL[,colnames(clusters_ALL) %in% colnames(residuals_ECs)]
clusters_ALL=clusters_ALL[,colnames(residuals_ECs)]

matrix1=clusters_ALL
matrix2=residuals_ECs

correlation_matrix <- cor(t(matrix1), t(matrix2), method="spearman")

cluster="cluster__11"

cor_cluster_10 <- correlation_matrix["cluster__10", ]
cor_cluster_11 <- correlation_matrix["cluster__11", ]
cor_cluster_3 <- correlation_matrix["cluster__3", ]
cor_cluster_13 <- correlation_matrix["cluster__13", ]
cor_cluster_5 <- correlation_matrix["cluster__5", ]
other_clusters <- correlation_matrix[rownames(correlation_matrix) != cluster, ]
mean_cor_other_clusters <- colMeans(other_clusters)

cor_diff <- cor_cluster_11 - mean_cor_other_clusters
cor_diff <- cor_cluster_10 - cor_cluster_11

names(cor_cluster_10) <- colnames(correlation_matrix)
names(cor_cluster_11) <- colnames(correlation_matrix)
names(cor_cluster_3) <- colnames(correlation_matrix)
names(cor_cluster_13) <- colnames(correlation_matrix)
names(cor_cluster_5) <- colnames(correlation_matrix)

names(cor_diff) <- colnames(correlation_matrix)
diff=as.data.frame(cor_diff)

cluster_10=as.data.frame(cor_cluster_10)
cluster_11=as.data.frame(cor_cluster_11)
cluster_3=as.data.frame(cor_cluster_3)
cluster_13=as.data.frame(cor_cluster_13)
cluster_5=as.data.frame(cor_cluster_5)

comp=cbind(cluster_10,cluster_11,cluster_3,cluster_13,cluster_5)

# Perform paired t-test
t_test_result <- t.test(cor_cluster_10, cor_cluster_11, paired = TRUE)
print(t_test_result)
# Perform Wilcoxon signed-rank test
wilcox_test_result <- wilcox.test(cor_cluster_10, cor_cluster_11, paired = TRUE)
print(wilcox_test_result)

