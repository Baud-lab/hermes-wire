library(stringr)

counts<- read.delim("~/NextFlow/Microbiome_Profiling/step1/00077E61FF_S250_L001.kai")
original_matrix=as.matrix(counts)
rownames(original_matrix)=original_matrix[,1]
original_matrix=as.matrix(original_matrix[,-1])

# Assuming original_matrix is your original matrix with row names

# Function to categorize rows based on the number of substrings
categorize_rows <- function(row_names) {
  num_substrings <- sapply(strsplit(row_names, "\\."), length)
  
  EC1_rows <- row_names[num_substrings == 1]
  EC2_rows <- row_names[num_substrings == 2]
  EC3_rows <- row_names[num_substrings == 3]
  EC4_rows <- row_names[num_substrings == 4]
  
  return(list(EC1 = EC1_rows, EC2 = EC2_rows, EC3 = EC3_rows, EC4 = EC4_rows))
}

# Categorize rows
categorized_rows <- categorize_rows(rownames(original_matrix))

# Create matrices based on categorized rows
EC1 <- as.matrix(original_matrix[categorized_rows$EC1, ])
EC2 <-  as.matrix(original_matrix[categorized_rows$EC2, ])
EC3 <-  as.matrix(original_matrix[categorized_rows$EC3, ])
EC4 <-  as.matrix(original_matrix[categorized_rows$EC4, ])

# Step 2: Remove the last substring from row names in EC4 and bind it with EC3
rnames=as.factor(paste(sapply(strsplit(rownames(EC4), "\\."), "[",1),
             sapply(strsplit(rownames(EC4), "\\."), "[",2),
             sapply(strsplit(rownames(EC4), "\\."), "[",3),
             sep="."))
new_matrix <- t(matrix(EC4, nrow = nrow(EC4), ncol = ncol(EC4)))
colnames(new_matrix) <- rnames


# Step 3: Remove the last substring from row names in EC3 and bind it with EC2
EC2 <- rbind(EC2, as.matrix(EC3[ -length(strsplit(rownames(EC3), "\\.")[[1]]),]))

# Step 4: Remove the last substring from row names in EC2 and bind it with EC1
EC1 <- rbind(EC1, as.matrix(EC2[ -length(strsplit(rownames(EC2), "\\.")[[1]]),]))


# Step 5: Sum the values for each row with the same name on EC1, EC2, and EC3
collapsed_EC1 <- sapply(by(EC1, rownames(EC1), colSums), identity)
collapsed_EC2 <- sapply(by(EC2, rownames(EC2), colSums), identity)
collapsed_EC3 <- sapply(by(EC3, rownames(EC3), colSums), identity)
