taxonomy_gtdb$Genome=gsub("GCA_","G",taxonomy_gtdb$Genome)
taxonomy_gtdb$Genome=gsub("GCF_","G",taxonomy_gtdb$Genome)
taxonomy_gtdb$Genome=gsub("\\.1","",taxonomy_gtdb$Genome)
taxonomy_16$GTDB=taxonomy_gtdb$Genome[match(taxonomy_16$Species, taxonomy_gtdb$Species)]
taxonomy_16=taxonomy_16[-1,]
taxonomy_16S=taxonomy_16[,c(1,2,3,10)]
write.table(taxonomy_16S,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_Greengenes2.txt",row.names=F,col.names=F,quote=F,sep="\t")


matrix = get(paste('filtered_trait',subset_id,sep='_'))
matrix = t(matrix)
colnames(matrix)=taxonomy$GTDB[match(colnames(matrix), taxonomy$Genome)]
collapsed_matrix=t(matrix)
collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
assign(paste('collapsed',rank,subset_id,sep='_'), collapsed_matrix)
