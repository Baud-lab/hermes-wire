statistics <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input_AIE/statistics.txt")
colnames(statistics)
colnames(statistics)=c("Sample",
                       "Raw_Reads",
                       "Reads_After_Trimming",
                       "Surviving_Rate_Trimming",
                       "Reads_Aligned_to_Host",
                       "Alignment_Rate",
                       "Non_Host_Reads",
                       "Mapped_Reads")
statistics$Surviving_Rate_Alignment=100*statistics$Non_Host_Reads/statistics$Reads_After_Trimming
statistics$Raw_Mapping_Ratio=100*statistics$Mapped_Reads/statistics$Non_Host_Reads
statistics=statistics[,colnames(new)]
write.table(statistics,file="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input_AIE/statistics.txt",sep="\t",row.names = F, col.names = T)
