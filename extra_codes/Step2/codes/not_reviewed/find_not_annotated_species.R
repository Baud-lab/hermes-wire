matrix_tax=filtered_prev_objs[[18]]
x=data.frame(species=rownames(matrix_tax),
             genus=taxonomy$Genus[match(rownames(matrix_tax),taxonomy$Species)],
             genome=taxonomy$Genome[match(rownames(matrix_tax),taxonomy$Species)])
x=x[x$genus %in% genera, ]

func=rownames(matrix)
tax=x$species
not_in_func <- setdiff(tax, func)
genomes=taxonomy$Genome[match(not_in_func,taxonomy$Species)]
write.table(genomes,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/genomes_v3.txt",quote=F,row.names = F,col.names = F)
