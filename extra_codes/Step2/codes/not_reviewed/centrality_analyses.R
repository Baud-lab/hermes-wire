prevotella=filtered_VCs[grepl("g__Prevotella",filtered_VCs$Genus),]
plot(prevotella$Heritability,prevotella$Closeness)
cor.test(prevotella$Heritability,prevotella$Closeness,method="spearman")

bacteroides=filtered_VCs[grepl("g__Bacteroides",filtered_VCs$Genus) & !grepl("g__Bacteroides_F",filtered_VCs$Genus),]
plot(bacteroides$Heritability,bacteroides$Closeness)
cor.test(bacteroides$Heritability,bacteroides$Closeness,method="spearman")

prevotella=filtered_VCs[grepl("g__Prevotella",filtered_VCs$Genus),]
plot(prevotella$Heritability,prevotella$Betweeness)
cor.test(prevotella$Heritability,prevotella$Betweeness,method="spearman")

bacteroides=filtered_VCs[grepl("g__Bacteroides",filtered_VCs$Genus) & !grepl("g__Bacteroides_F",filtered_VCs$Genus),]
plot(bacteroides$Heritability,bacteroides$Betweeness)
cor.test(bacteroides$Heritability,bacteroides$Betweeness,method="spearman")


prevotella=filtered_VCs[grepl("g__Prevotella",filtered_VCs$Genus),]
plot(prevotella$Heritability,prevotella$Degree)
cor.test(prevotella$Heritability,prevotella$Degree,method="spearman")

bacteroides=filtered_VCs[grepl("g__Bacteroides",filtered_VCs$Genus) & !grepl("g__Bacteroides_F",filtered_VCs$Genus),]
plot(bacteroides$Heritability,bacteroides$Degree)
cor.test(bacteroides$Heritability,bacteroides$Degree,method="spearman")

prevotella=filtered_VCs[grepl("g__CAG-873",filtered_VCs$Genus),]
plot(prevotella$Heritability,prevotella$Closeness)
cor.test(prevotella$Heritability,prevotella$Closeness,method="spearman")
