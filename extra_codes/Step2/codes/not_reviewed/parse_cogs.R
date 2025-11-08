cog.24.def <- read.delim("/nfs/users/abaud/fmorillo/paper_figures/cog-24.def.tab", header=FALSE)
enogs <- read.delim("/nfs/users/abaud/fmorillo/paper_figures/1_annotations.tsv", header=FALSE)
colnames(enogs)=c("V0","V1","V2","V3")
og_description <- read.delim2("/nfs/users/abaud/fmorillo/paper_figures/og_description.txt", header=FALSE)
og.to.category <- read.delim("/nfs/users/abaud/fmorillo/paper_figures/og-to-category.map", header=FALSE)
ogs=data.frame(V1=og_description$V1, V3=og_description$V2)
ogs$V2=og.to.category$V2[match(ogs$V1,og.to.category$V1)]
total=rbind(cog.24.def[,1:3],enogs[,-1],ogs)
significant_cogs_annotated_FINAL <- read_excel("paper_figures/orthogroup_summary_with_flags.xlsx")
significant_cogs_annotated_FINAL$`Official Name`=total$V3[match(significant_cogs_annotated_FINAL$Function,total$V1)]
significant_cogs_annotated_FINAL$Category=cog.24.def$V5[match(significant_cogs_annotated_FINAL$Function,cog.24.def$V1)]
significant_cogs_annotated_FINAL$Bioprocess=total$V2[match(significant_cogs_annotated_FINAL$Function,total$V1)]
write.table(significant_cogs_annotated_FINAL,file="~/paper_figures/significant_cogs_annotated_FINAL.tsv",sep="\t",quote=F,row.names=F,col.names=T)

og_description_filt=og_description[og_description$V1 %in% significant_cogs_annotated_FINAL$Function,]
og.to.category_filt=og.to.category[og.to.category$V1 %in% significant_cogs_annotated_FINAL$Function,]
