stats <- read.csv("/nfs/users/abaud/fmorillo/P50/microbiome/output/Woltka/stats.tsv")
stats$prop_gtdb=100*stats$aligned_to_db/stats$available_reads

load("/nfs/users/abaud/fmorillo/P50/microbiome/output/Woltka/good_samples.RData")
good_samples=sapply(strsplit(good_samples, "_"), "[",1)
stats=stats[stats$sample %in% good_samples,]

# Load non-host reads from qiita
qiita_after_filters <- read_excel("P50/microbiome/output/Woltka/qiita_after_filters.xlsx")
qiita_after_filters=qiita_after_filters[!grepl("BLANK",qiita_after_filters$filename),]
qiita_after_filters$type=sapply(strsplit(qiita_after_filters$filename, "_"), "[",4)
qiita_after_filters$filename=sapply(strsplit(qiita_after_filters$filename, "_"), "[",1)
qiita_after_filters=qiita_after_filters[!grepl("R2",qiita_after_filters$type),]
for(i in 1:length(qiita_after_filters$filename)){
  if(!nchar(qiita_after_filters$filename[i])>=10){qiita_after_filters$filename[i] <- paste0("000", qiita_after_filters$filename[i])}
}
qiita_after_filters=qiita_after_filters[qiita_after_filters$filename %in% stats$sample,]
qiita_after_filters <- qiita_after_filters[match(stats$sample, qiita_after_filters$filename), ]
stats$available_reads_qiita=qiita_after_filters$reads
stats$available_reads_qiita_2x=2*qiita_after_filters$reads
# Load genome counts table
input_matrix="/users/abaud/fmorillo/P50/microbiome/output/Woltka/191716_none.biom"
counts_table<-as.matrix(biom_data(read_biom(input_matrix)))
#### Customized!!!!
counts_table=counts_table[,!grepl("BLANK",colnames(counts_table))]
ids <- sapply(strsplit(colnames(counts_table), "\\."), function(x) tail(x, 1))
if(!all(nchar(ids)>=10)){
  for(i in 1:length(ids)){
    if(!nchar(ids[i])>=10){ids[i] <- paste0("000", ids[i])}
  }
}
colnames(counts_table)<-ids
counts_table=counts_table[,-c(which(colnames(counts_table)=="00077E9A45"), which(colnames(counts_table)=="00077E7BC3"))]
woltka_depth=data.frame(sample=colnames(counts_table),reads=colSums(counts_table))
woltka_depth <- woltka_depth[match(stats$sample, woltka_depth$sample), ]
stats$aligned_to_wol=woltka_depth$reads
stats$prop_wol=100*stats$aligned_to_wol/stats$available_reads_qiita_2x

gtdb=data.frame(value=stats$prop_gtdb,catalogue="GTDB with Kaiju")
wol=data.frame(value=stats$prop_wol,catalogue="WoL with Woltka")
db=rbind(gtdb,wol)


# Plot the histogram
pdf("/users/abaud/fmorillo/paper_figures/gtdb_vs_wol.pdf")
myplot=ggplot(db, aes(x = value, fill = catalogue)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(fill = "Profiling Method") +
  labs(title = "",
       x = "Mapping ratio (%)",
       y = "Frequency") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        legend.position = "bottom")
print(myplot)
dev.off()
