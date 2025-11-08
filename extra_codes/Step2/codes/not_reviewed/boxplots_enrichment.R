df=data.frame(species=rownames(presence_absence),M00153=presence_absence[,which(colnames(presence_absence)=="M00153")])
df$Presence[df$M00153==1]="Present"
df$Presence[df$M00153==0]="Absent"
df$heritability=filtered_VCs$Heritability[match(df$species,filtered_VCs$trait)]

order=c("Present", "Absent")
myplot<-ggplot(df, aes(x=factor(Presence, levels = order), y=heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Present", "Absent")),
              map_signif_level = TRUE)
print(myplot)
