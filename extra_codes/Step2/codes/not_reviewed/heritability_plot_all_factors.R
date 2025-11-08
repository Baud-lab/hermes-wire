# make plots (pwys)
cbPalette <- c("lightblue","lightcoral", "lightgreen")

filtered_VCs=all_VCs_full[!grepl("Non-significant",all_VCs_full$Significance),]
filtered_VCs=filtered_VCs[grepl(subset_id,filtered_VCs$subset_id),]
filtered_VCs=filtered_VCs[grepl(rank,filtered_VCs$rank),]
filtered_VCs <- filtered_VCs[order(filtered_VCs$var_Ad),]
filtered_VCs$LBL[filtered_VCs$Significance == "Significant (Bonferroni)"] <- "**"
filtered_VCs$LBL[filtered_VCs$Significance == "Significant (FDR)"] <- "*"
filtered_VCs$trait <- factor(as.character(filtered_VCs$trait),levels=filtered_VCs$trait[order(filtered_VCs$var_Ad)])

filtered_VCs_plot <- filtered_VCs[,c("trait","LBL","var_Ad","var_C","var_Ed","STE_Ad")]
filtered_VCs_plotL <- gather(filtered_VCs_plot,"Var.Exp","Var.Exp.NR", var_Ad:var_Ed,factor_key = T)
filtered_VCs_plotL$Var.Exp <- as.character(filtered_VCs_plotL$Var.Exp)
filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_Ed"] <- "Environment"
filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_C"] <- "Cohousing"
filtered_VCs_plotL$Var.Exp[filtered_VCs_plotL$Var.Exp == "var_Ad"] <- "Aggregate Host Genetics"
filtered_VCs_plotL$Var.Exp <- factor(as.character(filtered_VCs_plotL$Var.Exp),level = c("Environment","Cohousing","Aggregate Host Genetics"))

ggplot(filtered_VCs_plotL,aes(x=trait,y=Var.Exp.NR,fill=Var.Exp)) + 
  scale_fill_manual(values = cbPalette,name="Effect") +
  geom_col(col="grey", width=1,size=0.5) + 
  #geom_errorbar(ymin=filtered_VCs_plotL$STE_Ad,ymax=filtered_VCs_plotL$STE_Ad,width=0.25, linetype='solid') + 
  geom_text(data = filtered_VCs_plotL,
            aes(x = trait, y=STE_Ad,
                label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
            color="black", vjust=+0.75, angle = 0, hjust=-1,size=5) + ylim(-0.01,1.01) + 
  ylab('Variance Explained') + xlab('') + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  coord_flip()
