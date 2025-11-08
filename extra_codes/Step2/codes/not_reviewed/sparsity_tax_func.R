methods=c("Species","EC4")

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_No_Alpha_No_Net/Matrix_Processing/filtered_prev.RData")
filtered_prev_Species=filtered_prev_objs[[6]]
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_No_Alpha_No_Net/Matrix_Processing/filtered_prev.RData")
filtered_prev_EC4=filtered_prev_objs[[4]]

pdf("/users/abaud/fmorillo/paper_figures/tax_func_proportion_of_zeroes.pdf")

    proportion_zeroes_all=data.frame()
    for (method in methods){
      sample_harm=get(paste('filtered_prev',method,sep="_"))
      proportion_zeroes <- apply(sample_harm, 2, function(column) {
        sum(column == 0) / length(column)
      })
      proportion_zeroes=data.frame(sample=colnames(sample_harm),proportion=100*proportion_zeroes)
      proportion_zeroes$method=method
      proportion_zeroes_all=rbind(proportion_zeroes_all,proportion_zeroes)
    }
    p_value <- wilcox.test(proportion_zeroes_all$proportion[proportion_zeroes_all$method == "Species"], 
                           proportion_zeroes_all$proportion[proportion_zeroes_all$method == "EC4"])$p.value
    maxprop=max(proportion_zeroes_all$proportion)
    order=c("Species","EC4")
    myplot<-ggplot(proportion_zeroes_all, aes(x=factor(method, levels = order), y=proportion)) +
      geom_boxplot() +
      xlab("") +
      ylab("Proportion of zeroes (%)") +
      #ggtitle(paste0(title,": ",rank)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 14),
            legend.text=element_text(size=12),legend.title=element_text(size=15),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 60)) +
      geom_signif(comparisons = list(c("Species", "EC4")), 
                  annotations = sprintf("p = %.2e", p_value),
                  test = "wilcox.test",
                  y_position = (maxprop+10), map_signif_level = FALSE)
    print(myplot)
    myplot = ggplot(proportion_zeroes_all, aes(x = proportion, color = method, fill = method)) +
      geom_density(alpha = 0.5) +
      labs(fill = "Target", color = "Target") +
      labs(title = "",
           x = "Proportion of zeroes (%)",
           y = "Density") +
      scale_fill_manual(values = c("blue", "red")) +
      scale_color_manual(values = c("blue", "red")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            legend.position = "bottom")
    print(myplot)
    #print(paste0(title,": ",rank))
dev.off()
