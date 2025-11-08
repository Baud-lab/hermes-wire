ab_vs_prev_Genus_ALL <- ab_vs_prev_Genus_ALL %>% 
  mutate( bonf = case_when(
    Significance == "Significant (Bonferroni)" ~ "Significant (Bonferroni)",
    TRUE ~ "Non-significant traits"))



myplot<-ggplot(ab_vs_prev_Genus_ALL, aes(x=reorder(bonf,prevalence,FUN = median), y=prevalence)) +
  geom_boxplot() +
  xlab("") +
  ylab("Prevalence (%)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(unique(ab_vs_prev_Genus_ALL$bonf)),
              map_signif_level = TRUE)
print(myplot)

ggplot(aes(x = abundance, group= bonf, fill = bonf), data = ab_vs_prev_Genus_ALL) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_manual('', values = c('Non-significant traits' = 'lightcoral', 'Significant (Bonferroni)' = 'lightblue')) +
  labs(title = '', x = 'Abundance (%)', y = 'Frequence')
