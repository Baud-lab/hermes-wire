load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/old_Outputs/old_maternal/Output_Maternal/Heritability/heritability.RData")
maternal=all_VCs_full[all_VCs_full$subset_id=="ALL",]
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/old_Outputs/old_jul_24/Output_Shallow/Heritability/heritability.RData")
cage=all_VCs_full[all_VCs_full$subset_id=="ALL",]

cage=cage[rownames(cage) %in% rownames(maternal),]
maternal=maternal[rownames(maternal) %in% rownames(cage),]
maternal=maternal[rownames(cage),]
plot(cage$var_Ad,maternal$var_Ad,xlab="Without dam", ylab="With dam")
abline(lm(maternal$var_Ad ~ cage$var_Ad), col = "red")
cor_result=cor.test(cage$var_Ad,maternal$var_Ad,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
legend("topleft", legend = c(paste("R² =", round(r_value, 3)), "P < 0.01"), col = c("transparent", "transparent"), pch = 1)

x=data.frame(with_dam=maternal$var_Ad,without_dam=cage$var_Ad)

ggplot(x, aes(x = without_dam, y = with_dam)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
