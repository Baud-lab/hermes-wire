load("/users/abaud/fmorillo/paper_figures/enrichment_t.test_final.RData")

# PFAM
sig_pfam_cag485_high=correlations_total_CAG485_PFAMs$trait[correlations_total_CAG485_PFAMs$Significance!="Non-significant" & correlations_total_CAG485_PFAMs$Direction=="High heritability"]
sig_pfam_cag485_low=correlations_total_CAG485_PFAMs$trait[correlations_total_CAG485_PFAMs$Significance!="Non-significant" & correlations_total_CAG485_PFAMs$Direction=="Low heritability"]
sig_pfam_acetatifactor_high=correlations_total_Acetatifactor_PFAMs$trait[correlations_total_Acetatifactor_PFAMs$Significance!="Non-significant" & correlations_total_Acetatifactor_PFAMs$Direction=="High heritability"]
sig_pfam_acetatifactor_low=correlations_total_Acetatifactor_PFAMs$trait[correlations_total_Acetatifactor_PFAMs$Significance!="Non-significant" & correlations_total_Acetatifactor_PFAMs$Direction=="Low heritability"]
shared_pfam_high=intersect(sig_pfam_cag485_high,sig_pfam_acetatifactor_high)
shared_pfam_low=intersect(sig_pfam_cag485_low,sig_pfam_acetatifactor_low)

# Modules
sig_module_cag485_high=correlations_total_CAG485_KEGG_Module$trait[correlations_total_CAG485_KEGG_Module$Significance!="Non-significant" & correlations_total_CAG485_KEGG_Module$Direction=="High heritability"]
sig_module_cag485_low=correlations_total_CAG485_KEGG_Module$trait[correlations_total_CAG485_KEGG_Module$Significance!="Non-significant" & correlations_total_CAG485_KEGG_Module$Direction=="Low heritability"]
sig_module_acetatifactor_high=correlations_total_Acetatifactor_KEGG_Module$trait[correlations_total_Acetatifactor_KEGG_Module$Significance!="Non-significant" & correlations_total_Acetatifactor_KEGG_Module$Direction=="High heritability"]
sig_module_acetatifactor_low=correlations_total_Acetatifactor_KEGG_Module$trait[correlations_total_Acetatifactor_KEGG_Module$Significance!="Non-significant" & correlations_total_Acetatifactor_KEGG_Module$Direction=="Low heritability"]
shared_module_high=intersect(sig_module_cag485_high,sig_module_acetatifactor_high)
shared_module_low=intersect(sig_module_cag485_low,sig_module_acetatifactor_low)

#[Shared High] "M00004" (Pentose phosphate pathway) "M00007" (Pentose phosphate pathway, non-oxidative phase)

# COGs
sig_cog_cag485_high=correlations_total_CAG485_eggNOG_OGs$COG[correlations_total_CAG485_eggNOG_OGs$Significance!="Non-significant" & correlations_total_CAG485_eggNOG_OGs$Direction=="High heritability"]
sig_cog_cag485_low=correlations_total_CAG485_eggNOG_OGs$COG[correlations_total_CAG485_eggNOG_OGs$Significance!="Non-significant" & correlations_total_CAG485_eggNOG_OGs$Direction=="Low heritability"]
sig_cog_acetatifactor_high=correlations_total_Acetatifactor_eggNOG_OGs$COG[correlations_total_Acetatifactor_eggNOG_OGs$Significance!="Non-significant" & correlations_total_Acetatifactor_eggNOG_OGs$Direction=="High heritability"]
sig_cog_acetatifactor_low=correlations_total_Acetatifactor_eggNOG_OGs$COG[correlations_total_Acetatifactor_eggNOG_OGs$Significance!="Non-significant" & correlations_total_Acetatifactor_eggNOG_OGs$Direction=="Low heritability"]
shared_cog_high=intersect(sig_cog_cag485_high,sig_cog_acetatifactor_high)
shared_cog_low=intersect(sig_cog_cag485_low,sig_cog_acetatifactor_low)


##########


# PFAM
sig_pfam_cag485_high=all_cogs_CAG485_PFAMs$COG[all_cogs_CAG485_PFAMs$Significance!="Non-significant" & all_cogs_CAG485_PFAMs$Direction=="High heritability"]
sig_pfam_cag485_low=all_cogs_CAG485_PFAMs$COG[all_cogs_CAG485_PFAMs$Significance!="Non-significant" & all_cogs_CAG485_PFAMs$Direction=="Low heritability"]
sig_pfam_acetatifactor_high=all_cogs_Acetatifactor_PFAMs$COG[all_cogs_Acetatifactor_PFAMs$Significance!="Non-significant" & all_cogs_Acetatifactor_PFAMs$Direction=="High heritability"]
sig_pfam_acetatifactor_low=all_cogs_Acetatifactor_PFAMs$COG[all_cogs_Acetatifactor_PFAMs$Significance!="Non-significant" & all_cogs_Acetatifactor_PFAMs$Direction=="Low heritability"]
shared_pfam_high=intersect(sig_pfam_cag485_high,sig_pfam_acetatifactor_high)
shared_pfam_low=intersect(sig_pfam_cag485_low,sig_pfam_acetatifactor_low)

#########

# Modules
sig_module_cag485_high=all_cogs_CAG485_KEGG_Module$COG[all_cogs_CAG485_KEGG_Module$Significance!="Non-significant" & all_cogs_CAG485_KEGG_Module$Direction=="High heritability"]
sig_module_cag485_low=all_cogs_CAG485_KEGG_Module$COG[all_cogs_CAG485_KEGG_Module$Significance!="Non-significant" & all_cogs_CAG485_KEGG_Module$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_KEGG_Module$COG[all_cogs_Prevotella_KEGG_Module$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_Module$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_KEGG_Module$COG[all_cogs_Prevotella_KEGG_Module$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_Module$Direction=="Low heritability"]
shared_module_high=intersect(sig_module_cag485_high,sig_module_Prevotella_high)
shared_module_low=intersect(sig_module_cag485_low,sig_module_Prevotella_low)

# KOs
sig_module_cag485_high=all_cogs_CAG485_KEGG_ko$COG[all_cogs_CAG485_KEGG_ko$Significance!="Non-significant" & all_cogs_CAG485_KEGG_ko$Direction=="High heritability"]
sig_module_cag485_low=all_cogs_CAG485_KEGG_ko$COG[all_cogs_CAG485_KEGG_ko$Significance!="Non-significant" & all_cogs_CAG485_KEGG_ko$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_KEGG_ko$COG[all_cogs_Prevotella_KEGG_ko$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_ko$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_KEGG_ko$COG[all_cogs_Prevotella_KEGG_ko$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_ko$Direction=="Low heritability"]
shared_ko_high=intersect(sig_module_cag485_high,sig_module_Prevotella_high)
shared_ko_low=intersect(sig_module_cag485_low,sig_module_Prevotella_low)

# COGs
sig_module_cag485_high=all_cogs_CAG485_eggNOG_OGs$COG[all_cogs_CAG485_eggNOG_OGs$Significance!="Non-significant" & all_cogs_CAG485_eggNOG_OGs$Direction=="High heritability"]
sig_module_cag485_low=all_cogs_CAG485_eggNOG_OGs$COG[all_cogs_CAG485_eggNOG_OGs$Significance!="Non-significant" & all_cogs_CAG485_eggNOG_OGs$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_eggNOG_OGs$COG[all_cogs_Prevotella_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Prevotella_eggNOG_OGs$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_eggNOG_OGs$COG[all_cogs_Prevotella_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Prevotella_eggNOG_OGs$Direction=="Low heritability"]
shared_cog_high=intersect(sig_module_cag485_high,sig_module_Prevotella_high)
shared_cog_low=intersect(sig_module_cag485_low,sig_module_Prevotella_low)

##########

# Modules
sig_module_Bacteroides_high=all_cogs_Bacteroides_KEGG_Module$COG[all_cogs_Bacteroides_KEGG_Module$Significance!="Non-significant" & all_cogs_Bacteroides_KEGG_Module$Direction=="High heritability"]
sig_module_Bacteroides_low=all_cogs_Bacteroides_KEGG_Module$COG[all_cogs_Bacteroides_KEGG_Module$Significance!="Non-significant" & all_cogs_Bacteroides_KEGG_Module$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_KEGG_Module$COG[all_cogs_Prevotella_KEGG_Module$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_Module$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_KEGG_Module$COG[all_cogs_Prevotella_KEGG_Module$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_Module$Direction=="Low heritability"]
shared_module_high=intersect(sig_module_Bacteroides_high,sig_module_Prevotella_high)
shared_module_low=intersect(sig_module_Bacteroides_low,sig_module_Prevotella_low)

# KOs
sig_module_Bacteroides_high=all_cogs_Bacteroides_KEGG_ko$COG[all_cogs_Bacteroides_KEGG_ko$Significance!="Non-significant" & all_cogs_Bacteroides_KEGG_ko$Direction=="High heritability"]
sig_module_Bacteroides_low=all_cogs_Bacteroides_KEGG_ko$COG[all_cogs_Bacteroides_KEGG_ko$Significance!="Non-significant" & all_cogs_Bacteroides_KEGG_ko$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_KEGG_ko$COG[all_cogs_Prevotella_KEGG_ko$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_ko$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_KEGG_ko$COG[all_cogs_Prevotella_KEGG_ko$Significance!="Non-significant" & all_cogs_Prevotella_KEGG_ko$Direction=="Low heritability"]
shared_ko_high=intersect(sig_module_Bacteroides_high,sig_module_Prevotella_high)
shared_ko_low=intersect(sig_module_Bacteroides_low,sig_module_Prevotella_low)

# COGs
sig_module_Bacteroides_high=all_cogs_Bacteroides_eggNOG_OGs$COG[all_cogs_Bacteroides_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Bacteroides_eggNOG_OGs$Direction=="High heritability"]
sig_module_Bacteroides_low=all_cogs_Bacteroides_eggNOG_OGs$COG[all_cogs_Bacteroides_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Bacteroides_eggNOG_OGs$Direction=="Low heritability"]
sig_module_Prevotella_high=all_cogs_Prevotella_eggNOG_OGs$COG[all_cogs_Prevotella_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Prevotella_eggNOG_OGs$Direction=="High heritability"]
sig_module_Prevotella_low=all_cogs_Prevotella_eggNOG_OGs$COG[all_cogs_Prevotella_eggNOG_OGs$Significance!="Non-significant" & all_cogs_Prevotella_eggNOG_OGs$Direction=="Low heritability"]
shared_cog_high=intersect(sig_module_Bacteroides_high,sig_module_Prevotella_high)
shared_cog_low=intersect(sig_module_Bacteroides_low,sig_module_Prevotella_low)
