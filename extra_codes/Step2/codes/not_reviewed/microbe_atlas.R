# 3.1 Load tables
env_info <- read.delim("/users/abaud/data/secondary/indexes/microbeatlas/samples.env.info", header=F)  
otu_info <- read.delim("/users/abaud/data/secondary/indexes/microbeatlas/otus.99.allinfo", header=F)

# 3.2 Summarize counts per OTU
# otu_info has columns: TaxonID, TotalSamples, SoilCounts, AquaticCounts, AnimalCounts, PlantCounts :contentReference[oaicite:6]{index=6}

otu_env <- otu_info %>%
  mutate(
    P_soil   = SoilCounts   / TotalSamples,
    P_aquatic= AquaticCounts/ TotalSamples,
    P_animal = AnimalCounts / TotalSamples,
    P_plant  = PlantCounts  / TotalSamples
  )

# 3.3 Map your assembly list to species names
library(rentrez)
assembly_list <- c("GCA_000007325.1", "GCF_910593845.1")
taxids <- sapply(assembly_list, function(acc) {
  esum <- entrez_summary(db="assembly", id=acc)
  esum$TaxId
})
species <- sapply(taxids, function(tid) {
  tx <- entrez_summary(db="taxonomy", id=tid)
  tx$ScientificName
})

# 3.4 Join species to OTU table by name (assuming OTU table has a Name column)
results <- data.frame(Assembly=assembly_list, TaxId=taxids, Species=species) %>%
  left_join(otu_env, by = c("Species" = "Name"))

# 3.5 Extract the environment proportions
final <- results %>%
  select(Assembly, Species, P_soil, P_aquatic, P_animal, P_plant)
print(final)
