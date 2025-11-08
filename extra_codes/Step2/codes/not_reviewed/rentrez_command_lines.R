
script="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/code/extra_codes/rentrez.py"
input="/users/abaud/fmorillo/paper_figures/species_gtdb.csv"
accession_column_name="Genome"
out="/users/abaud/fmorillo/paper_figures/gtdb_classified.csv"
python $script --input_csv $input --accession_column $accession_column_name --output_csv $out


#######


script="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/code/extra_codes/rentrez.py"
input="/users/abaud/fmorillo/paper_figures/species_wol.csv"
accession_column_name="Genome"
out="/users/abaud/fmorillo/paper_figures/wol_classified.csv"
python $script --input_csv $input --accession_column $accession_column_name --output_csv $out


#########

script="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/code/extra_codes/rentrez.py"
input="/users/abaud/fmorillo/paper_figures/species_test.csv"
accession_column_name="Genome"
out="/users/abaud/fmorillo/paper_figures/test_classified.csv"
python $script --input_csv $input --accession_column $accession_column_name --output_csv $out
