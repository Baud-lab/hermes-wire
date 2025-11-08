#!/bin/bash
#SBATCH --job-name=gtdb_fasta
#SBATCH --open-mode=append
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --qos=short

### Install constainers ###

cd singularity_containers

singularity pull docker://quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0
singularity pull docker://quay.io/biocontainers/fastx_toolkit:0.0.14--h503566f_13
singularity pull docker://quay.io/biocontainers/kaiju:1.9.2--h5b5514e_0
cd ../

### Step 1: Make the concatenated fasta file ###

GTDB="../input/gtdb207_demo/"

echo "Started"
#Extract files
cd ${GTDB}
tar xvzf protein_fastas.tar.gz
echo "Extraction done"
cd protein_fastas
mv * ../
cd ../
rm -rf protein_fastas
echo "Number of files:"
echo $(ls | grep "protein.faa" | wc -l)

#Get all genome IDs from the gtdb fasta files
ls *.faa | xargs -n1 basename | sed 's/GB_//' | sed 's/RS_//' | sed 's/_protein.faa//' | sort > gtdb_taxids.tsv
echo "GTDB unique taxids file done"

#Replace the protein headers of the files by the GTDB ID (file name)
##Put all catalogue sequences in just one line
mkdir gtdb_ids
for INPUT in $(ls *.faa | sort)
do
OUTPUT=$(echo $INPUT | xargs -n1 basename)
cat ${INPUT} | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > gtdb_ids/${OUTPUT}.temp
done
echo "Temp files with GTDB IDs done"
##Assign GTDB IDs to the catalogue
for INPUT in $(ls gtdb_ids/*.temp | sort)
do
OUTPUT=$(echo $INPUT | xargs -n1 basename | sed 's/_protein.faa.temp//' )
cat ${INPUT} | perl -pe 'if(/\>/){s/\n/\t/; s/>//; s/\n//}' | awk -F '\t' -v id="$OUTPUT" '{$1=id; print}' | awk '{print ">"$1"\n"$2}' > gtdb_ids/${OUTPUT}_gtdb.faa
rm ${INPUT}
done
echo "Fasta file with GTDB IDs done"

#Concatenate files
cat gtdb_ids/*faa > gtdb_gene_catalogue.faa
echo "Files concatenated"
rm -rf gtdb_ids *protein.faa gtdb_taxids.tsv

#Split sequences
singularity run ../singularity_containers/quay.io-biocontainers-seqkit-2.6.1--h9ee0642_0.img
PERL="../codes/split.pl"
seqkit seq -M 10000 -g gtdb_gene_catalogue.faa -o under_10k.fasta
seqkit seq -m 10001 -g gtdb_gene_catalogue.faa -o over_10k.fasta
exit
cat over_10k.fasta | perl ${PERL} > split_sequences.fasta
cat split_sequences.fasta under_10k.fasta > gtdb_kaiju.tmp
echo "split sequences done"

# Remove characters that should not be on the sequences and spaces + Put in lines with 60
singularity run ../singularity_containers/fastx_toolkit_0.0.14--h503566f_13.sif
cat gtdb_kaiju.tmp | perl -lpe 's/[^ACDEFGHIKLMNPQRSTVWY]//g unless /^>/' | sed 's/ //g' | fasta_formatter -w 60 > gtdb_kaiju.faa
exit
echo "fasta file edition done"
rm under_10k.fasta over_10k.fasta split_sequences.fasta gtdb_kaiju.tmp gtdb_gene_catalogue.faa

### Step 2: Build Kaiju index ###

singularity run ../singularity_containers/kaiju_1.9.2--h5b5514e_0.sif
kaiju-mkbwt -a ACDEFGHIKLMNPQRSTVWY -o gtdb_index gtdb_kaiju.faa
kaiju-mkfmi ${DIR}gtdb_index
exit
echo "Index was created"
rm gtdb_index.bwt gtdb_index.sa gtdb_kaiju.faa

