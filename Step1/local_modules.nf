/*
* Create a channel for tool options
*/

def getParameters(pars_tools_file) {
	def pars_tools = file(pars_tools_file)
	if( !pars_tools.exists() ) exit 1, "Missing tools options config: '$pars_tools'"

	def progPars = [:]
	def allLines  = pars_tools.readLines()

	for( line : allLines ) {
    	def list = line.split("\t")
    	if (list.length <3) {
			 error "ERROR!!! Tool option file has to be tab separated\n" 
		}
    	if (!(list[0] =~ /#/ )) {
			progPars["${list[0]}--${list[1]}"] = list[2].replace("\"", "").replace('$projectDir', "${projectDir}").replace('${projectDir}', "${projectDir}")
    	}  
	}	
	return(progPars)
}

// Processes

// 1. Taxonomic mapping

// 2. Functional mapping

// 2.2 Functional mapping
process prromenade_index {
    
    tag { "${reference}" }
    container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0"

    memory '200 GB'
    cpus 8
    time '48h'
    
    input:
    path(reference)

    output:
    path("${reference}.*.bt2l")
    
    """
    bowtie2-build --threads ${task.cpus} -f ${reference} ${reference}
    """
}

process prromenade_classifier {

	tag { pair_id }
	container 'morillofe/seqkit-prromenade:1.4'
	
	memory '50 GB'
  cpus 8
  time '6h'
    
  input:
    tuple val(pair_id), path(reads), path(index), path(ec_numbers)
  
  output:
    tuple val(pair_id), path("${pair_id}.out")
	
    script:
    def indexname = "${index[0]}".takeWhile{ it != '.' }
    """
    # 1. Concatenate fastqc pairs
    zcat ${reads[0]} ${reads[1]} > ${pair_id}_concat.fq
    # 2. Turn fastqc into fasta
    seqkit fq2fa ${pair_id}_concat.fq -o ${pair_id}.fa
    rm ${pair_id}_concat.fq
    # 3. Run PRROMenade
    prromenade-classifier  -i ${indexname} ${pair_id}.fa -t ${task.cpus} -o ${pair_id}.bin
    rm ${pair_id}.fa
    # 4. Edit final report so it can have the mappings
    cat ${pair_id}.bin | awk -F '\\t' -v OFS='\\t' '{gsub(/:[^;]*/,"",\$2); gsub(/:/,"",\$2); gsub(/;/,",",\$2); if (\$2 != "1") print}' > ${pair_id}.tmp
    rm ${pair_id}.bin
    awk 'BEGIN { FS="\\t"; OFS="\\t" }
        NR==FNR { map[\$1] = \$2; next }
        {
            split(\$2, arr, ",");
            for (i = 1; i <= length(arr); i++) {
                arr[i] = map[arr[i]];
            }
            \$2 = join(arr, ",");
            print;
        }
        function join(arr, sep, start, end, result) {
            start = start ? start : 1;
            end = end ? end : length(arr);
            for (i = start; i <= end; i++) {
                result = result ? result sep arr[i] : arr[i];
            }
            return result;
        }
    ' ${ec_numbers} ${pair_id}.tmp > temp_file && mv temp_file ${pair_id}.tmp
    cat ${pair_id}.tmp | awk -F '\\t' -v OFS='\\t' '{if (\$2 != "0") print}' > ${pair_id}.out
    rm ${pair_id}.tmp
    """

}

process humann_profiling {

	tag { pair_id }
	container 'quay.io/biocontainers/humann:3.8--pyh7cba7a3_0'
	
	memory '50 GB'
  cpus 8
  time '24h'
    
  input:
    tuple val(pair_id), path(reads), val(index), val(utility), val(classification)
  
  output:
    tuple val(pair_id),
    path("${pair_id}.out")
	
    script:
    """
    # 1. Concatenate fastqc pairs
    zcat ${reads[0]} ${reads[1]} > ${pair_id}.fastq
    echo "Concatenation Done"
    
    # 2. Run Humann3 Alignment
    mkdir ${pair_id}_profiles
    humann -t ${task.cpus} \
      ${params.EXTRAPARS} \
      --protein-database ${index} \
      --input ${pair_id}.fastq  \
      --output .
    rm ${pair_id}.fastq
    cat ${pair_id}_genefamilies.tsv | grep -v "#" > ${pair_id}_uniref90.tsv
    echo "Uniref90 alignment ended"
    
    # 3. Regroup to the different functional classification types
    ## 3.1 MetaCyc Pathways
    cat ${pair_id}_pathabundance.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | grep -v "UNINTEGRATED" | sed 's/^/MetaCycPath__/' | sed 's/ /_/g' | awk -F '[:\t]' '{print \$1 "\t" \$3}' > ${pair_id}.MetaCyc_Pathways
    echo "MetaCyc Pathways Ended"
    ## 3.2 MetaCyc Reactions
    humann_regroup_table --input ${pair_id}_uniref90.tsv -g uniref90_rxn --output ${pair_id}_metacyc_rxn_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_metacyc_rxn_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/MetaCycRxn__/' > ${pair_id}.MetaCyc_Reactions
    rm ${pair_id}_metacyc_rxn_tmp.tsv
    echo "MetaCyc Reactions Ended"
    ## 3.3 GO
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_go_uniref90.txt.gz --output ${pair_id}_GO_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_GO_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/GO__/' > ${pair_id}.GO
    rm ${pair_id}_GO_tmp.tsv
    echo "GO Ended"
    ## 3.4 KO
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_ko_uniref90.txt.gz --output ${pair_id}_KO_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_KO_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/KO__/' > ${pair_id}.KO
    rm ${pair_id}_KO_tmp.tsv
    echo "KO Ended"
    ## 3.5 EC
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_level4ec_uniref90.txt.gz --output ${pair_id}_EC_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_EC_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/EC__/' > ${pair_id}.EC
    rm ${pair_id}_EC_tmp.tsv
    echo "EC Ended"
    ## 3.6 Pfam
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_pfam_uniref90.txt.gz --output ${pair_id}_Pfam_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_Pfam_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/Pfam__/' > ${pair_id}.Pfam
    rm ${pair_id}_Pfam_tmp.tsv
    echo "Pfam Ended"
    ## 3.7 EggNOG
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_eggnog_uniref90.txt.gz --output ${pair_id}_EggNOG_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_EggNOG_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/EggNOG__/' > ${pair_id}.EggNOG
    rm ${pair_id}_EggNOG_tmp.tsv
    echo "EggNOG Ended"
    ## 3.8 CAZy
    humann_regroup_table --input ${pair_id}_uniref90.tsv -c ${utility}/map_cazy_uniref90.txt.gz --output ${pair_id}_CAZy_tmp.tsv --ungrouped N --protected N
    cat ${pair_id}_CAZy_tmp.tsv | grep -v "#" | grep -v "|" | grep -v "UNMAPPED" | sed 's/^/CAZy__/' > ${pair_id}.CAZy
    rm ${pair_id}_CAZy_tmp.tsv
    echo "CAZy Ended"
    """

}

// 3. Final counts table and reports

// 3.1 Parsing read counts with adjusted counts (Unique + Multiple )
process parse_kaiju {

	  tag ("${id}")
    container 'biocorecrg/almalinux-perlbrew-pyenv3'
    memory = '5G'
    cpus 1
    time '12h'
    
    input:
      tuple val(id), path(kaiju_res)
      val(prof)
	
    output:
	    tuple val(id), path("${id}.kai"), emit: kai
	    tuple val(id), path("${id}.kstats"), emit: kstats
	
    script:
    """
		parseKaijuOut.py -n ${id} -s ${id}.kstats -i ${kaiju_res} -o ${id}.kai -p ${prof}
    """

}

process parse_kaiju_norm {

	  tag ("${id}")
    container 'biocorecrg/almalinux-perlbrew-pyenv3'
    memory = '5G'
    cpus 1
    time '6h'
    
    input:
      tuple val(id), path(kaiju_res), path(genome_stats)
      val(prof)
	
    output:
	    tuple val(id), path("${id}.kai"), emit: kai
	    tuple val(id), path("${id}.kstats"), emit: kstats
	
    script:
    """
		parseKaijuOut.py -n ${id} -g ${genome_stats} -s ${id}.kstats -i ${kaiju_res} -o ${id}.kai -p ${prof}
    """
}

// 3.2 Make final counts table
process makeKaiTable {
  
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') } 
	container 'biocorecrg/almalinux-perlbrew-pyenv3'
	memory = '5G'
	cpus 1
  time '6h'
    
	  input:
      path(kai_files)
      val(norm)
	
    output:
  	  path("out.table"), emit: out
  	  path("norm.table"), emit: norm, optional: true
	
    script:
    def pars  = ""
    if (norm == "YES") {
    	pars = "-n norm.table"
    }
    def inputlist = kai_files.join(',')
    """
	  makeTable.py -f ${inputlist} -o out.table ${pars} -k
    """

}

// 3.3 Make table with mapping statistics
process makeKaijuTable {
	
    container 'biocorecrg/almalinux-perlbrew-pyenv3'
    memory = '5G'
    cpus 1
    time '6h'
    
    input:
      path(kaiju_files)
	
    output:
	    path("kaiju_mqc.txt")
	
    script:
    """
   	echo '# id: kaiju
    # plot_type: table
    # section_name: Kaiju alignment statistics
    # description: Number of reads aligned at least once to a genome
    Sample	reads	aligned' > kaiju_mqc.txt;
    cat ${kaiju_files} | grep -v "aligned" >> kaiju_mqc.txt
    """

}

// 3.4 Make QC Report
process get_stats {
    container 'biocorecrg/almalinux-perlbrew-pyenv3'
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }
    memory = '10G'
    cpus 1
    time '6h'
    
    input:
      path(multiqc_data)
	
    output:
	    path("statistics.txt")
	
    script:
    """
		parseMultiQC.py -i ${multiqc_data} -o statistics.txt
    """
}

// 3.5 Make EC count tables
process makeECTables {
  
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') } 
	container 'biocorecrg/almalinux-perlbrew-pyenv3'
	memory = '1G'
	cpus 1
  time '6h'
    
	  input:
      path(out_table)
	
    output:
  	  path("EC*")
	
    script:
    """
    # Step 1
    awk -F'\\t' 'NR > 1 {
        count = split(\$1, parts, /\\./);
        if (count == 4) print >> "EC4.tmp";
        else if (count == 3) print >> "EC3.tmp";
        else if (count == 2) print >> "EC2.tmp";
        else if (count == 1) print >> "EC1.tmp";
    }' ${out_table}
    
    # Step 2: Collapse EC3
    awk -F'\\t' '{sub(/\\.[^.]*\$/, "", \$1); print \$0 > "EC3.tmp2"}' EC4.tmp
    cat EC3.tmp EC3.tmp2 > EC3.tmp3
    awk '{for (i=2; i<=NF; i++) arr[\$1,i]+=\$i; if (!seen[\$1]++) order[++count]=\$1} END {for (i=1; i<=count; i++) { key=order[i]; printf "%s", key; for (j=2; j<=NF; j++) printf "\\t%s", arr[key, j]; print ""}}' EC3.tmp3 > EC3.tmp4

    # Step 3: Collapse EC2
    awk -F'\\t' '{sub(/\\.[^.]*\$/, "", \$1); print \$0 > "EC2.tmp2"}' EC3.tmp4
    cat EC2.tmp EC2.tmp2 > EC2.tmp3
    awk '{for (i=2; i<=NF; i++) arr[\$1,i]+=\$i; if (!seen[\$1]++) order[++count]=\$1} END {for (i=1; i<=count; i++) { key=order[i]; printf "%s", key; for (j=2; j<=NF; j++) printf "\\t%s", arr[key, j]; print ""}}' EC2.tmp3 > EC2.tmp4
    
    # Step 4: Collapse EC1
    awk -F'\\t' '{sub(/\\.[^.]*\$/, "", \$1); print \$0 > "EC1.tmp2"}' EC2.tmp4
    cat EC1.tmp EC1.tmp2 > EC1.tmp3
    awk '{for (i=2; i<=NF; i++) arr[\$1,i]+=\$i; if (!seen[\$1]++) order[++count]=\$1} END {for (i=1; i<=count; i++) { key=order[i]; printf "%s", key; for (j=2; j<=NF; j++) printf "\\t%s", arr[key, j]; print ""}}' EC1.tmp3 > EC1.tmp4
    
    # Step 5: Final files
    cat ${out_table} | head -1 > EC4.out
    cat ${out_table} | head -1 > EC3.out
    cat ${out_table} | head -1 > EC2.out
    cat ${out_table} | head -1 > EC1.out
    cat EC4.tmp | sed 's/^/EC4__/' >> EC4.out
    cat EC3.tmp4 | sed 's/^/EC3__/' >> EC3.out
    cat EC2.tmp4 | sed 's/^/EC2__/' >> EC2.out
    cat EC1.tmp4 | sed 's/^/EC1__/' >> EC1.out
    
    rm *tmp*
    """

}