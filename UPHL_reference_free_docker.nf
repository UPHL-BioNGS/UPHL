#!/usr/bin/env nextflow

println("UPHL Reference Free Workflow v.20200605")

//# nextflow run ~/sandbox/UPHL/UPHL_reference_free_docker.nf
// nextflow run ~/sandbox/UPHL/UPHL_reference_free_docker.nf -c ~/sandbox/UPHL/URF_scripts/singularity.nextflow.config

params.outdir = workflow.launchDir
params.genome_sizes = workflow.projectDir + "/URF_scripts/genome_sizes.txt"
params.center = 'UPHL'
params.seqyclean_contaminant_file = "/Adapters_plus_PhiX_174.fasta"
params.mash_refseq_sketch = "/db/RefSeqSketchesDefaults.msh"
params.blast_db = "/home/Bioinformatics/Data/blastdb/"
params.blast_db_container = "/blast/blastdb"
params.abricate_db = ["ecoh", "vfdb", "serotypefinder", "ncbi" ]

list_of_known_organisms=params.genome_sizes.readLines()

maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${maxcpus}")
if ( maxcpus < 5 ) {
  medcpus = maxcpus
} else {
  medcpus = 5
}

maxmem = Math.round(Runtime.runtime.totalMemory() / 10241024)
println("The maximum amount of memory used in this workflow is ${maxmem}")

Channel
  .fromFilePairs(["${params.outdir}/Sequencing_reads/Raw/*_R{1,2}_001.fastq.gz", "${params.outdir}/Sequencing_reads/Raw/*_{1,2}.fastq" ], size: 2 )
  .map{ reads -> [reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]] }
  .ifEmpty{ exit 1, println("No paired fastq or fastq.gz files were found at ${params.outdir}/Sequencing_reads/Raw") }
  .into { fastq_reads; fastq_reads2; fastq_reads3 }

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p Sequencing_reads/QCed logs/seqyclean'

  input:
  set val(sample), file(reads) from fastq_reads

  output:
  tuple sample, file("Sequencing_reads/QCed/${sample}_clean_PE{1,2}.fastq") into clean_reads, clean_reads2, clean_reads3, clean_reads4, clean_reads5, clean_reads6
  file("Sequencing_reads/QCed/${sample}_clean_SE.fastq")
  file("Sequencing_reads/QCed/${sample}_clean_SummaryStatistics.{txt,tsv}")
  file("logs/seqyclean/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file

    seqyclean -minlen 25 -qual -c !{params.seqyclean_contaminant_file} -1 !{reads[0]} -2 !{reads[1]} -o Sequencing_reads/QCed/!{sample}_clean 2>> $err_file >> $log_file
  '''
}

fastq_reads2
  .combine(clean_reads2, by: 0)
  .into { raw_clean_reads; raw_clean_reads2 }

process shovill {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus maxcpus
  memory maxmem + ' GB'

  beforeScript 'mkdir -p shovill ALL_assembled logs/shovill'

  input:
  set val(sample), file(reads) from clean_reads

  output:
  tuple sample, file("shovill/${sample}/contigs.fa") into contigs, contigs2, contigs4, contigs5, contigs6, contigs7, contigs8
  file("shovill/${sample}/contigs.gfa")
  file("shovill/${sample}/shovill.{corrections,log}")
  file("shovill/${sample}/spades.fasta")
  file("ALL_assembled/${sample}_contigs.fa") into contigs3
  file("logs/shovill/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/shovill/!{sample}.!{workflow.sessionId}.log
    err_file=logs/shovill/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    shovill --version >> $log_file

    shovill --cpu !{task.cpus} --ram !{task.memory} --outdir shovill/!{sample} --R1 !{reads[0]} --R2 !{reads[1]} --force 2>> $err_file >> $log_file
    cp shovill/!{sample}/contigs.fa ALL_assembled/!{sample}_contigs.fa
  '''
}

process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p fastqc logs/fastqc'

  input:
  set val(sample), file(raw), file(clean) from raw_clean_reads

  output:
  path("fastqc/")
  tuple sample, env(raw_1), env(raw_2), env(clean_1), env(clean_2) into fastqc_results
  file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/fastqc/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastqc/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc --outdir fastqc --threads !{task.cpus} !{sample}*.fastq* 2>> $err_file >> $log_file

    raw_1=$(unzip -p fastqc/!{raw[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=$(unzip -p fastqc/!{raw[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    clean_1=$(unzip -p fastqc/!{clean[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    clean_2=$(unzip -p fastqc/!{clean[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
  '''
}

process mash_sketch {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p mash logs/mash_sketch'

  input:
  set val(sample), file(reads) from clean_reads3

  output:
  tuple sample, file("mash/${sample}.msh") into sketches
  tuple sample, env(mash_genome_size) into mash_genome_size
  file("logs/mash_sketch/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/mash_sketch/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_sketch/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version)" >> $log_file

    cat !{reads[0]} !{reads[1]} | mash sketch -m 2 -o mash/!{sample} - 2>> $err_file >> $log_file
    mash_genome_size=$(grep 'Estimated genome size:' $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' | xargs printf "%.0f")
  '''
}

process mash_dist {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus medcpus

  beforeScript 'mkdir -p mash logs/mash_dist'

  input:
  set val(sample), file(sketch) from sketches

  output:
  tuple sample, file("mash/${sample}_mashdist.txt")
  file("logs/mash_dist/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(mash_simple) into mash_organism, mash_organism2, mash_organism3, mash_organism4
  tuple sample, env(mash_full), env(mash_simple) into mash_results

  shell:
  '''
    log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version)" >> $log_file

    mash dist -p !{task.cpus} -v 0 !{params.mash_refseq_sketch} !{sketch} | sort -gk3 > mash/!{sample}_mashdist.txt 2>> $err_file

    if [ ! -s "mash/!{sample}_mashdist.txt" ]
    then
      echo "No p-value=0 results were found for !{sample}, expanding to the top 100 hits" | tee -a $log_file $err_file > /dev/null
      mash dist -p !{task.cpus} !{params.mash_refseq_sketch} !{sketch} | sort -gk3 | head -n 100 > mash/!{sample}_mashdist.txt 2>> $err_file
    fi

    mash_full=$(cat mash/!{sample}_mashdist.txt | head -n 1 | cut -f 1)
    mash_simple=$(echo $mash_full | cut -f 8 -d "-" | sed 's/^_//g' | cut -f 1,2 -d "_" | cut -f 1 -d "." )
  '''
}

contigs
  .combine(mash_organism, by: 0)
  .set { contig_organism }

process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus maxcpus

  beforeScript "mkdir -p prokka/${sample} ALL_gff logs/prokka"

  input:
  set val(sample), file(contig), val(mash) from contig_organism

  output:
  file("Prokka/${sample}/${sample}.gff")
  file("ALL_gff/${sample}.gff")
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/prokka/!{sample}.!{workflow.sessionId}.log
    err_file=logs/prokka/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    prokka -v >> $log_file

    mash_result=(echo !{mash} | sed 's/_/ /g')

    prokka --cpu !{task.cpus} \
      --compliant \
      --centre !{params.center} \
      --mincontiglen 500 \
      --outdir prokka/!{sample} \
      --locustag locus_tag \
      --prefix !{sample} \
      --genus ${mash_result[0]} \
      --species ${mash_result[1]} \
      --force !{contig} 2>> $err_file >> $log_file

    cp prokka/!{sample}/!{sample}.gff ALL_gff/!{sample}.gff

    exit 1
  '''
}

process quast {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus medcpus

  beforeScript 'mkdir -p quast logs/quast'

  input:
  set val(sample), file(contig) from contigs4

  output:
  file("quast/${sample}/icarus.html")
  path("quast/${sample}/{basic_stats,icarus_viewers}")
  file("quast/${sample}/report.{html,pdf,tex,txt,tsv}")
  file("quast/${sample}/transposed_report.{tex,tsv,txt}")
  file("logs/quast/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/quast/!{sample}.!{workflow.sessionId}.log
    err_file=logs/quast/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    quast.py --version >> $log_file

    quast.py !{contig} --output-dir quast/!{sample} --threads !{task.cpus} 2>> $err_file | tee -a $log_file
  '''
}

process cg_pipeline_shuffle {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p Sequencing_reads/shuffled logs/cg_pipeline_shuffle'

  input:
  set val(sample), file(raw), file(clean) from raw_clean_reads2

  output:
  tuple sample, file("Sequencing_reads/shuffled/${sample}_raw_shuffled.fastq.gz"), file("Sequencing_reads/shuffled/${sample}_clean_shuffled.fastq.gz") into shuffled
  file("logs/cg_pipeline_shuffle/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/cg_pipeline_shuffle/!{sample}.!{workflow.sessionId}.log
    err_file=logs/cg_pipeline_shuffle/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null

    run_assembly_shuffleReads.pl -gz !{raw[0]} !{raw[1]} > Sequencing_reads/shuffled/!{sample}_raw_shuffled.fastq.gz 2>> $err_file
    run_assembly_shuffleReads.pl -gz !{clean[0]} !{clean[1]} > Sequencing_reads/shuffled/!{sample}_clean_shuffled.fastq.gz 2>> $err_file
  '''
}

shuffled
  .combine(mash_genome_size, by:0)
  .combine(mash_organism3, by:0)
  .set { shuffled_mash }

process cg_pipeline {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p cg_pipeline logs/cg_pipeline'

  input:
  set val(sample), file(raw), file(clean), val(genome_size), val(organism) from shuffled_mash

  output:
  file("cg_pipeline/${sample}.{raw,clean}.out.txt")
  tuple sample, env(raw_cov), env(clean_cov) into cg_pipeline_results
  file("logs/cg_pipeline/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  if (list_of_known_organisms.findAll { it.contains(organism) })
    '''
      log_file=logs/cg_pipeline/!{sample}.!{workflow.sessionId}.log
      err_file=logs/cg_pipeline/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null

      organism=$(echo !{organism} | sed 's/ /_/g')
      genome_length=$(grep $organism !{params.genome_sizes} | grep -v "#" | head -n 1 | cut -f 2 -d ":" | awk '{ print $0 "e+06" }' )

      echo "The genome was found for !{organism} in !{params.genome_sizes}" >> $log_file

      run_assembly_readMetrics.pl !{raw} --fast --numcpus !{task.cpus} -e $genome_length 2>> $err_file > cg_pipeline/!{sample}.raw.out.txt
      run_assembly_readMetrics.pl !{clean} --fast --numcpus !{task.cpus} -e $genome_length 2>> $err_file > cg_pipeline/!{sample}.clean.out.txt

      raw_cov=$(tail -n 1 cg_pipeline/!{sample}.raw.out.txt | head -n 1 | awk '{print $9 }')
      clean_cov=$(tail -n 1 cg_pipeline/!{sample}.clean.out.txt | head -n 1 | awk '{print $9 }')
    '''
  else
    '''
      log_file=logs/cg_pipeline/!{sample}.!{workflow.sessionId}.log
      err_file=logs/cg_pipeline/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null

      echo "The genome was not found for !{organism} in !{params.genome_sizes}" >> $log_file

      run_assembly_readMetrics.pl !{raw} --fast --numcpus 1 -e !{genome_size} 2>> $err_file > cg_pipeline/!{sample}.raw.out.txt
      run_assembly_readMetrics.pl !{clean} --fast --numcpus 1 -e !{genome_size} 2>> $err_file > cg_pipeline/!{sample}.clean.out.txt

      raw_cov=$(tail -n 1 cg_pipeline/!{sample}.raw.out.txt | head -n 1 | awk '{print $9 }')
      clean_cov=$(tail -n 1 cg_pipeline/!{sample}.clean.out.txt | head -n 1 | awk '{print $9 }')
    '''
}

clean_reads6
  .combine(mash_organism4, by:0)
  .set { clean_mash }

process seqsero2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus medcpus

  beforeScript 'mkdir -p seqsero2/${sample} logs/seqsero2'

  input:
  set val(sample), file(reads), val(mash) from clean_mash

  output:
  file("seqsero2/${sample}/data_log.txt")
  file("seqsero2/${sample}/SeqSero_result.{tsv,txt}")
  file("seqsero2/${sample}/${sample}_clean_PE1.fastq_H_and_O_and_specific_genes.fasta_mem.fasta")
  file("seqsero2/${sample}/blasted_output.xml")
  file("seqsero2/${sample}/Extracted_antigen_alleles.fasta")
  file("seqsero2/${sample}/SeqSero_log.txt")
  tuple sample, env(antigenic_profile), env(predicted_serotype), env(notes) into seqsero2_results
  file("logs/seqsero2/${sample}.${workflow.sessionId}.log")
  file("logs/seqsero2/${sample}.${workflow.sessionId}.err")

  when:
    mash == 'Salmonella_enterica_'

  shell:
  if (mash == 'Salmonella_enterica_' )
  '''
    rm -rf seqsero2/${sample}/*

    log_file=logs/seqsero2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqsero2/!{sample}.!{workflow.sessionId}.err


    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    SeqSero2_package.py --version >> $log_file
    SeqSero2_package.py --check >> $log_file

    SeqSero2_package.py -p !{task.cpus} -n !{sample} -d seqsero2/!{sample} -t 2 -i !{reads[0]} !{reads[1]}

    antigenic_profile=$(grep "Predicted antigenic profile:"	seqsero2/!{sample}/SeqSero_result.txt | cut -f 2)
    predicted_serotype=$(grep "Predicted serotype:"	seqsero2/!{sample}/SeqSero_result.txt | cut -f 2)
    notes=$(grep "Notes:" seqsero2/!{sample}/SeqSero_result.txt | cut -f 2)
  '''
  else
  '''
    antigenic_profile='not_salmonella'
    predicted_serotype='not_salmonella'
    notes='not_salmonella'
  '''
}

contigs3
  .combine(mash_organism2, by:0)
  .set { contigs_organism }

process abricate {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus medcpus

  beforeScript "mkdir -p abricate/${db} logs/abricate"

  input:
  each db from params.abricate_db
  set val(sample), file(contigs), val(organism) from contigs_organism

  output:
  file("abricate/${db}/${db}.${sample}.out.tab") into abricate
  tuple sample, ( db : env(abricate_result)) into abricate_results
  file("logs/abricate/${sample}.${db}.${workflow.sessionId}.{log,err}")

  shell:
  if( db == 'ecoh'  &&  organism == "Escherichia_coli" )
  '''
    log_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.log
    err_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file
    abricate --list >> $log_file

    abricate --db !{db} --threads !{task.cpus} --minid 80 --mincov 80 !{contigs} > abricate/!{db}/!{db}.!{sample}.out.tab 2>> $err_file

    O_group=($(grep !{sample} abricate/!{db}/!{db}.!{sample}.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "O"  ))
    H_group=($(grep !{sample} abricate/!{db}/!{db}.!{sample}.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "H"  ))
    O=$(echo ${O_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_O" ]; then abricate_serotype_O="none"; fi
    H=$(echo ${H_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_H" ]; then abricate_serotype_H="none"; fi
    echo $O_group $O $H_group $H

    exit 1
  '''
  else if ( db == 'serotypefinder'  &&  organism == "Escherichia_coli" )
  '''
    log_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.log
    err_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file
    abricate --list >> $log_file

    abricate --db !{db} --threads !{task.cpus} --minid 80 --mincov 80 !{contigs} > abricate/!{db}/!{db}.!{sample}.out.tab 2>> $err_file

    O_group=($(grep !{sample} abricate/!{db}/!{db}.!{sample}.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "O"   ))
    H_group=($(grep !{sample} abricate/!{db}/!{db}.!{sample}.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "H" ))
    O=$(echo ${O_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_O" ]; then abricate_serotype_O="none"; fi
    H=$(echo ${H_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_H" ]; then abricate_serotype_H="none"; fi
    echo $O_group $O $H_group $H
        exit 1
  '''
  else
  '''
    log_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.log
    err_file=logs/abricate/!{sample}.!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file
    abricate --list >> $log_file

    abricate --db !{db} --threads !{task.cpus} --minid 80 --mincov 80 !{contigs} > abricate/!{db}/!{db}.!{sample}.out.tab 2>> $err_file

    abricate_result=$(grep !{sample} abricate/!{db}/!{db}.!{sample}.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | tr ' ' '_')
    if [ -z "$abricate_result" ]; then abricate_result="not_found"; fi
    echo $abricate_result
  '''
}

process abricate_summary {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p abricate/summary logs/abricate_summary'

  input:
  each db from params.abricate_db
  file(abricate) from abricate.collect()

  output:
  file("abricate/summary/${db}.abricate_summary.txt")
  file("abricate/${db}.stxeae.all.txt") optional true
  file("abricate/${db}.all.txt")
  file("logs/abricate_summary/summary.${db}.${workflow.sessionId}.{log,err}")

  shell:
  if( db == 'vfdb')
  '''
    log_file=logs/abricate_summary/summary.!{db}.!{workflow.sessionId}.log
    err_file=logs/abricate_summary/summary.!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file

    abricate --summary !{db}*out.tab > abricate/summary/!{db}.abricate_summary.txt 2>> $err_file
    cat !{db}*out.tab > abricate/!{db}.all.txt
    cat abricate/!{db}.all.txt | grep -e 'eae' -e 'stx' > abricate/!{db}.stxeae.all.txt
  '''
  else
  '''
    log_file=logs/abricate_summary/summary.!{db}.!{workflow.sessionId}.log
    err_file=logs/abricate_summary/summary.!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file

    abricate --summary !{db}*out.tab > abricate/summary/!{db}.abricate_summary.txt 2>> $err_file
    cat !{db}*out.tab > abricate/!{db}.all.txt
  '''
}

process bwa_index {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p bwa logs/bwa_index'

  input:
  set val(sample), file(contigs) from contigs6

  output:
  tuple sample, file("bwa/contigs.fa.{amb,ann,bwt,pac,sa}") into index
  file("logs/bwa_index/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/bwa_index/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa_index/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file

    bwa index !{contigs} 2>> $err_file >> $log_file

    mv contigs.fa.* bwa/.
  '''
}

process blastn {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus medcpus

  beforeScript 'mkdir -p blast logs/blastn'

  input:
  set val(sample), file(contigs) from contigs7

  output:
  tuple sample, file("blast/${sample}.tsv") into blast
  file("logs/blastn/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/blastn/!{sample}.!{workflow.sessionId}.log
    err_file=logs/blastn/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    blastn -version >> $log_file
    echo "The blastdb location is $BLASTDB" >> $log_file

    blastn -query !{contigs} \
      -out blast/!{sample}.tsv \
      -num_threads !{task.cpus} \
      -db !{params.blast_db_container}/nt \
      -outfmt '6 qseqid staxids bitscore std' \
      -max_target_seqs 10 \
      -max_hsps 1 \
      -evalue 1e-25 2>> $err_file >> $log_file
   '''
}

index
  .combine(clean_reads5, by:0)
  .combine(contigs8, by:0)
  .set { index_clean_contigs }

process bwa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus maxcpus

  beforeScript 'mkdir -p bwa logs/bwa'

  input:
  set val(sample), file(index), file(reads), file(contigs) from index_clean_contigs

  output:
  tuple sample, file("bwa/${sample}.sorted.bam"), file("bwa/${sample}.sorted.bam.bai") into bams
  file("logs/bwa/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/bwa/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    samtools --version >> {output.log}.log

    bwa mem -t !{task.cpus} !{contigs} !{reads[0]} !{reads[1]} 2>> $err_file | \
      samtools sort -o bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file

    samtools index bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file

    exit 1
  '''
}

blast
  .combine(bams, by:0)
  .combine(contigs2, by:0)
  .set { blastn_bwa_contigs }

process blobtools_create {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p blobtools logs/blobtools_create'

  input:
  set val(sample), file(blast), file(bwa), file(contigs) from blastn_bwa_contigs

  output:
  file("blobtools/${sample}.${sample}.sorted.bam.cov")
  tuple sample, file("blobtools/${sample}.blobDB.json") into create, create2
  file("logs/blobtools_create/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/blobtools_create/!{sample}.!{workflow.sessionId}.log
    err_file=logs/blobtools_create/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "blobtools version $(blobtools -v)" >> $log_file

    blobtools create -o blobtools/!{sample} -i !{contigs} -b !{bam} -t !{blast} 2>> $err_file >> $log_file

    exit 1
  '''
}

process blobtools_view {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p logs/blobtools_view'

  input:
  set val(sample), file(json) from create

  output:
  tuple sample, file("blobtools/${sample}.blobDB.table.txt") into view
  file("logs/blobtools_view/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/blobtools_view/!{sample}.!{workflow.sessionId}.log
    err_file=logs/blobtools_view/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "blobtools version $(blobtools -v)" >> $log_file

    blobtools view -i !{json} -o blobtools/ 2>> $err_file >> $log_file

    exit 1
  '''
}

create2
  .combine(view, by:0)
  .set { create_view }

process blobtools_plot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p logs/blobtools_plot'

  input:
  set val(sample), file(json), file(view) from create_view

  output:
  file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png")
  file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png")
  file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt")
  tuple sample, env(result) into blobtools_results
  file("logs/blobtools_plot/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/blobtools_plot/!{sample}.!{workflow.sessionId}.log
    err_file=logs/blobtools_plot/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "blobtools version $(blobtools -v)" >> $log_file

    blobtools plot -i !{cov} -o blobtools/ -r species --format png 2>> $err_file >> $log_file

    result=$(grep -v ^"#" $out/blobtools/$sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt | grep -v ^"all" | head -n 1 | tr ' ' '_' | cut -f 1,13)

    exit 1
  '''
}

process mlst {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p mlst logs/mlst'

  input:
  tuple val(sample), file(contigs) from contigs5

  output:
  file("mlst/mlst.${sample}.txt")
  tuple sample, env(result) into mlst_results
  file("logs/mlst/mlst.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/mlst/mlst.!{workflow.sessionId}.log
    err_file=logs/mlst/mlst.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    mlst --version >> $log_file

    mlst !{contigs} > mlst/mlst.!{sample}.txt 2>> $err_file

    result=$(cat mlst/mlst.!{sample}.txt | awk '{ print "MLST" $2 ",PubMLST" $3 }' )
  '''
}

fastqc_results
  .combine(mash_results, by:0)
  .combine(cg_pipeline_results, by:0)
  .combine(seqsero2_results, by:0)
  .combine(abricate_results, by:0)
  .combine(blobtools_results, by:0)
  .combine(mlst_results, by:0)
  .set { workflow_results }

workflow_results.view()

// process multiqc {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p multiqc logs/multiqc'
//
//   input:
//
//   output:
//   tuple sample, file("mash/${sample}_mashdist.txt")
//   file("logs/multiqc/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/multiqc/multiqc.!{workflow.sessionId}.log
//     err_file=logs/multiqc/multiqc.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     multiqc --version >> $log_file
//
//     multiqc -f \
//       --outdir ./logs \
//       --cl_config "prokka_fn_snames: True"  \
//       !{params.outdir}/abricate_results/summary \
//       !{params.outdir}/blobtools \
//       !{params.outdir}/cg-pipeline \
//       !{params.outdir}/fastqc \
//       !{params.outdir}/mash \
//       !{params.outdir}/prokka \
//       !{params.outdir}/quast \
//       !{params.outdir}/seqsero2 \
//       !{params.outdir}/Sequencing_reads/QCed/*tsv \
//       2>> $err_file | tee -a $log_file
//
//     exit 1
//   '''
// }
