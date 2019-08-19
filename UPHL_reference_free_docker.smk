import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.0.2019.08.19")

base_directory=workflow.basedir + "/URF_scripts"
output_directory=os.getcwd()

SAMPLE, MIDDLE, EXTENSION = glob_wildcards('Sequencing_reads/Raw/{sample, [^_]+}_{middle}.f{extension}')
DATABASE = [ 'ncbi', 'serotypefinder', 'vfdb' ]

rule all:
    input:
        # copying files over
        expand("Sequencing_reads/Raw/{sample}_{middle}.f{extension}", zip, sample=SAMPLE, middle=MIDDLE, extension=EXTENSION),
        # running seqyclean
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
        # running FastQC
        "fastqc/fastqc.complete",
        # running shovill
        expand("shovill_result/{sample}/contigs.fa", sample=SAMPLE),
        expand("ALL_assembled/{sample}_contigs.fa", sample=SAMPLE),
        # mash results
        expand("mash/{sample}_mashdist.txt", sample=SAMPLE),
        "mash/mash_results.txt",
        # prokka results
        expand("Prokka/{sample}/{sample}.gff", sample=SAMPLE),
        expand("ALL_gff/{sample}.gff", sample=SAMPLE),
        # quast results
        expand("quast/{sample}/report.tsv", sample=SAMPLE),
        # seqsero results
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE),
        "SeqSero/Seqsero_serotype_results.txt",
        # cg-pipeline results
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE,raw_or_clean=['raw', 'clean']),
        "cg-pipeline/cg-pipeline-summary.txt",
        # abricate results
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results/{database}/{database}.summary.csv", database=DATABASE),
        # blobtools results - which requires bwa, samtools, and blast before blobtools can be run
        expand("shovill_result/{sample}/contigs.fa.sa", sample=SAMPLE),
        expand("bwa/{sample}.sorted.bam", sample=SAMPLE),
        expand("bwa/{sample}.sorted.bam.bai", sample=SAMPLE),
        expand("blast/{sample}.tsv", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.json", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.table.txt", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png", sample=SAMPLE),
        "blobtools/blobtools_results.txt",
        # file summary
        "results_for_multiqc/File_heatmap.csv",
        "run_results.txt"
    singularity:
        "docker://ewels/multiqc:1.7"
    params:
        output_directory=output_directory,
        base_directory=base_directory
    shell:
        "mkdir -p results_for_multiqc ; "
        "mkdir -p logs/all ; "
        "date >> logs/all/all.log ; " # time stamp
        "multiqc --version >> logs/all/all.log ; "
        "wget -nc https://raw.githubusercontent.com/StaPH-B/UPHL/master/URF_scripts/multiqc_config_URF_snakemake_docker.yaml -O multiqc_config.yaml 2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "multiqc -f --outdir {params.output_directory}/logs --cl_config \"prokka_fn_snames: True\" {params.output_directory}/results_for_multiqc "#2>> logs/all/all.err | tee -a logs/all/all.log || true ; "

def get_read1(wildcards):
    read1=glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "*_R1_001.fastq.gz") + glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "_1.fastq")
    return(''.join(read1))

def get_read2(wildcards):
    read2=glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "*_R2_001.fastq.gz") + glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "_2.fastq")
    return(''.join(read2))

rule seqyclean:
    input:
        read1=get_read1,
        read2=get_read2
    output:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq",
        se="Sequencing_reads/QCed/{sample}_clean_SE.fastq",
        sstxt="Sequencing_reads/QCed/{sample}_clean_SummaryStatistics.txt",
        sstsv="Sequencing_reads/QCed/{sample}_clean_SummaryStatistics.tsv",
        log=temp("logs/seqyclean/{sample}")
    threads:
        1
    singularity:
        "docker://staphb/seqyclean:1.10.09"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"seqyclean version: $(seqyclean -h | grep Version)\" >> {output.log}.log ; " # log version
        "seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean "
        "2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output}"

rule fastqc:
    input:
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
    output:
        file="fastqc/fastqc.complete",
        log=temp("logs/fastqc/fastqc")
    threads:
        1
    singularity:
        "docker://dukegcb/fastqc:0.11.4"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "fastqc --version >> {output.log}.log ; " # log version
        "fastqc --outdir fastqc --threads {threads} Sequencing_reads/*/*.fastq* 2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output}"

rule shovill:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    threads:
        48
    output:
        file="shovill_result/{sample}/contigs.fa",
        final="ALL_assembled/{sample}_contigs.fa",
        log=temp("logs/shovill/{sample}")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "shovill --version >> {output.log}.log ; " # logging shovill version
        "RAM=$(free -m --giga | grep \"Mem:\" | awk '{{ print ($2*0.8) }}' | cut -f 1 -d \".\") ; " # getting available RAM
        "echo \"Using $RAM RAM and {threads} cpu for shovill\" >> {output.log}.log ; " # logging the amount of RAM
        "shovill --cpu {threads} --ram $RAM --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force 2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}" # Duplicating files

rule mash_sketch:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        file="mash/{sample}.msh",
        log=temp("logs/mash/{sample}_sketch")
    threads:
        1
    singularity:
        "docker://staphb/mash:2.1"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"mash version: $(mash --version)\" >> {output.log}.log ; " # logging version
        "cat {input.read1} {input.read2} | "
        "mash sketch -m 2 -o mash/{wildcards.sample} - 2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output}"

rule mash_dist:
    input:
        rules.mash_sketch.output.file
    output:
        file="mash/{sample}_mashdist.txt",
        log=temp("logs/mash/{sample}_dist")
    threads:
        1
    singularity:
        "docker://staphb/mash:2.1"
    shell:
        "date 2>> {output.log}.err | tee -a {output.log}.log ; " # time stamp
        "echo \"mash version: $(mash --version)\" >> {output.log}.log ; " # logging version
        "mash dist -p {threads} -v 0 /db/RefSeqSketchesDefaults.msh {input} | sort -gk3 > {output.file} "
        "2>> {output.log}.err || true ; touch {output} "

rule mash_multiqc:
    input:
        expand("mash/{sample}_mashdist.txt", sample=SAMPLE)
    output:
        file="mash/mash_results.txt",
        log=temp("logs/mash/multiqc")
    threads:
        1
    shell:
        "date >> {output.log}.log ; " # time stamp
        "organisms=($(cat mash/*_mashdist.txt | awk '{{ if ( $4 == 0 ) print $1 }}' | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" | sort | uniq -c | sort -rhk 1,1 | awk '{{ print $2 }}' )) ; "
        "echo \"The organisms found in this run: ${{organisms[@]}}\" >> {output.log}.log ; "
        "header=\"Sample\" ; "
        """
        for organism in ${{organisms[@]}}
        do
            header=$(echo \"$header\\t$organism\" )
        done
        """
        "echo -e \"$header\" > mash/mash_results.txt ; "
        "mash_results=($(ls mash/*_mashdist.txt | sed 's!.*/!!' | cut -d \"_\" -f 1 | cut -d '.' -f 1 | sort | uniq )) ; "
        """
        for mash_result in ${{mash_results[@]}}
        do
            mash_line=$mash_result
            for organism in ${{organisms[@]}}
            do
                if [ -z \"$(grep $organism mash/$mash_result*_mashdist.txt | head -n 2 )\" ]
                then
                    number=\"0\"
                else
                    number=$(grep $organism mash/$mash_result*_mashdist.txt | awk '{{ if ( $4 == 0 ) print $1 }}' | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" | sort | uniq -c | grep $organism | awk '{{ print $1 }}' )
                fi
                mash_line=$(echo \"$mash_line\t$number\" )
            done
            echo -e \"$mash_line\" >> mash/mash_results.txt
        done
        """
        "touch {output} "

rule prokka:
    input:
        contig_file=rules.shovill.output.file,
        mash_file=rules.mash_dist.output.file
    threads:
        48
    output:
        file="Prokka/{sample}/{sample}.gff",
        final="ALL_gff/{sample}.gff",
        log=temp("logs/prokka/{sample}")
    singularity:
        "docker://staphb/prokka:1.14.0"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "prokka -v >> {output.log}.log ; " # logging version
        "mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" | sed 's/_/ /g' )) || mash_result=('none', 'none') ; "
        "prokka --cpu {threads} --compliant --centre --URF --mincontiglen 500 --outdir Prokka/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample} --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} 2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}" # duplicating files

rule quast:
    input:
        rules.shovill.output.final
    output:
        file="quast/{sample}/report.tsv",
        log=temp("logs/quast/{sample}")
    threads:
        1
    singularity:
        "docker://staphb/quast:5.0.2"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "quast.py --version >> {output.log}.log ; " # logging version
        "quast.py {input} --output-dir quast/{wildcards.sample} --threads {threads} "
        "2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output} ; "

rule CG_pipeline_shuffle_raw:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        file="Sequencing_reads/shuffled/{sample}_raw_shuffled.fastq.gz",
        log=temp("logs/cg_pipeline/{sample}_shuffle_raw")
    threads:
        1
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log}.log ; " # time stamp, no version
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.log}.err "
        "|| true ; touch {output}"

rule CG_pipeline_shuffle_clean:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        file="Sequencing_reads/shuffled/{sample}_clean_shuffled.fastq.gz",
        log=temp("logs/cg_pipeline/{sample}_shuffle_clean")
    threads:
        1
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log}.log ; " # time stamp, no version
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.log}.err "
        "|| true ; touch {output} "

rule CG_pipeline:
    input:
        shuffled_fastq="Sequencing_reads/shuffled/{sample}_{raw_or_clean}_shuffled.fastq.gz",
        mash_error=rules.mash_sketch.output.log,
        mash_file=rules.mash_dist.output.file
    output:
        file="cg-pipeline/{sample}.{raw_or_clean}.out.txt",
        log=temp("logs/cg_pipeline/{sample}_{raw_or_clean}")
    threads:
        48
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log}.log ; " # time stamp, no version
        "wget -nc https://raw.githubusercontent.com/StaPH-B/UPHL/master/URF_scripts/genome_sizes.txt -O logs/genome_sizes.txt 2>> {output.log}.err || true ; "
        "mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" )) || true ; "
        "genome_length=$(grep $mash_result logs/genome_sizes.txt | grep -v \"#\" | head -n 1 | cut -f 2 -d \":\" | awk '{{ print $0 \"e+06\" }}' ) || genome_length=$(grep 'Estimated genome size:' {input.mash_error}.err | cut -f 4 -d \" \" ) || genome_length=\"0\" ; "
        "echo \"The genome length for {wildcards.sample} is $genome_length\" >> {output.log}.log ; "
        "run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e $genome_length 2>> {output.log}.err > {output.file} "
        "|| true ; touch {output}"

rule CG_pipeline_multiqc:
    input:
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE, raw_or_clean=['raw', 'clean'])
    output:
        file="cg-pipeline/cg-pipeline-summary.txt",
        log=temp("logs/cg_pipeline/multiqc")
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory=output_directory
    shell:
        "date >> {output.log}.log ; " # time stamp, no version
        "grep \"avgReadLength\" cg-pipeline/*.out.txt | sort | uniq | head -n 1 | cut -f 2- -d ':' > cg-pipeline/cg-pipeline-summary.txt || true ; "
        "grep -v \"avgReadLength\" cg-pipeline/*.out.txt | cut -f 2- -d ':' | sort | uniq >> cg-pipeline/cg-pipeline-summary.txt || true ; "
        "touch {output} "

rule seqsero:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read2
    output:
        file="SeqSero/{sample}/Seqsero_result.txt",
        final="SeqSero/{sample}.Seqsero_result.txt",
        log=temp("logs/seqsero/{sample}")
    threads:
        1
    singularity:
        "docker://staphb/seqsero:1.0.1"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} 2>> {output.log}.err >> {output.log}.log "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}"

rule seqsero_multiqc:
    input:
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE)
    output:
        file="SeqSero/Seqsero_serotype_results.txt",
        log=temp("logs/seqsero/multiqc")
    params:
        base_directory=base_directory,
        output_directory=output_directory
    threads:
        1
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo -e \"Sample\tInput_files\tO_antigen_prediction\tH1_antigen_prediction(fliC)\tH2_antigen_prediction(fljB)\tPredicted_antigenic_profile\tPredicted_serotype(s)\" > SeqSero/Seqsero_serotype_results_all.txt ; "
        "RESULTS=$(ls SeqSero/*/Seqsero_result.txt) ; "
        """
        for result in ${{RESULTS[@]}}
        do
          SAMPLE=$(head -n 1 $result | awk '{{ print $3 }}' | awk -F \"_\" '{{ print $1 }}' )
          seqsero_inputfile=$(grep \"Input files\"                 $result | cut -f 2 | tr ' ' '_' )
          seqsero_Oantipred=$(grep \"O antigen prediction\"        $result | cut -f 2 | tr ' ' '_' )
          seqsero_H1antpred=$(grep \"H1 antigen prediction(fliC)\" $result | cut -f 2 | tr ' ' '_' )
          seqsero_H2antpred=$(grep \"H2 antigen prediction(fljB)\" $result | cut -f 2 | tr ' ' '_' )
          seqsero_antigenic=$(grep \"Predicted antigenic profile\" $result | cut -f 2 | tr ' ' '_' )
          seqsero_serotypes=$(grep \"Predicted serotype(s)\"       $result | cut -f 2 | tr ' ' '_' )
          echo -e \"$SAMPLE\t$seqsero_inputfile\t$seqsero_Oantipred\t$seqsero_H1antpred\t$seqsero_H2antpred\t$seqsero_antigenic\t$seqsero_serotypes\" >> SeqSero/Seqsero_serotype_results_all.txt || true
        done
        """
        "grep -v \"O--\" SeqSero/Seqsero_serotype_results_all.txt > SeqSero/Seqsero_serotype_results.txt || true ; "
        "touch {output}"

rule abricate:
    input:
        rules.shovill.output.file
    output:
        file="abricate_results/{database}/{database}.{sample}.out.tab",
        log=temp("logs/abricate/{sample}.{database}")
    threads:
        5
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "abricate --version >> {output.log}.log ; " # version of abricate
        "abricate --list >> {output.log}.log ; " # date of databases
        "abricate --db {wildcards.database} --threads {threads} {input} > {output.file} 2>> {output.log}.err "
        "|| true ; touch {output}"

rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
    output:
        file="abricate_results/{database}/{database}.summary.txt",
        log=temp("logs/abricate/{database}_summary")
    threads:
        1
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "abricate --version >> {output.log}.log ; " # version of abricate
        "abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output.file} 2>> {output.log}.err "
        "|| true ; touch {output}"

rule abricate_multiqc:
    input:
        rules.abricate_summary.output.file
    output:
        file="abricate_results/{database}/{database}.summary.csv",
        log=temp("logs/abricate/{database}_multiqc")
    threads:
        1
    shell:
        "date 2>> {output.log}.err | tee -a {output.log}.log ; "
        "cat {input} | "
        "sed 's/#//g' | sed 's/.tab//g' | sed \"s/{wildcards.database}.//g\" | "
        "awk '{{ sub(\"^.*/\", \"\", $1); print}}' | "
        "awk '{{ for (i=1;i<=NF;i++) if ($i ~ \";\" )gsub(\";.*$\",\"\",$i)g ; else continue}}{{print $0}}' | "
        "awk '{{ $2=\"\" ; print $0 }}' | sed 's/\\t/,/g' | sed 's/ /,/g' | "
        "sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' "
        "> {output.file} 2>> {output.log}.err "
        "|| true ; touch {output} "

rule bwa_index:
    input:
        rules.shovill.output.file
    output:
        index="shovill_result/{sample}/contigs.fa.sa",
        log=temp("logs/bwa/{sample}_index")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"bwa $(bwa 2>&1 | grep Version )\" >> {output.log}.log ; " # version of bwa
        "bwa index {input} 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule blastn:
    input:
        rules.shovill.output.file
    output:
        tsv="blast/{sample}.tsv",
        log=temp("logs/blastn/{sample}")
    threads:
        10
    singularity:
        "docker://ncbi/blast:2.9.0"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "blastn -version >> {output.log}.log ; " # version of blastn
        "echo \"The blastdb location is $BLASTDB\" >> {output.log}.log ; "
        "blastn -query {input} -out {output.tsv} -num_threads {threads} -db /blast/blastdb/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 2>> {output.log}.log | tee -a {output.log}.err || true ; "
        "touch {output}"

rule bwa:
    input:
        contig=rules.shovill.output.file,
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2,
        index=rules.bwa_index.output.index,
        test=rules.blastn.output.tsv
    threads:
        48
    output:
        bam="bwa/{sample}.sorted.bam",
        bai="bwa/{sample}.sorted.bam.bai",
        log=temp("logs/bwa/{sample}")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"bwa $(bwa 2>&1 | grep Version )\" >> {output.log}.log ; " # version of bwa
        "samtools --version >> {output.log}.log ; " # version of samtools
        "if [ -s \"{input.test}\" ] ; then bwa mem -t {threads} {input.contig} {input.read1} {input.read2} 2>> {output.log}.err | samtools sort -o {output.bam} 2>> {output.log}.err > {output.bam} || true ; fi ; " # the if statement is for those without a blast nt database
        "samtools index {output.bam} 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule blobtools_create:
    input:
        contig=rules.shovill.output.file,
        blast=rules.blastn.output.tsv,
        bam=rules.bwa.output.bam
    output:
        cov="blobtools/{sample}.{sample}.sorted.bam.cov",
        json="blobtools/{sample}.blobDB.json",
        log=temp("logs/blobtools/{sample}_create")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"blobtools version $(blobtools -v)\" >> {output.log}.log ; " # version of blobtools
        "blobtools create -o blobtools/{wildcards.sample} -i {input.contig} -b {input.bam} -t {input.blast} 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule blobtools_view:
    input:
        rules.blobtools_create.output.json,
    output:
        txt="blobtools/{sample}.blobDB.table.txt",
        log=temp("logs/blobtools/{sample}_view")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"blobtools version $(blobtools -v)\" >> {output.log}.log ; " # version of blobtools
        "blobtools view -i {input} -o blobtools/ 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule blobtools_plot:
    input:
        table=rules.blobtools_view.output.txt,
        json=rules.blobtools_create.output.json
    output:
        png="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png",
        covpng="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png",
        txt="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt",
        log=temp("logs/blobtools/{sample}_plot")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "echo \"blobtools version $(blobtools -v)\" >> {output.log}.log ; " # version of blobtools
        "blobtools plot -i {input.json} -o blobtools/ -r species --format png 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule blobtools_multiqc:
    input:
        expand("blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt", sample=SAMPLE)
    output:
        file="blobtools/blobtools_results.txt",
        mapping="blobtools/blobtools_mapping.txt",
        log=temp("logs/blobtools/multiqc")
    threads:
        1
    shell:
        "date >> {output.log}.log ; " # time stamp
        "organisms=($(cut -f 1 blobtools/*blobplot.stats.txt | grep -v \"all\" | grep -v ^\"#\" | tr ' ' '_' | sort | uniq -c | sort -rhk 1,1 | awk '{{ print $2 }}' )) ; "
        "echo \"The organisms found in this run: ${{organisms[@]}}\" >> {output.log}.log ; "
        "header=\"Sample\" ; "
        """
        for organism in ${{organisms[@]}}
        do
            header=$(echo \"$header\\t$organism\" )
        done
        """
        "echo -e \"$header\" > {output.file} ; "
        "echo -e \"Sample\tmapped_reads\" > {output.mapping} ; "
        "blobtools_results=($(ls blobtools/*blobplot.stats.txt | sed 's!.*/!!' | cut -d \"_\" -f 1 | cut -d '.' -f 1 | sort | uniq )) ; "
        """
        for blobtools_result in ${{blobtools_results[@]}}
        do
            blobtools_line=$blobtools_result
            for organism in ${{organisms[@]}}
            do
                if [ -z \"$(cut -f 1 blobtools/$blobtools_result*blobplot.stats.txt | tr ' ' '_' | grep $organism | head -n 1 )\" ]
                then
                    number=\"0\"
                else
                    number=$(cat blobtools/$blobtools_result*blobplot.stats.txt | tr ' ' '_' | grep $organism | cut -f 13 | head -n 1 | sed 's/%//g' )
                fi
                blobtools_line=$(echo \"$blobtools_line\t$number\" )
            done
            blobtools_mapping=$(cat blobtools/$blobtools_result*blobplot.stats.txt | tr ' ' '_' | grep all | cut -f 13 | head -n 1 | sed 's/%//g')
            echo -e \"$blobtools_result\t$blobtools_mapping\" >> blobtools/blobtools_mapping.txt
            echo -e \"$blobtools_line\" >> blobtools/blobtools_results.txt
        done
        """
        "touch {output} "

rule multiqc_prep:
    input:
        # copying files over
        expand("Sequencing_reads/Raw/{sample}_{middle}.f{extension}", zip, sample=SAMPLE, middle=MIDDLE, extension=EXTENSION),
        # running seqyclean
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
        # running FastQC
        "fastqc/fastqc.complete",
        # running shovill
        expand("shovill_result/{sample}/contigs.fa", sample=SAMPLE),
        expand("ALL_assembled/{sample}_contigs.fa", sample=SAMPLE),
        # mash results
        expand("mash/{sample}_mashdist.txt", sample=SAMPLE),
        "mash/mash_results.txt",
        # prokka results
        expand("Prokka/{sample}/{sample}.gff", sample=SAMPLE),
        expand("ALL_gff/{sample}.gff", sample=SAMPLE),
        # quast results
        expand("quast/{sample}/report.tsv", sample=SAMPLE),
        # seqsero results
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE),
        "SeqSero/Seqsero_serotype_results.txt",
        # cg-pipeline results
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE,raw_or_clean=['raw', 'clean']),
        "cg-pipeline/cg-pipeline-summary.txt",
        # abricate results
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results/{database}/{database}.summary.csv", database=DATABASE),
        # blobtools results
        expand("blobtools/{sample}.blobDB.json", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.table.txt", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png", sample=SAMPLE),
        "blobtools/blobtools_mapping.txt",
        "blobtools/blobtools_results.txt",
    output:
        log=temp("logs/all/all"),
        file="results_for_multiqc/File_heatmap.csv",
        final="results_for_multiqc/run_results_summary.txt",
        run_results="run_results.txt"

    params:
        output_directory=output_directory,
        base_directory=base_directory
    shell:
        "date 2>> {output.log}.err | tee -a {output.log}.log ; " # time stamp
        # getting the results in the right places
        "{params.base_directory}/check_multiqc_docker.sh                      {params.output_directory}                       2>> {output.log}.err | tee -a {output.log}.log ; "
        "ln -s {params.output_directory}/Prokka*/*/*txt                       {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/SeqSero/Seqsero_serotype_results.txt {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/mash/mash_results.txt                {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/blobtools/blobtools_results.txt      {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/blobtools/blobtools_mapping.txt      {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/cg-pipeline/cg-pipeline-summary.txt  {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/abricate_results/*/*.summary.csv     {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/quast                                {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/logs/File_heatmap.csv                {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/logs/raw_clean_coverage.txt          {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/logs/raw_clean_scatter.csv           {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
#        "ln -s {params.output_directory}/run_file_summary.txt                 {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "ln -s {params.output_directory}/fastqc                               {params.output_directory}/results_for_multiqc/. 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        # formatting for multiqc
        "cat {output.run_results} | sed 's/simple_mash_result/A.simple_mash_result/g' | "
        "sed 's/simple_seqsero_result/B.simple_seqsero_result/g' | "
        "sed 's/abricate_serotype_O/C.abricate_serotype_O/g' | "
        "sed 's/abricate_serotype_H/D.abricate_serotype_H/g' | "
        "sed 's/fastqc_raw_reads_2/E.fastqc_raw_reads_2/g' | "
        "sed 's/fastqc_clean_reads_PE2/F.fastqc_clean_reads_PE2/g' | "
        "sed 's/cg_raw_coverage/G.cg_raw_coverage/g' | "
        "sed 's/cg_cln_coverage/H.cg_cln_coverage/g' | "
        "sed 's/ncbi/J.ncbi_antibiotic_resistence_genes/g' | "
        "sed 's/stxeae_result/I.stx_and_eae_virulence_factor_result/g' | "
        "sed 's/blobtools/K.blobtools/g' > results_for_multiqc/run_results_summary.txt || true ; "
        "touch {output}"
