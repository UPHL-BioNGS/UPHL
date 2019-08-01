import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.2019.08.05")

base_directory=workflow.basedir + "/URF_scripts"
output_directory=os.getcwd()

SAMPLE, MIDDLE, EXTENSION = glob_wildcards('Sequencing_reads/Raw/{sample, [^_]+}_{middle}.f{extension}')
DATABASE = [ 'ncbi', 'serotypefinder' ]

rule all:
    input:
        "logs/final.txt"
    singularity:
        "docker://ewels/multiqc:1.7"
    params:
        output_directory=output_directory,
        base_directory=base_directory
    shell:
        "date >> logs/all/all.log ; " # time stamp
        "multiqc --version 2>> logs/all/all.err | tee -a logs/all/all.log ; "
        "wget https://raw.githubusercontent.com/StaPH-B/UPHL/master/URF_scripts/multiqc_config_URF_snakemake.yaml 2>> logs/all/all.err | tee -a logs/all/all.log ; "
        "mv multiqc_config_URF_snakemake.yaml multiqc_config.yaml 2>> logs/all/all.err | tee -a logs/all/all.log ; "
        "multiqc -f --outdir {params.output_directory}/logs --cl_config \"prokka_fn_snames: True\" {params.output_directory}/results_for_multiqc 2>> logs/all/all.err | tee -a logs/all/all.log || true ; "

rule moving:
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
        expand("quast/{sample}/report.txt", sample=SAMPLE),
        # seqsero results
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE),
        "SeqSero/Seqsero_serotype_results.txt",
        # cg-pipeline results
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE,raw_or_clean=['raw', 'clean']),
        "cg-pipeline/cg-pipeline-summary.txt",
        # abricate results
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("logs/abricate_results/{database}.summary.csv", database=DATABASE),
    output:
        "logs/final.txt"
    params:
        output_directory=output_directory,
        base_directory=base_directory
    shell:
        # getting the results in the right places
        "mkdir -p logs/all ; "
        "mkdir -p results_for_multiqc ; "
        "echo -e \"made it here!\n\" ; "
        "{params.base_directory}/check_multiqc_docker.sh                      {params.output_directory}                            2>> logs/all/all.err | tee -a logs/all/all.log ; "
        "echo -e \"made it here!a\n\" ; "
        "ln -s {params.output_directory}/fastqc                               {params.output_directory}/results_for_multiqc/fastqc 2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!1\n\" ; "
        "ln -s {params.output_directory}/Prokka*/*/*txt                       {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!2\n\" ; "
        "ln -s {params.output_directory}/SeqSero/Seqsero_serotype_results.txt {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!3\n\" ; "
        "ln -s {params.output_directory}/mash/mash_results.txt                {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!4\n\" ; "
        "ln -s {params.output_directory}/cg-pipeline/cg-pipeline-summary.txt  {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!5\n\" ; "
        "ln -s {params.output_directory}/logs/abricate_results/*.summary.csv  {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!6\n\" ; "
        "ln -s {params.output_directory}/quast                                {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!7\n\" ; "
        "ln -s {params.output_directory}/logs/File_heatmap.csv                {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!8\n\" ; "
        "ln -s {params.output_directory}/logs/raw_clean_coverage.txt          {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!9\n\" ; "
        "ln -s {params.output_directory}/logs/raw_clean_scatter.csv           {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        "echo -e \"made it here!10\n\" ; "
        "ln -s {params.output_directory}/run_file_summary.txt                 {params.output_directory}/results_for_multiqc/.      2>> logs/all/all.err | tee -a logs/all/all.log || true ; "
        # formatting for multiqc
        "cat run_results_summary.txt | sed 's/simple_mash_result/A.simple_mash_result/g' | sed 's/simple_seqsero_result/B.simple_seqsero_result/g' | "
        "sed 's/abricate_serotype_O/C.abricate_serotype_O/g' | sed 's/abricate_serotype_H/D.abricate_serotype_H/g' | sed 's/fastqc_raw_reads_2/E.fastqc_raw_reads_2/g' | "
        "sed 's/fastqc_clean_reads_PE2/F.fastqc_clean_reads_PE2/g' | sed 's/cg_raw_coverage/G.cg_raw_coverage/g' | sed 's/cg_cln_coverage/H.cg_cln_coverage/g' | "
        "sed 's/ncbi/J.ncbi_antibiotic_resistence_genes/g' | sed 's/stxeae_result/I.stx_and_eae_virulence_factor_result/g' > results_for_multiqc/run_results_summary.txt || true ; "
        "touch {output}"


def get_read1(wildcards):
    read1=glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "*_R1_001.fastq.gz") + glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "_1.fastq")
    return(''.join(read1))

def get_read2(wildcards):
    read2=glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "*_R2_001.fastq.gz") + glob.glob("Sequencing_reads/Raw/" + wildcards.sample + "_2.fastq")
    return(''.join(read2))

def get_reads(wildcards):
    reads=glob.glob("Sequencing_reads/Raw/*fastq*")
    return(reads)

rule seqyclean:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq",
        se="Sequencing_reads/QCed/{sample}_clean_SE.fastq",
        sstxt="Sequencing_reads/QCed/{sample}_clean_SummaryStatistics.txt",
        sstsv="Sequencing_reads/QCed/{sample}_clean_SummaryStatistics.tsv",
        err="logs/seqyclean/{sample}.err",
        log="logs/seqyclean/{sample}.log"
    threads:
        1
    singularity:
        "docker://staphb/seqyclean:1.10.09"
    shell:
        "date >> {output.log} ; " # time stamp
        "echo \"seqyclean version: $(seqyclean -h | grep Version)\" >> {output.log} ; " # log version
        "seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean "
        "2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output}"

rule fastqc:
    input:
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
    output:
        file="fastqc/fastqc.complete",
        err="logs/fastqc/fastqc.err",
        log="logs/fastqc/fastqc.log"
    threads:
        1
    singularity:
        "docker://dukegcb/fastqc:0.11.4"
    shell:
        "date >> {output.log} ; " # time stamp
        "fastqc --version >> {output.log} ; " # log version
        "fastqc --outdir fastqc --threads {threads} Sequencing_reads/*/*.fastq* 2>> {output.err} | tee -a {output.log} "
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
        log="logs/shovill/{sample}.log",
        err="logs/shovill/{sample}.err"
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        "date >> {output.log} ; " # time stamp
        "shovill --version >> {output.log} ; " # logging shovill version
        "RAM=$(free -m --giga | grep \"Mem:\" | awk '{{ print ($2*0.8) }}' | cut -f 1 -d \".\") ; " # getting available RAM
        "echo \"Using $RAM RAM and {threads} cpu for shovill\" >> {output.log} ; " # logging the amount of RAM
        "shovill --cpu {threads} --ram $RAM --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force "
        "2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}" # Duplicating files

rule mash_sketch:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        file="mash/{sample}.msh",
        log="logs/mash_sketch/{sample}.log",
        err="logs/mash_sketch/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/mash:2.1"
    shell:
        "date >> {output.log} ; " # time stamp
        "echo \"mash version: $(mash --version)\" >> {output.log} ; " # logging version
        "cat {input.read1} {input.read2} | "
        "mash sketch -m 2 -o mash/{wildcards.sample} - 2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output.file}"

rule mash_dist:
    input:
        rules.mash_sketch.output.file
    output:
        file="mash/{sample}_mashdist.txt",
        log="logs/mash_dist/{sample}.log",
        err="logs/mash_dist/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/mash:2.1"
    shell:
        "date 2>> {output.err} | tee -a {output.log} ; " # time stamp
        "echo \"mash version: $(mash --version)\" >> {output.log} ; " # logging version
        "mash dist -p {threads} -v 0 /db/RefSeqSketchesDefaults.msh {input} | sort -gk3 > {output.file} "
        "2>> {output.err} || true ; touch {output} "

rule mash_multiqc:
    input:
        expand("mash/{sample}_mashdist.txt", sample=SAMPLE)
    output:
        file="mash/mash_results.txt",
        log="logs/mash_pipeline_multiqc/log.log",
        err="logs/mash_pipeline_multiqc/log.err"
    threads:
        1
    shell:
        "date >> {output.log} ; " # time stamp
        "organisms=($(cat mash/*_mashdist.txt | awk '{{ if ( $4 == 0 ) print $1 }}' | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" | sort | uniq -c | sort -rhk 1,1 | awk '{{ print $2 }}' )) ; "
        "echo \"The organisms found in this run: ${{organisms[@]}}\" >> {output.log} ; "
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
        log="logs/prokka/{sample}.log",
        err="logs/prokka/{sample}.err"
    singularity:
        "docker://staphb/prokka:1.13"
    shell:
        "date >> {output.log} ; " # time stamp
        "prokka -v >> {output.log} ; " # logging version
        "mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" | sed 's/_/ /g' )) || mash_result=('none', 'none') ; "
        "prokka --cpu {threads} --compliant --centre --URF --mincontiglen 500 --outdir Prokka/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample} --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} 2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}" # duplicating files

rule quast:
    input:
        rules.shovill.output.file
    output:
        file="quast/{sample}/report.txt",
        log="logs/quast/{sample}.log",
        err="logs/quast/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/quast:5.0.2"
    shell:
        "date >> {output.log} ; " # time stamp
        "quast.py --version >> {output.log} ; " # logging version
        "quast.py {input} --output-dir quast/{wildcards.sample} --threads {threads} "
        "2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} "

rule CG_pipeline_shuffle_raw:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        file="Sequencing_reads/shuffled/{sample}_raw_shuffled.fastq.gz",
        log="logs/cg_pipeline_shuffle_raw/{sample}.log",
        err="logs/cg_pipeline_shuffle_raw/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log} ; " # time stamp, no version
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.err} "
        "|| true ; touch {output}"

rule CG_pipeline_shuffle_clean:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        file="Sequencing_reads/shuffled/{sample}_clean_shuffled.fastq.gz",
        log="logs/cg_pipeline_shuffle_clean/{sample}.log",
        err="logs/cg_pipeline_shuffle_clean/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log} ; " # time stamp, no version
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.err} "
        "|| true ; touch {output} "

rule CG_pipeline:
    input:
        shuffled_fastq="Sequencing_reads/shuffled/{sample}_{raw_or_clean}_shuffled.fastq.gz",
        mash_error=rules.mash_sketch.output.err,
        mash_file=rules.mash_dist.output.file
    output:
        file="cg-pipeline/{sample}.{raw_or_clean}.out.txt",
        log="logs/cg_pipeline/{sample}.{raw_or_clean}.log",
        err="logs/cg_pipeline/{sample}.{raw_or_clean}.err"
    threads:
        48
    singularity:
        "docker://staphb/lyveset:2.0.1"
    shell:
        "date >> {output.log} ; " # time stamp, no version
        "if [ ! -f \"logs/genome_sizes.txt\" ] ; then wget https://raw.githubusercontent.com/StaPH-B/UPHL/master/genome_sizes.txt ; mv genome_sizes.txt logs/. ; fi ; "
        "mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d \"-\" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d \"_\" | cut -f 1 -d \".\" )) || true ; "
        "genome_length=$(grep $mash_result logs/genome_sizes.txt | grep -v \"#\" | head -n 1 | cut -f 2 -d \":\" | awk '{{ print $0 \"e+06\" }}' ) || genome_length=$(grep 'Estimated genome size:' {input.mash_error} | cut -f 4 -d \" \" ) || genome_length=\"0\" ; "
        "echo \"The genome length for {wildcards.sample} is $genome_length\" >> {output.log} ; "
        "run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e $genome_length 2>> {output.err} > {output.file} "
        "|| true ; touch {output}"

rule CG_pipeline_multiqc:
    input:
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE, raw_or_clean=['raw', 'clean'])
    output:
        file="cg-pipeline/cg-pipeline-summary.txt",
        log="logs/cg_pipeline_multiqc/log.log",
        err="logs/cg_pipeline_multiqc/log.err"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory=output_directory
    shell:
        "date >> {output.log} ; " # time stamp, no version
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
        log="logs/seqsero/{sample}.log",
        err="logs/seqsero/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/seqsero:1.0.1"
    shell:
        "date >> {output.log} ; " # time stamp
        "SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} 2>> {output.err} >> {output.log} "
        "|| true ; touch {output} ; "
        "cp {output.file} {output.final}"

rule seqsero_multiqc:
    input:
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE)
    output:
        file="SeqSero/Seqsero_serotype_results.txt",
        log="logs/seqsero_multiqc/log.log",
        err="logs/seqsero_multiqc/log.err"
    params:
        base_directory=base_directory,
        output_directory=output_directory
    threads:
        1
    shell:
        "date >> {output.log} ; " # time stamp
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
        log="logs/abricate/{sample}.{database}.log",
        err="logs/abricate/{sample}.{database}.err"
    threads:
        5
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        "date >> {output.log} ; " # time stamp
        "abricate --version >> {output.log} ; " # version of abricate
        "abricate --list >> {output.log} ; " # date of databases
        "abricate --db {wildcards.database} --threads {threads} {input} > {output.file} 2>> {output.err} "
        "|| true ; touch {output}"

rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
    output:
        file="logs/abricate_results/{database}.summary.txt",
        log="logs/abricate_summary/{database}.log",
        err="logs/abricate_summary/{database}.err"
    threads:
        1
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        "date >> {output.log} ; " # time stamp
        "abricate --version >> {output.log} ; " # version of abricate
        "abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output.file} 2>> {output.err} "
        "|| true ; touch {output}"

rule abricate_multiqc:
    input:
        rules.abricate_summary.output.file
    output:
        file="logs/abricate_results/{database}.summary.csv",
        log="logs/abricate_multiqc/{database}.log",
        err="logs/abricate_multiqc/{database}.err"
    threads:
        1
    shell:
        "date 2>> {output.err} | tee -a {output.log} ; "
        "cat {input} | "
        "sed 's/#//g' | sed 's/.tab//g' | sed \"s/{wildcards.database}.//g\" | "
        "awk '{{ sub(\"^.*/\", \"\", $1); print}}' | "
        "awk '{{ for (i=1;i<=NF;i++) if ($i ~ \";\" )gsub(\";.*$\",\"\",$i)g ; else continue}}{{print $0}}' | "
        "awk '{{ $2=\"\" ; print $0 }}' | sed 's/\\t/,/g' | sed 's/ /,/g' | "
        "sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' "
        "> {output.file} 2>> {output.err} "
        "|| true ; touch {output} "
