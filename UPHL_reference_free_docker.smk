import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.2019.06.27")

base_directory=workflow.basedir + "/URF_scripts"
output_directory=os.getcwd()

SAMPLE, MIDDLE, EXTENSION = glob_wildcards('Sequencing_reads/Raw/{sample, [^_]+}_{middle}.f{extension}')
DATABASE = ['argannot', 'card', 'ecoh', 'ecoli_vf', 'ncbi', 'plasmidfinder', 'vfdb' ]

rule all:
    input:
        # copying files over
        expand("Sequencing_reads/Raw/{sample}_{middle}.f{extension}", zip, sample=SAMPLE, middle=MIDDLE, extension=EXTENSION),
        # running seqyclean
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
        #"results_for_multiqc/seqyclean_summary.txt",
        # running FastQC
        "fastqc/fastqc.complete",
        # running shovill
        expand("shovill_result/{sample}/contigs.fa", sample=SAMPLE),
        expand("ALL_assembled/{sample}_contigs.fa", sample=SAMPLE),
        # mash results
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE),
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
    params:
        output_directory=output_directory,
        base_directory=base_directory
    threads:
        48
    run:
        # getting the Summary
        shell("{params.base_directory}/check_multiqc.sh {params.output_directory} 2>> logs/all/all.err | tee -a logs/all/all.log || true "),
        shell("ln -s {params.output_directory}/fastqc {params.output_directory}/results_for_multiqc/fastqc 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp Prokka*/*/*txt results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp SeqSero/Seqsero_serotype_results.txt results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp mash/mash_results.txt results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp cg-pipeline/cg-pipeline-summary.txt results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp logs/abricate_results/*.summary.csv results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("ln -s {params.output_directory}/quast {params.output_directory}/results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp logs/File_heatmap.csv  results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp logs/raw_clean_coverage.txt  results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cp logs/raw_clean_scatter.csv  results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        shell("cat run_results_summary.txt | sed 's/simple_mash_result/A.simple_mash_result/g' | sed 's/simple_seqsero_result/B.simple_seqsero_result/g' | sed 's/abricate_serotype_O/C.abricate_serotype_O/g' | sed 's/abricate_serotype_H/D.abricate_serotype_H/g' | sed 's/fastqc_raw_reads_2/E.fastqc_raw_reads_2/g' | sed 's/fastqc_clean_reads_PE2/F.fastqc_clean_reads_PE2/g' | sed 's/cg_raw_coverage/G.cg_raw_coverage/g' | sed 's/cg_cln_coverage/H.cg_cln_coverage/g' | sed 's/ncbi/J.ncbi_antibiotic_resistence_genes/g' | sed 's/stxeae_result/I.stx_and_eae_virulence_factor_result/g' > results_for_multiqc/run_results_summary.txt || true ")
        shell("cp run_file_summary.txt results_for_multiqc/. 2>> logs/all/all.err | tee -a logs/all/all.log || true ")
        # running multiqc
        shell("which multiqc 2>> logs/all/all.err | tee -a logs/all/all.log")
        shell("multiqc --version 2>> logs/all/all.err | tee -a logs/all/all.log")
        shell("cp {params.base_directory}/multiqc_config_URF_snakemake.yaml multiqc_config.yaml 2>> logs/all/all.err | tee -a logs/all/all.log"),
        shell("multiqc -f --outdir {params.output_directory}/logs --cl_config \"prokka_fn_snames: True\" {params.output_directory}/results_for_multiqc 2>> logs/all/all.err | tee -a logs/all/all.log || true"),

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
        "docker://staphb/seqyclean:latest"
    shell:
        "seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean "
        "2>> {output.err} | tee -a {output.log} "
        #"|| true ; touch {output}"

#Seqyclean Module pending
#rule seqyclean_multiqc:
#    input:
#        expand("Sequencing_reads/QCed/{sample}_clean_SummaryStatistics.tsv", sample=SAMPLE)
#    output:
#        file="results_for_multiqc/seqyclean_summary.txt",
#        err="logs/seqyclean_multiqc/log.err",
#        log="logs/seqyclean_multiqc/log.log"
#    threads:
#        1
#    shell:
#        "grep \"Version\" Sequencing_reads/QCed/*clean_SummaryStatistics.tsv | head -n 1 | cut -f 2- -d \':\' | cut -f 2- | awk '{{ print \"Sample,\" $0 }}' | sed 's/PE1PE2/PE1,PE2/g' | sed 's/ /,/g' | sed 's/\\t/,/g' | sed 's/,,/,/g'  2>> logs/seqyclean_multiqc/log.err | tee -a {output.file} ; "
#        "grep \"Sequencing_reads\" Sequencing_reads/QCed/*_clean_SummaryStatistics.tsv | cut -f 2- -d \':\' | sort | uniq | cut -f 2- | awk \'{{ print $1 \" \" $0 }}\' | awk '{{ sub(\"^.*/\", \"\", $1); print $0 }}' | awk '{{ sub(\"_.*fastq.*,\", \"\", $1 ); print $0 }}' | sed 's/ /,/g' | sed 's/\\t/,/g' | sed 's/,,/,/g' 2>> logs/seqyclean_multiqc/log.err | tee -a {output.file} ; "
#        "touch {output}"

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
        "docker://dukegcb/fastqc:latest"
    shell:
        "fastqc --outdir fastqc --threads {threads} Sequencing_reads/*/*.fastq* 2>> {output.err} | tee -a {output.log} || true ; touch {output}"

rule shovill:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    threads:
        48
#TBD    resources:
#        200
    output:
        file="shovill_result/{sample}/contigs.fa",
        log="logs/shovill/{sample}.log",
        err="logs/shovill/{sample}.err"
    singularity:
        "docker://staphb/shovill:latest"
    shell:
        "shovill --cpu {threads} --ram 200 --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force "
        "2>> {output.err} | tee -a {output.log} "
        #"|| true ; touch {output} "

rule shovill_move:
    input:
        rules.shovill.output
    threads:
        1
    output:
        file="ALL_assembled/{sample}_contigs.fa",
        err="logs/shovill_move/{sample}.err",
        log="logs/shovill_move/{sample}.log"
    shell:
        "cp {input} {output.file} 2>> {output.err} | tee -a {output.log} ; "
        "touch {output}"

rule mash_cat:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read1,
        rules.seqyclean.output.se
    output:
        file="mash/{sample}.clean_all.fastq",
        err="logs/mash_cat/{sample}.err"
    threads:
        1
    shell:
        "cat {input} > {output.file} 2>> {output.err}"

rule mash_sketch:
    input:
        rules.mash_cat.output
    output:
        file="mash/{sample}.clean_all.fastq.msh",
        log="logs/mash_sketch/{sample}.log",
        err="logs/mash_sketch/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/mash:latest"
    shell:
        "mash sketch -m 2 {input} 2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output.file}"

rule mash_dist:
    input:
        rules.mash_sketch.output
    output:
        file="mash/{sample}.clean_all.fastq.msh.distance.txt",
        err="logs/mash_dist/{sample}.err"
    singularity:
        "docker://staphb/mash:latest"
    shell:
        "mash dist -p {threads} /db/RefSeqSketchesDefaults.msh {input} > {output.file} "
        "2>> {output.err} || true ; touch {output} "

rule mash_sort:
    input:
        rules.mash_dist.output
    output:
        file="mash/{sample}.clean_all.fastq.msh.distance.sorted.txt",
        err="logs/mash_sort/{sample}.err"
    threads:
        1
    shell:
        "sort -gk3 {input} > {output.file} 2>> {output.err} "

rule mash_multiqc:
    input:
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE)
    output:
        file="mash/mash_results.txt",
        log="logs/mash_pipeline_multiqc/log.log",
        err="logs/mash_pipeline_multiqc/log.err"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory= output_directory
    shell:
        "{params.base_directory}/mash_multiqc.sh {params.output_directory} 2>> {output.err} | tee -a {output.log}"

rule prokka:
    input:
        contig_file=rules.shovill_move.output,
        mash_file=rules.mash_sort.output
    threads:
        48
    output:
        file="Prokka/{sample}/{sample}.gff",
        log="logs/prokka/{sample}.log",
        err="logs/prokka/{sample}.err"
    singularity:
        "docker://staphb/prokka:latest"
    shell:
        """
        if [ -s \"{input.mash_file}\" ]
        then
        mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | awk -F \"-.-\" '{{ print $NF }}' | sed 's/.fna//g' | awk -F \"_\" '{{ print $1 \" \" $2 }}' ))
        prokka --cpu {threads} --compliant --centre --URF --mincontiglen 500 --outdir Prokka/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample} --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} 2>> {output.err} | tee -a {output.log} || true
        fi
        touch {output}
        """

rule prokka_move:
    input:
        rules.prokka.output
    output:
        file="ALL_gff/{sample}.gff",
        log="logs/prokka_move/{sample}.log",
        err="logs/prokka_move/{sample}.err"
    threads:
        1
    shell:
        "cp {input} {output.file} 2>> {output.err} | tee -a {output.log}"

rule quast:
    input:
        rules.shovill_move.output
    output:
        file="quast/{sample}/report.txt",
        log="logs/quast/{sample}.log",
        err="logs/quast/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/quast:latest"
    shell:
        "quast.py {input} --output-dir quast/{wildcards.sample} --threads {threads} "
        "2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} "

rule CG_pipeline_shuffle_raw:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        file="Sequencing_reads/shuffled/{sample}_raw_shuffled.fastq.gz",
        err="logs/cg_pipeline_shuffle_raw/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/lyveset:latest"
    shell:
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.err}"

rule CG_pipeline_shuffle_clean:
    input:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    output:
        file="Sequencing_reads/shuffled/{sample}_clean_shuffled.fastq.gz",
        err="logs/cg_pipeline_shuffle_clean/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/lyveset:latest"
    shell:
        "run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.err} "
        "|| true ; touch {output} "

rule CG_pipeline:
    input:
        shuffled_fastq="Sequencing_reads/shuffled/{sample}_{raw_or_clean}_shuffled.fastq.gz",
        quast_file="quast/{sample}/report.txt"
    output:
        file="cg-pipeline/{sample}.{raw_or_clean}.out.txt",
        err="logs/cg_pipeline/{sample}.{raw_or_clean}.err"
    threads:
        48
    params:
        base_directory=base_directory
    singularity:
        "docker://staphb/lyveset:latest"
    shell:
        "run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e $(grep 'Total length (>= 0 bp)' {input.quast_file} | grep -v \"All statistics\" | sed 's/Total length (>= 0 bp)//g' | sed 's/ //g' ) 2>> {output.err} > {output.file} "
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
        "{params.base_directory}/cgpipeline_multiqc.sh {params.output_directory} "
        "2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} "

rule seqsero:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read2
    output:
        results="SeqSero/{sample}/Seqsero_result.txt",
        datalog="SeqSero/{sample}/data_log.txt",
        log="logs/seqsero/{sample}.log",
        err="logs/seqsero/{sample}.err"
    threads:
        1
    singularity:
        "docker://staphb/seqsero:latest"
    shell:
        "SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} 2>> {output.err} | tee -a {output.log} "
        "|| true ; touch {output} "

rule seqsero_move:
    input:
        rules.seqsero.output.results
    output:
        file="SeqSero/{sample}.Seqsero_result.txt",
        log="logs/seqsero_move/{sample}.log",
        err="logs/seqsero_move/{sample}.err"
    threads:
        1
    shell:
        "cp {input} {output.file} 2>> {output.err} | tee -a {output.log}"

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
        "{params.base_directory}/seqsero_multiqc.sh {params.output_directory} "
        "2>> {output.err} | tee -a {output.log}"

rule abricate:
    input:
        rules.shovill_move.output
    output:
        file="abricate_results/{database}/{database}.{sample}.out.tab",
        err="logs/abricate/{sample}.{database}.err"
    threads:
        5
    singularity:
        "docker://staphb/abricate:latest"
    shell:
        "abricate --db {wildcards.database} --threads {threads} {input} > {output.file} 2>> {output.err} "
        "|| true ; touch {output}"

rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
    output:
        file="logs/abricate_results/{database}.summary.txt",
        err="logs/abricate_summary/{database}.err"
    threads:
        1
    singularity:
        "docker://staphb/abricate:latest"
    shell:
        "abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output.file} 2>> {output.err} "
        "|| true ; touch {output}"

rule abricate_multiqc:
    input:
        rules.abricate_summary.output
    output:
        file="logs/abricate_results/{database}.summary.csv",
        err="logs/abricate_multiqc/{database}.err"
    threads:
        1
    shell:
        "cat {input} | "
        "sed 's/#//g' | sed 's/.tab//g' | sed \"s/{wildcards.database}.//g\" | "
        "awk '{{ sub(\"^.*/\", \"\", $1); print}}' | "
        "awk '{{ for (i=1;i<=NF;i++) if ($i ~ \";\" )gsub(\";.*$\",\"\",$i)g ; else continue}}{{print $0}}' | "
        "awk '{{ $2=\"\" ; print $0 }}' | sed 's/\\t/,/g' | sed 's/ /,/g' | "
        "sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' "
        "> {output.file} 2>> {output.err} "
        "|| true ; touch {output} "
