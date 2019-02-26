import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.2019.2.21")

base_directory=workflow.basedir
output_directory=os.getcwd()
seqyclean_adaptors="/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta"
mash_sketches="/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh"

SAMPLE, MIDDLE, EXTENSION = glob_wildcards('Sequencing_reads/Raw/{sample, [^_]+}_{middle}.f{extension}')
DATABASE = ['argannot', 'resfinder', 'card', 'plasmidfinder', 'vfdb', 'ecoli_vf', 'ncbi', 'ecoh', 'serotypefinder', 'cpd']

rule all:
    input:
        # copying files over
        expand("Sequencing_reads/Raw/{sample}_{middle}.f{extension}", zip, sample=SAMPLE, middle=MIDDLE, extension=EXTENSION),
        # running seqyclean
        expand("Sequencing_reads/QCed/{sample}_clean_PE1.fastq", sample=SAMPLE),
        expand("Sequencing_reads/QCed/{sample}_clean_PE2.fastq", sample=SAMPLE),
        "Sequencing_reads/Logs/seqyclean_summary.txt",
        # running FastQC
        expand("fastqc/{sample}_clean_{end}_fastqc.zip", sample=SAMPLE, end=['PE1', 'PE2', 'SE']),
        "fastqc/raw.complete",
        # running shovill
        expand("ALL_assembled/{sample}_contigs.fa",sample=SAMPLE),
        expand("ALL_assembled_plasmids/{sample}_plasmidcontigs.fa",sample=SAMPLE),
        # mash results
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE),
        "mash/mash_results.txt",
        # prokka results
        expand("ALL_gff/{sample}.gff", sample=SAMPLE),
        expand("ALL_gff_plasmids/{sample}_plasmid.gff", sample=SAMPLE),
        # quast results
        expand("quast/{sample}/report.txt", sample=SAMPLE),
        expand("quast_plasmids/{sample}_plasmid/report.txt", sample=SAMPLE),
        # seqsero results
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE),
        "SeqSero/Seqsero_serotype_results.txt",
        # cg-pipeline results
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE,raw_or_clean=['raw', 'clean']),
        "cg-pipeline/cg-pipeline-summary.txt",
        # abricate results
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results/{database}/{database}.summary.csv", database=DATABASE),
        expand("abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results/plasmids/{database}/{database}.plasmids.summary.csv", database=DATABASE),
    params:
        output_directory=output_directory,
        base_directory=workflow.basedir
    threads:
        48
    run:
        # creating a table from the benchmarks
        shell("{params.base_directory}/benchmark_multiqc.sh {params.output_directory}"),
        # getting the Summary
        shell("{params.base_directory}/check_multiqc.sh {params.output_directory}"),
        # running multiqc
        shell("which multiqc")
        shell("multiqc --version")
        shell("cp {params.base_directory}/multiqc_config_URF_snakemake.yaml multiqc_config.yaml"),
        shell("multiqc -f --ignore abricate_results/ --ignore beforeshovillupdate/ --outdir {params.output_directory}/logs --cl_config \"prokka_fn_snames: True\" {params.output_directory}"),

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
        "Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        "Sequencing_reads/QCed/{sample}_clean_PE2.fastq",
        "Sequencing_reads/QCed/{sample}_clean_SE.fastq",
        "Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.txt",
        "Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.tsv"
    params:
        seqyclean_adaptors
    threads:
        1
    log:
        "logs/seqyclean/{sample}.log"
    benchmark:
        "logs/benchmark/seqyclean/{sample}.log"
    run:
        shell("which seqyclean")
        shell("seqyclean | grep \"Version\"")
        shell("seqyclean -minlen 25 -qual -c {params} -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean")
        shell("mv Sequencing_reads/QCed/{wildcards.sample}_clean_SummaryStatistics* Sequencing_reads/Logs/.")

rule seqyclean_multiqc:
    input:
        expand("Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.tsv", sample=SAMPLE)
    output:
        "Sequencing_reads/Logs/seqyclean_summary.txt"
    log:
        "logs/seqyclean_multiqc/log.log"
    benchmark:
        "logs/benchmark/seqyclean_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory= workflow.basedir,
        output_directory= output_directory
    shell:
        "{params.base_directory}/seqyclean_multiqc.sh {params.output_directory}"

rule fastqc:
    input:
        "Sequencing_reads/QCed/{sample}_clean_{end}.fastq"
    output:
        "fastqc/{sample}_clean_{end}_fastqc.zip"
    log:
        "logs/fastqc/{sample}_clean_{end}.log"
    benchmark:
        "logs/benchmark/fastqc/{sample}_clean_{end}.log"
    threads:
        1
    run:
        shell("which fastqc")
        shell("fastqc --version")
        shell("fastqc --outdir fastqc --threads {threads} {input} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule fastqc_raw:
    input:
        get_reads
    output:
        "fastqc/raw.complete"
    log:
        "logs/fastqc/raw.log"
    benchmark:
        "logs/benchmark/fastqc/raw.log"
    threads:
        1
    run:
        shell("which fastqc")
        shell("fastqc --version")
        shell("fastqc --outdir fastqc --threads {threads} Sequencing_reads/Raw/*fastq* || true")
        shell("touch {output}")

rule shovill:
    input:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    threads:
        48
    log:
        "logs/shovill/{sample}.log"
    benchmark:
        "logs/benchmark/shovill/{sample}.log"
    output:
        "shovill_result/{sample}/contigs.fa",
    run:
        shell("which shovill")
        shell("shovill --version")
        shell("shovill --cpu {threads} --ram 200 --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule shovill_move:
    input:
        "shovill_result/{sample}/contigs.fa"
    log:
        "logs/shovill_move/{sample}.log"
    benchmark:
        "logs/benchmark/shovill_move/{sample}.log"
    threads:
        1
    output:
        "ALL_assembled/{sample}_contigs.fa"
    shell:
        "cp {input} {output}"

rule plasmidshovill:
    input:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    threads:
        48
    log:
        "logs/plasmidshovill/{sample}.log"
    benchmark:
        "logs/benchmark/plasmidshovill/{sample}.log"
    output:
        "shovill_result_plasmid/{sample}/contigs.fa",
    run:
        shell("which shovill")
        shell("shovill --version")
        shell("shovill --cpu {threads} --ram 200 --opts \"--plasmid\" --outdir shovill_result_plasmid/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule plasmidshovill_move:
    input:
        rules.plasmidshovill.output
    log:
        "logs/plasmidshovill_move/{sample}.log"
    benchmark:
        "logs/benchmark/plasmidshovill_move/{sample}.log"
    threads:
        1
    output:
        "ALL_assembled_plasmids/{sample}_plasmidcontigs.fa"
    shell:
        "cp {input} {output}"

rule mash_cat:
    input:
        "Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        "Sequencing_reads/QCed/{sample}_clean_PE2.fastq",
        "Sequencing_reads/QCed/{sample}_clean_SE.fastq"
    output:
        "mash/{sample}.clean_all.fastq"
    log:
        "logs/mash_cat/{sample}.log"
    benchmark:
        "logs/benchmark/mash_cat/{sample}.log"
    threads:
        1
    shell:
        "cat {input} > {output}"

rule mash_sketch:
    input:
        rules.mash_cat.output
    output:
        "mash/{sample}.clean_all.fastq.msh"
    log:
        "logs/mash_sketch/{sample}.log"
    benchmark:
        "logs/benchmark/mash_sketch/{sample}.log"
    threads:
        1
    run:
        shell("which mash")
        shell("mash --version")
        shell("mash sketch -m 2 {input} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule mash_dist:
    input:
        rules.mash_sketch.output
    output:
        "mash/{sample}.clean_all.fastq.msh.distance.txt"
    threads:
        48
    params:
        mash_sketches
    log:
        "logs/mash_dist/{sample}.log"
    benchmark:
        "logs/benchmark/mash_dist/{sample}.log"
    run:
        shell("which mash")
        shell("mash --version")
        shell("mash dist -p {threads} {params} {input} > {output} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule mash_sort:
    input:
        rules.mash_dist.output
    output:
        "mash/{sample}.clean_all.fastq.msh.distance.sorted.txt"
    log:
        "logs/mash_sort/{sample}.log"
    benchmark:
        "logs/benchmark/mash_sort/{sample}.log"
    threads:
        1
    shell:
        "sort -gk3 {input} > {output} "

rule mash_multiqc:
    input:
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE)
    output:
        "mash/mash_results.txt"
    log:
        "logs/mash_pipeline_multiqc/log.log"
    benchmark:
        "logs/benchmark/mash_pipeline_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory= workflow.basedir,
        output_directory= output_directory
    shell:
        "{params.base_directory}/mash_multiqc.sh {params.output_directory}"

rule prokka:
    input:
        contig_file=rules.shovill_move.output,
        mash_file=rules.mash_sort.output
    threads:
        48
    output:
        "Prokka/{sample}/{sample}.gff",
    log:
        "logs/prokka/{sample}.log"
    benchmark:
        "logs/benchmark/prokka/{sample}.log"
    shell:
        """
        which prokka
        prokka --version
        if [ -s \"{input.mash_file}\" ]
        then
            mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | awk -F \"-.-\" '{{ print $NF }}' | sed 's/.fna//g' | awk -F \"_\" '{{ print $1 \" \" $2 }}' ))
            prokka --cpu {threads} --compliant --centre --UPHL --mincontiglen 500 --outdir Prokka/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample} --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} || true
        fi
        if [ ! -f {output} ] ; then touch {output} ; fi
        """

rule prokka_move:
    input:
        rules.prokka.output
    output:
        "ALL_gff/{sample}.gff"
    log:
        "logs/prokka_move/{sample}.log"
    benchmark:
        "logs/benchmark/prokka_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output}"

rule plasmidprokka:
    input:
        contig_file=rules.plasmidshovill_move.output,
        mash_file=rules.mash_sort.output
    threads:
        48
    output:
        "Prokka_plasmids/{sample}/{sample}_plasmid.gff",
    log:
        "logs/plasmidprokka/{sample}.log"
    benchmark:
        "logs/benchmark/plasmidprokka/{sample}.log"
    shell:
        """
        which prokka
        prokka --version
        if [ -s \"{input.mash_file}\" ]
        then
            mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | awk -F \"-.-\" '{{ print $NF }}' | sed 's/.fna//g' | awk -F \"_\" '{{ print $1 \" \" $2 }}' ))
            prokka --cpu {threads} --compliant --centre --UPHL --mincontiglen 500 --outdir Prokka_plasmids/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample}_plasmid --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} || true
        fi
        if [ ! -f {output} ] ; then touch {output} ; fi
        """

rule plasmidprokka_move:
    input:
        rules.plasmidprokka.output
    output:
        "ALL_gff_plasmids/{sample}_plasmid.gff"
    log:
        "logs/plasmidprokka_move/{sample}.log"
    benchmark:
        "logs/benchmark/plasmidprokka_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output}"

rule quast:
    input:
        rules.shovill_move.output
    output:
        "quast/{sample}/report.txt",
    log:
        "logs/quast/{sample}.log"
    benchmark:
        "logs/benchmark/quast/{sample}.log"
    threads:
        1
    run:
        shell("which quast")
        shell("quast --version")
        shell("quast {input} --output-dir quast/{wildcards.sample} --threads {threads} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule plasmidquast:
    input:
        rules.plasmidshovill_move.output
    output:
        "quast_plasmids/{sample}_plasmid/report.txt",
    log:
        "logs/plasmidquast/{sample}.log"
    benchmark:
        "logs/benchmark/plasmidquast/{sample}.log"
    threads:
        1
    run:
        shell("which quast")
        shell("quast --version")
        shell("quast {input} --output-dir quast_plasmids/{wildcards.sample}_plasmid --threads {threads} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule GC_pipeline_shuffle_raw:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        "Sequencing_reads/shuffled/{sample}_raw_shuffled.fastq.gz"
    log:
        "logs/gc_pipeline_shuffle_raw/{sample}.log"
    benchmark:
        "logs/benchmark/gc_pipeline_shuffle_raw/{sample}.log"
    threads:
        1
    run:
        shell("which run_assembly_shuffleReads.pl")
        shell("run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output}")

rule GC_pipeline_shuffle_clean:
    input:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    output:
        "Sequencing_reads/shuffled/{sample}_clean_shuffled.fastq.gz"
    log:
        "logs/gc_pipeline_shuffle_clean/{sample}.log"
    benchmark:
        "logs/benchmark/gc_pipeline_shuffle_clean/{sample}.log"
    threads:
        1
    run:
        shell("which run_assembly_shuffleReads.pl")
        shell("run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output}")

rule GC_pipeline:
    input:
        shuffled_fastq="Sequencing_reads/shuffled/{sample}_{raw_or_clean}_shuffled.fastq.gz",
        quast_file="quast/{sample}/report.txt"
    output:
        "cg-pipeline/{sample}.{raw_or_clean}.out.txt"
    threads:
        48
    log:
        "logs/gc_pipeline/{sample}.{raw_or_clean}.log"
    benchmark:
        "logs/benchmark/gc_pipeline/{sample}.{raw_or_clean}.log"
    params:
        base_directory=workflow.basedir
    run:
        shell("which run_assembly_readMetrics.pl")
        if "PNUSAS" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 5000000 > {output} || true")
        elif "PNUSAE" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 5000000 > {output} || true")
        elif "PNUSAC" in wildcards.sample :
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 1600000 > {output} || true")
        elif "PNUSAL" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 3000000 > {output} || true")
        else:
            shell("{params.base_directory}/genome_length_cg.sh {input.shuffled_fastq} {input.quast_file} {threads} {output} || true")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi")

rule GC_pipeline_multiqc:
    input:
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE, raw_or_clean=['raw', 'clean'])
    output:
        "cg-pipeline/cg-pipeline-summary.txt"
    log:
        "logs/gc_pipeline_multiqc/log.log"
    benchmark:
        "logs/benchmark/gc_pipeline_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory=workflow.basedir,
        output_directory=output_directory
    shell:
        "{params.base_directory}/cgpipeline_multiqc.sh {params.output_directory}"

rule seqsero:
    input:
        "Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        "Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    output:
        "SeqSero/{sample}/Seqsero_result.txt",
        "SeqSero/{sample}/data_log.txt"
    log:
        "logs/seqsero/{sample}.log"
    benchmark:
        "logs/benchmark/seqsero/{sample}.log"
    threads:
        1
    run:
        shell("which SeqSero.py")
        shell("SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input}")

rule seqsero_move:
    input:
        "SeqSero/{sample}/Seqsero_result.txt"
    output:
        "SeqSero/{sample}.Seqsero_result.txt"
    log:
        "logs/seqsero_move/{sample}.log"
    benchmark:
        "logs/benchmark/seqsero_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output}"

rule seqsero_multiqc:
    input:
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE)
    output:
        "SeqSero/Seqsero_serotype_results.txt"
    log:
        "logs/seqsero_multiqc/log.log"
    benchmark:
        "logs/benchmark/seqsero_multiqc/benchmark.log"
    params:
        base_directory=workflow.basedir,
        output_directory=output_directory
    threads:
        1
    shell:
        "{params.base_directory}/seqsero_multiqc.sh {params.output_directory}"

rule abricate:
    input:
        rules.shovill_move.output
    output:
        "abricate_results/{database}/{database}.{sample}.out.tab"
    log:
        "logs/abricate/{sample}.{database}.log"
    benchmark:
        "logs/benchmark/abricate_{database}/{sample}.log"
    threads:
        48
    run:
        shell("which abricate")
        shell("abricate --version")
        shell("abricate --list")
        shell("abricate --db {wildcards.database} --threads {threads} {input} > {output} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule abricate_combine:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE)
    output:
        summary="abricate_results/{database}/{database}.summary.txt",
        results="logs/{database}.summary.txt"
    log:
        "logs/abricate_combine_{database}/{database}.log"
    benchmark:
        "logs/benchmark/abricate_combine/{database}.log"
    threads:
        1
    run:
        shell("which abricate")
        shell("abricate --version")
        shell("abricate --summary abricate_results/{wildcards.database}/{wildcards.database}*tab > {output.summary}")
        shell("cp {output.summary} {output.results}")

rule abricate_multiqc:
    input:
        rules.abricate_combine.output.summary
    output:
        "abricate_results/{database}/{database}.summary.csv"
    log:
        "logs/abricate_multiqc/{database}.log"
    benchmark:
        "logs/benchmark/abricate_multiqc/{database}.log"
    threads:
        1
    shell:
        "cat {input} | "
        "sed 's/#//g' | sed 's/.tab//g' | sed \"s/{wildcards.database}.//g\" | "
        "awk '{{ sub(\"^.*/\", \"\", $1); print}}' | "
        "awk '{{ for (i=1;i<=NF;i++) if ($i ~ \";\" )gsub(\";.*$\",\"\",$i)g ; else continue}}{{print $0}}' | "
        "awk '{{ $2=\"\" ; print $0 }}' | sed 's/\t/,/g' | sed 's/ /,/g' | "
        "sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' "
        "> {output}"

rule plasmidabricate:
    input:
        rules.plasmidshovill_move.output
    output:
        "abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab"
    log:
        "logs/plasmidabricate/{sample}.{database}.log"
    benchmark:
        "logs/benchmark/plasmidabricate_{database}/{sample}.log"
    threads:
        48
    run:
        shell("which abricate")
        shell("abricate --version")
        shell("abricate --list")
        shell("abricate --db {wildcards.database} --threads {threads} {input} > {output} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule plasmidabricate_combine:
    input:
        expand("abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab", sample=SAMPLE, database=DATABASE)
    output:
        summary="abricate_results_plasmids/{database}/{database}.plasmids.summary.txt",
        results="logs/{database}.plasmids.summary.txt"
    log:
        "logs/plasmidabricate_combine_{database}/{database}.log"
    benchmark:
        "logs/benchmark/plasmidabricate_combine/{database}.log"
    threads:
        1
    run:
        shell("which abricate")
        shell("abricate --version")
        shell("abricate --summary abricate_results/plasmids/{wildcards.database}/{wildcards.database}*tab > {output.summary}")
        shell("cp {output.summary} {output.results}")

rule plasmidabricate_multiqc:
    input:
        rules.plasmidabricate_combine.output.summary
    output:
        "abricate_results/plasmids/{database}/{database}.plasmids.summary.csv"
    log:
        "logs/plasmidabricate_multiqc/{database}.log"
    benchmark:
        "logs/benchmark/plasmidabricate_multiqc/{database}.log"
    threads:
        1
    shell:
        "cat {input} | "
        "sed 's/#//g' | sed 's/.tab//g' | sed \"s/{wildcards.database}.//g\" | "
        "awk '{{ sub(\"^.*/\", \"\", $1); print}}' | "
        "awk '{{ for (i=1;i<=NF;i++) if ($i ~ \";\" )gsub(\";.*$\",\"\",$i)g ; else continue}}{{print $0}}' | "
        "awk '{{ $2=\"\" ; print $0 }}' | sed 's/\t/,/g' | sed 's/ /,/g' | "
        "sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' "
        "> {output}"
