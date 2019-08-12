import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.2019.8.14")

base_directory=workflow.basedir + "/URF_scripts"
output_directory=os.getcwd()
seqyclean_adaptors="/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta"
mash_sketches="/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh"

SAMPLE, MIDDLE, EXTENSION = glob_wildcards('Sequencing_reads/Raw/{sample, [^_]+}_{middle}.f{extension}')
DATABASE = ['ncbi', 'vfdb', 'serotypefinder']
#DATABASE = ['argannot', 'bacmet2', 'card', 'ecoh', 'ecoli_vf', 'ncbi', 'plasmidfinder', 'vfdb', 'serotypefinder', 'cpd']

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
        expand("shovill_result/{sample}/contigs.fa", sample=SAMPLE),
        expand("shovill_result_plasmid/{sample}/contigs.fa", sample=SAMPLE),
        expand("ALL_assembled/{sample}_contigs.fa", sample=SAMPLE),
        expand("ALL_assembled_plasmids/{sample}_plasmidcontigs.fa", sample=SAMPLE),
        expand("results_for_multiqc/{sample}_flash.hist", sample=SAMPLE),
        # mash results
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE),
        "mash/mash_results.txt",
        # mash screen
#        expand("mash/{sample}_mashscreen.txt", sample=SAMPLE),
#        "mash/mash_screen_results.txt",
        # prokka results
        expand("Prokka/{sample}/{sample}.gff", sample=SAMPLE),
        expand("Prokka_plasmids/{sample}/{sample}_plasmid.gff", sample=SAMPLE),
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
        expand("abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab", sample=SAMPLE, database=DATABASE),
        expand("logs/abricate_results/{database}.summary.csv", database=DATABASE),
        #blobtools
        expand("bwa/{sample}.sorted.bam", sample=SAMPLE),
        expand("blast/{sample}.tsv", sample=SAMPLE),
#        expand("blobtools/{sample}.blobDB.json", sample=SAMPLE),
#        expand("blobtools/{sample}.blobDB.table.txt", sample=SAMPLE),
#        expand("blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png", sample=SAMPLE),
    params:
        output_directory=output_directory,
        base_directory=base_directory
    log:
        out="logs/all/all.log",
        err="logs/all/all.err"
    threads:
        48
    run:
        # creating a table from the benchmarks
        shell("{params.base_directory}/benchmark_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"),
        # getting the Summary
        shell("{params.base_directory}/check_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"),
        # getting all multiqc files into logs
        shell("mkdir -p results_for_multiqc 2>> {log.err} | tee -a {log.out}")
        shell("cp mash/mash*results.txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out}")
        shell("cp cg-pipeline/cg-pipeline-summary.txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out}")
        shell("cp SeqSero/Seqsero_serotype_results*.txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out}")
        shell("cp Sequencing_reads/Logs/seqyclean_summary.txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out}")
        shell("ln -s {params.output_directory}/fastqc {params.output_directory}/results_for_multiqc/fastqc 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp Prokka*/*/*txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("ln -s {params.output_directory}/quast {params.output_directory}/results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("ln -s {params.output_directory}/quast_plasmids {params.output_directory}/results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp logs/abricate_results/*.summary.csv results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp logs/benchmark_summary.csv results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp logs/File_heatmap.csv  results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp logs/raw_clean_coverage.txt  results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cp logs/raw_clean_scatter.csv  results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        shell("cat run_results_summary.txt | sed 's/simple_mash_result/A.simple_mash_result/g' | sed 's/simple_seqsero_result/B.simple_seqsero_result/g' | sed 's/abricate_serotype_O/C.abricate_serotype_O/g' | sed 's/abricate_serotype_H/D.abricate_serotype_H/g' | sed 's/fastqc_raw_reads_2/E.fastqc_raw_reads_2/g' | sed 's/fastqc_clean_reads_PE2/F.fastqc_clean_reads_PE2/g' | sed 's/cg_raw_coverage/G.cg_raw_coverage/g' | sed 's/cg_cln_coverage/H.cg_cln_coverage/g' | sed 's/ncbi/J.ncbi_antibiotic_resistence_genes/g' | sed 's/stxeae_result/I.stx_and_eae_virulence_factor_result/g' > results_for_multiqc/run_results_summary.txt || true ")
        shell("cp run_file_summary.txt results_for_multiqc/. 2>> {log.err} | tee -a {log.out} || true ")
        # running multiqc
        shell("which multiqc 2>> {log.err} | tee -a {log.out}")
        shell("multiqc --version 2>> {log.err} | tee -a {log.out}")
        shell("cp {params.base_directory}/multiqc_config_URF_snakemake.yaml multiqc_config.yaml 2>> {log.err} | tee -a {log.out}"),
        shell("multiqc -f --outdir {params.output_directory}/logs --cl_config \"prokka_fn_snames: True\" {params.output_directory}/results_for_multiqc 2>> {log.err} | tee -a {log.out} || true"),

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
        redse="Sequencing_reads/QCed/{sample}_clean_SE.fastq",
        sstxt="Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.txt",
        sstsv="Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.tsv"
    params:
        seqyclean_adaptors
    threads:
        1
    log:
        out="logs/seqyclean/{sample}.log",
        err="logs/seqyclean/{sample}.err"
    benchmark:
        "logs/benchmark/seqyclean/{sample}.log"
    run:
        shell("which seqyclean 2>> {log.err} | tee -a {log.out}")
        shell("seqyclean | grep \"Version\" 2>> {log.err} | tee -a {log.out}")
        shell("seqyclean -minlen 25 -qual -c {params} -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean 2>> {log.err} | tee -a {log.out}")
        shell("mv Sequencing_reads/QCed/{wildcards.sample}_clean_SummaryStatistics* Sequencing_reads/Logs/. 2>> {log.err} | tee -a {log.out}")

rule seqyclean_multiqc:
    input:
        expand("Sequencing_reads/Logs/{sample}_clean_SummaryStatistics.tsv", sample=SAMPLE)
    output:
        "Sequencing_reads/Logs/seqyclean_summary.txt"
    log:
        out="logs/seqyclean_multiqc/log.log",
        err="logs/seqyclean_multiqc/log.err"
    benchmark:
        "logs/benchmark/seqyclean_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory= output_directory
    shell:
        "{params.base_directory}/seqyclean_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"

rule fastqc:
    input:
        "Sequencing_reads/QCed/{sample}_clean_{end}.fastq"
    output:
        "fastqc/{sample}_clean_{end}_fastqc.zip"
    log:
        out="logs/fastqc/{sample}_clean_{end}.log",
        err="logs/fastqc/{sample}_clean_{end}.err"
    threads:
        1
    run:
        shell("which fastqc 2>> {log.err} | tee -a {log.out}")
        shell("fastqc --version 2>> {log.err} | tee -a {log.out}")
        shell("fastqc --outdir fastqc --threads {threads} {input} 2>> {log.err} | tee -a {log.out} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule fastqc_raw:
    input:
        glob.glob("Sequencing_reads/Raw/*fastq*")
#        get_reads
    output:
        "fastqc/raw.complete"
    log:
        out="logs/fastqc/raw.log",
        err="logs/fastqc/raw.err"
    threads:
        1
    run:
        shell("which fastqc 2>> {log.err} | tee -a {log.out}")
        shell("fastqc --version 2>> {log.err} | tee -a {log.out}")
        shell("fastqc --outdir fastqc --threads {threads} Sequencing_reads/Raw/*fastq* 2>> {log.err} | tee -a {log.out} || true")
        shell("touch {output}")

rule shovill:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    threads:
        48
    log:
        out="logs/shovill/{sample}.log",
        err="logs/shovill/{sample}.err"
    benchmark:
        "logs/benchmark/shovill/{sample}.log"
    output:
        "shovill_result/{sample}/contigs.fa",
    run:
        shell("which shovill 2>> {log.err} | tee -a {log.out}")
        shell("shovill --version 2>> {log.err} | tee -a {log.out}")
        shell("shovill --cpu {threads} --ram 200 --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force 2>> {log.err} | tee -a {log.out} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule shovill_move:
    input:
        rules.shovill.output
    log:
        out="logs/shovill_move/{sample}.log",
        err="logs/shovill_move/{sample}.err"
    benchmark:
        "logs/benchmark/shovill_move/{sample}.log"
    threads:
        1
    output:
        "ALL_assembled/{sample}_contigs.fa"
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule plasmidshovill:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    threads:
        10
    log:
        out="logs/plasmidshovill/{sample}.log",
        err="logs/plasmidshovill/{sample}.err"
    benchmark:
        "logs/benchmark/plasmidshovill/{sample}.log"
    output:
        contig="shovill_result_plasmid/{sample}/contigs.fa",
        flash="shovill_result_plasmid/{sample}/flash.hist"
    run:
        shell("which shovill 2>> {log.err} | tee -a {log.out}")
        shell("shovill --version 2>> {log.err} | tee -a {log.out}")
        shell("shovill --cpu {threads} --ram 200 --opts \"--plasmid\" --outdir shovill_result_plasmid/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force 2>> {log.err} | tee -a {log.out} || true ; touch {output}")

rule plasmidshovill_move:
    input:
        rules.plasmidshovill.output.contig
    log:
        out="logs/plasmidshovill_move/{sample}.log",
        err="logs/plasmidshovill_move/{sample}.err"
    benchmark:
        "logs/benchmark/plasmidshovill_move/{sample}.log"
    threads:
        1
    output:
        "ALL_assembled_plasmids/{sample}_plasmidcontigs.fa"
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule mash_cat:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read2,
        rules.seqyclean.output.redse
    output:
        "mash/{sample}.clean_all.fastq"
    log:
        err="logs/mash_cat/{sample}.err"
    benchmark:
        "logs/benchmark/mash_cat/{sample}.log"
    threads:
        1
    shell:
        "cat {input} > {output} 2>> {log.err}"

rule mash_sketch:
    input:
        rules.mash_cat.output
    output:
        "mash/{sample}.clean_all.fastq.msh"
    log:
        out="logs/mash_sketch/{sample}.log",
        err="logs/mash_sketch/{sample}.err"
    benchmark:
        "logs/benchmark/mash_sketch/{sample}.log"
    threads:
        1
    run:
        shell("which mash 2>> {log.err} | tee -a {log.out}")
        shell("mash --version 2>> {log.err} | tee -a {log.out}")
        shell("mash sketch -m 2 {input} 2>> {log.err} | tee -a {log.out} || true ; touch {output} ")

rule mash_dist:
    input:
        rules.mash_sketch.output
    output:
        "mash/{sample}.clean_all.fastq.msh.distance.txt"
    threads:
        5
    params:
        mash_sketches
    log:
        out="logs/mash_dist/{sample}.log",
        err="logs/mash_dist/{sample}.err"
    benchmark:
        "logs/benchmark/mash_dist/{sample}.log"
    run:
        shell("which mash 2>> {log.err} | tee -a {log.out}")
        shell("mash --version 2>> {log.err} | tee -a {log.out}")
        shell("mash dist -p {threads} {params} {input} > {output} 2>> {log.err} || true ; touch {output} ")

rule mash_sort:
    input:
        rules.mash_dist.output
    output:
        "mash/{sample}.clean_all.fastq.msh.distance.sorted.txt"
    log:
        err="logs/mash_sort/{sample}.err"
    benchmark:
        "logs/benchmark/mash_sort/{sample}.log"
    threads:
        1
    shell:
        "sort -gk3 {input} > {output} 2>> {log.err} "

rule mash_multiqc:
    input:
        expand("mash/{sample}.clean_all.fastq.msh.distance.sorted.txt", sample=SAMPLE)
    output:
        "mash/mash_results.txt"
    log:
        out="logs/mash_pipeline_multiqc/log.log",
        err="logs/mash_pipeline_multiqc/log.err"
    benchmark:
        "logs/benchmark/mash_pipeline_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory= output_directory
    shell:
        "{params.base_directory}/mash_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"

rule mash_screen:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        "mash/{sample}_mashscreen.txt"
    log:
        out="logs/mash_screen/{sample}.log",
        err="logs/mash_screen/{sample}.err"
    benchmark:
        "logs/benchmark/mash_screen/{sample}.log"
    threads:
        5
    params:
        mash_sketches
    shell:
        "mash screen -p {threads} {params} -w -v 0 -i 0.95 {input.read1} {input.read2} > {output} || true ; touch {output}"

rule mash_screen_multiqc:
    input:
        expand("mash/{sample}_mashscreen.txt", sample=SAMPLE)
    output:
        "mash/mash_screen_results.txt"
    log:
        out="logs/mash_screen_multiqc/log.log",
        err="logs/mash_screen_multiqc/log/err"
    benchmark:
        "logs/benchmarks/mash_screen_multiqc/log.log"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory= output_directory
    shell:
        "{params.base_directory}/mash_screen_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"

rule prokka:
    input:
        contig_file=rules.shovill_move.output,
        mash_file=rules.mash_sort.output
    threads:
        48
    output:
        "Prokka/{sample}/{sample}.gff",
    log:
        out="logs/prokka/{sample}.log",
        err="logs/prokka/{sample}.err"
    benchmark:
        "logs/benchmark/prokka/{sample}.log"
    shell:
        """
        which prokka 2>> {log.err} | tee -a {log.out}
        prokka --version 2>> {log.err} | tee -a {log.out}
        if [ -s \"{input.mash_file}\" ]
        then
        mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | awk -F \"-.-\" '{{ print $NF }}' | sed 's/.fna//g' | awk -F \"_\" '{{ print $1 \" \" $2 }}' ))
        prokka --cpu {threads} --compliant --centre --UPHL --mincontiglen 500 --outdir Prokka/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample} --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} 2>> {log.err} | tee -a {log.out} || true
        fi
        if [ ! -f {output} ] ; then touch {output} ; fi
        """

rule prokka_move:
    input:
        rules.prokka.output
    output:
        "ALL_gff/{sample}.gff"
    log:
        out="logs/prokka_move/{sample}.log",
        err="logs/prokka_move/{sample}.err"
    benchmark:
        "logs/benchmark/prokka_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule plasmidprokka:
    input:
        contig_file=rules.plasmidshovill_move.output,
        mash_file=rules.mash_sort.output
    threads:
        48
    output:
        "Prokka_plasmids/{sample}/{sample}_plasmid.gff",
    log:
        out="logs/plasmidprokka/{sample}.log",
        err="logs/plasmidprokka/{sample}.err"
    benchmark:
        "logs/benchmark/plasmidprokka/{sample}.log"
    shell:
        """
        which prokka 2>> {log.err} | tee -a {log.out}
        prokka --version 2>> {log.err} | tee -a {log.out}
        if [ -s \"{input.mash_file}\" ]
        then
        mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | awk -F \"-.-\" '{{ print $NF }}' | sed 's/.fna//g' | awk -F \"_\" '{{ print $1 \" \" $2 }}' ))
        prokka --cpu {threads} --compliant --centre --UPHL --mincontiglen 500 --outdir Prokka_plasmids/{wildcards.sample} --locustag locus_tag --prefix {wildcards.sample}_plasmid --genus ${{mash_result[0]}} --species ${{mash_result[1]}} --force {input.contig_file} 2>> {log.err} | tee -a {log.out} || true
        fi
        if [ ! -f {output} ] ; then touch {output} ; fi
        """

rule plasmidprokka_move:
    input:
        rules.plasmidprokka.output
    output:
        "ALL_gff_plasmids/{sample}_plasmid.gff"
    log:
        out="logs/plasmidprokka_move/{sample}.log",
        err="logs/plasmidprokka_move/{sample}.err"
    benchmark:
        "logs/benchmark/plasmidprokka_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule quast:
    input:
        rules.shovill_move.output
    output:
        "quast/{sample}/report.txt",
    log:
        out="logs/quast/{sample}.log",
        err="logs/quast/{sample}.err"
    benchmark:
        "logs/benchmark/quast/{sample}.log"
    threads:
        1
    run:
        shell("which quast 2>> {log.err} | tee -a {log.out}")
        shell("quast --version 2>> {log.err} | tee -a {log.out}")
        shell("quast {input} --output-dir quast/{wildcards.sample} --threads {threads} 2>> {log.err} | tee -a {log.out} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule plasmidquast:
    input:
        rules.plasmidshovill_move.output
    output:
        "quast_plasmids/{sample}_plasmid/report.txt",
    log:
        out="logs/plasmidquast/{sample}.log",
        err="logs/plasmidquast/{sample}.err"
    benchmark:
        "logs/benchmark/plasmidquast/{sample}.log"
    threads:
        1
    run:
        shell("which quast 2>> {log.err} | tee -a {log.out}")
        shell("quast --version 2>> {log.err} | tee -a {log.out}")
        shell("quast {input} --output-dir quast_plasmids/{wildcards.sample}_plasmid --threads {threads} 2>> {log.err} | tee -a {log.out} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule CG_pipeline_shuffle_raw:
    input:
        read1= get_read1,
        read2= get_read2
    output:
        "Sequencing_reads/shuffled/{sample}_raw_shuffled.fastq.gz"
    log:
        out="logs/cg_pipeline_shuffle_raw/{sample}.log",
        err="logs/cg_pipeline_shuffle_raw/{sample}.err"
    benchmark:
        "logs/benchmark/cg_pipeline_shuffle_raw/{sample}.log"
    threads:
        1
    run:
        shell("which run_assembly_shuffleReads.pl 2>> {log.err} | tee -a {log.out}")
        shell("run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output} 2>> {log.err}")

rule CG_pipeline_shuffle_clean:
    input:
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2
    output:
        "Sequencing_reads/shuffled/{sample}_clean_shuffled.fastq.gz"
    log:
        out="logs/cg_pipeline_shuffle_clean/{sample}.log",
        err="logs/cg_pipeline_shuffle_clean/{sample}.err"
    benchmark:
        "logs/benchmark/cg_pipeline_shuffle_clean/{sample}.log"
    threads:
        1
    run:
        shell("which run_assembly_shuffleReads.pl 2>> {log.err} | tee -a {log.out}")
        shell("run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output} 2>> {log.err}")

rule CG_pipeline:
    input:
        shuffled_fastq="Sequencing_reads/shuffled/{sample}_{raw_or_clean}_shuffled.fastq.gz",
        quast_file="quast/{sample}/report.txt",
        mash_file=rules.mash_sort.output
    output:
        "cg-pipeline/{sample}.{raw_or_clean}.out.txt"
    threads:
        10
    log:
        out="logs/cg_pipeline/{sample}.{raw_or_clean}.log",
        err="logs/cg_pipeline/{sample}.{raw_or_clean}.err"
    benchmark:
        "logs/benchmark/cg_pipeline/{sample}.{raw_or_clean}.log"
    params:
        base_directory=base_directory
    run:
        shell("which run_assembly_readMetrics.pl 2>> {log.err} | tee -a {log.out}")
        if "PNUSAS" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 5000000 > {output} 2>> {log.err} || true")
        elif "PNUSAE" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 5000000 > {output} 2>> {log.err} || true")
        elif "PNUSAC" in wildcards.sample :
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 1600000 > {output} 2>> {log.err} || true")
        elif "PNUSAL" in wildcards.sample:
            shell("run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e 3000000 > {output} 2>> {log.err} || true")
        else:
            shell("{params.base_directory}/genome_length_cg.sh {wildcards.sample} {threads} {output} {input.shuffled_fastq} 2>> {log.err} | tee -a {log.out} || true")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi")


rule CG_pipeline_multiqc:
    input:
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE, raw_or_clean=['raw', 'clean'])
    output:
        "cg-pipeline/cg-pipeline-summary.txt"
    log:
        out="logs/cg_pipeline_multiqc/log.log",
        err="logs/cg_pipeline_multiqc/log.err"
    benchmark:
        "logs/benchmark/cg_pipeline_multiqc/benchmark.log"
    threads:
        1
    params:
        base_directory=base_directory,
        output_directory=output_directory
    shell:
        "{params.base_directory}/cgpipeline_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"

rule seqsero:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read2
    output:
        result="SeqSero/{sample}/Seqsero_result.txt",
        datlog="SeqSero/{sample}/data_log.txt"
    log:
        out="logs/seqsero/{sample}.log",
        err="logs/seqsero/{sample}.err"
    benchmark:
        "logs/benchmark/seqsero/{sample}.log"
    threads:
        1
    run:
        shell("which SeqSero.py 2>> {log.err} | tee -a {log.out}")
        shell("SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} 2>> {log.err} | tee -a {log.out} || true ; touch {output} ")

rule seqsero_move:
    input:
        rules.seqsero.output.result
    output:
        "SeqSero/{sample}.Seqsero_result.txt"
    log:
        out="logs/seqsero_move/{sample}.log",
        err="logs/seqsero_move/{sample}.err"
    benchmark:
        "logs/benchmark/seqsero_move/{sample}.log"
    threads:
        1
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule seqsero_multiqc:
    input:
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE)
    output:
        "SeqSero/Seqsero_serotype_results.txt"
    log:
        out="logs/seqsero_multiqc/log.log",
        err="logs/seqsero_multiqc/log.err"
    params:
        base_directory=base_directory,
        output_directory=output_directory
    threads:
        1
    shell:
        "{params.base_directory}/seqsero_multiqc.sh {params.output_directory} 2>> {log.err} | tee -a {log.out}"

rule abricate:
    input:
        rules.shovill_move.output
    output:
        "abricate_results/{database}/{database}.{sample}.out.tab"
    log:
        out="logs/abricate/{sample}.{database}.log",
        err="logs/abricate/{sample}.{database}.err"
    benchmark:
        "logs/benchmark/abricate_{database}/{sample}.log"
    threads:
        5
    run:
        shell("which abricate 2>> {log.err} | tee -a {log.out}")
        shell("abricate --version 2>> {log.err} | tee -a {log.out}")
        shell("abricate --list 2>> {log.err} | tee -a {log.out}")
        shell("abricate --db {wildcards.database} --threads {threads} {input} > {output} 2>> {log.err} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule plasmidabricate:
    input:
        rules.plasmidshovill_move.output
    output:
        "abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab"
    log:
        out="logs/plasmidabricate/{sample}.{database}.log",
        err="logs/plasmidabricate/{sample}.{database}.err"
    benchmark:
        "logs/benchmark/plasmidabricate_{database}/{sample}.log"
    threads:
        5
    run:
        shell("which abricate 2>> {log.err} | tee -a {log.out}")
        shell("abricate --version 2>> {log.err} | tee -a {log.out}")
        shell("abricate --list 2>> {log.err} | tee -a {log.out}")
        shell("abricate --db {wildcards.database} --threads {threads} {input} > {output} 2>> {log.err} || true ")
        shell("if [ ! -f {output} ] ; then touch {output} ; fi ")

rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results_plasmids/{database}/{database}.{sample}.plasmids.out.tab", sample=SAMPLE, database=DATABASE)
    output:
        "logs/abricate_results/{database}.summary.txt"
    log:
        out="logs/abricate_summary/{database}.log",
        err="logs/abricate_summary/{database}.err"
    threads:
        1
    run:
        shell("which abricate 2>> {log.err} | tee -a {log.out}")
        shell("abricate --version 2>> {log.err} | tee -a {log.out}")
        shell("abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output} 2>> {log.err}")

rule abricate_multiqc:
    input:
        rules.abricate_summary.output
    output:
        "logs/abricate_results/{database}.summary.csv"
    log:
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
        "> {output} 2>> {log.err}"

rule flash_move:
    input:
        rules.plasmidshovill.output.flash
    log:
        out="logs/flash_move/{sample}.log",
        err="logs/flash_move/{sample}.err"
    benchmark:
        "logs/benchmark/flash_move/{sample}.log"
    threads:
        1
    output:
        "results_for_multiqc/{sample}_flash.hist"
    shell:
        "cp {input} {output} 2>> {log.err} | tee -a {log.out}"

rule bwa_index:
    input:
        rules.shovill.output
    output:
        "shovill_result/{sample}/contigs.fa.sa"
    shell:
        "bwa index {input}"

rule bwa:
    input:
        contig=rules.shovill.output,
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2,
        index=rules.bwa_index.output
    threads:
        48
    output:
        bam="bwa/{sample}.sorted.bam",
        bai="bwa/{sample}.sorted.bam.bai"
    shell:
        "bwa mem -t {threads} {input.contig} {input.read1} {input.read2} | samtools sort -o {output.bam} ; "
        "samtools index {output.bam} "

rule blastn:
    input:
        rules.shovill.output
    output:
        "blast/{sample}.tsv"
    threads:
        10
    shell:
        "blastn -query {input} -out {output} -num_threads {threads} -db /home/Bioinformatics/Data/blastdb/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25"

rule blobtools_create:
    input:
        contig=rules.shovill.output,
        blast=rules.blastn.output,
        bam=rules.bwa.output.bam
    output:
        cov="blobtools/{sample}.{sample}.sorted.bam.cov",
        json="blobtools/{sample}.blobDB.json"
    threads:
        1
    shell:
        "blobtools create -o blobtools/{wildcards.sample} -i {input.contig} -b {input.bam} -t {input.blast} || true ; "
        "touch {output}"

rule blobtools_view:
    input:
        rules.blobtools_create.output.json,
    output:
        "blobtools/{sample}.blobDB.table.txt"
    threads:
        1
    shell:
        "blobtools view -i {input} -o blobtools/ || true ; "
        "touch {output}"

rule blobtools_plot:
    input:
        table=rules.blobtools_view.output,
        json=rules.blobtools_create.output.json
    output:
        "blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png",
        "blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png",
        "blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt"
    threads:
        1
    shell:
        "blobtools plot -i {input.json} -o blobtools/ -r species --format png || true ; "
        "touch {output}"
