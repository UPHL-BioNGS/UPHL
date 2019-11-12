import fnmatch
import os
import glob
import shutil
from os.path import join
print("UPHL reference free pipeline v.0.2019.11.08")
# This version removes all custom scripts

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
        # prokka results
        expand("Prokka/{sample}/{sample}.gff", sample=SAMPLE),
        expand("ALL_gff/{sample}.gff", sample=SAMPLE),
        # quast results
        expand("quast/{sample}/report.tsv", sample=SAMPLE),
        # seqsero results
        expand("SeqSero/{sample}.Seqsero_result.txt", sample=SAMPLE),
        # cg-pipeline results
        expand("cg-pipeline/{sample}.{raw_or_clean}.out.txt", sample=SAMPLE,raw_or_clean=['raw', 'clean']),
        # abricate results
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
        expand("abricate_results/summary/{database}.abricate_summary.txt", database=DATABASE),
        # blobtools results - which requires bwa, samtools, and blast before blobtools can be run
        expand("shovill_result/{sample}/contigs.fa.sa", sample=SAMPLE),
        expand("bwa/{sample}.sorted.bam", sample=SAMPLE),
        expand("bwa/{sample}.sorted.bam.bai", sample=SAMPLE),
        expand("blast/{sample}.tsv", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.json", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.table.txt", sample=SAMPLE),
        expand("blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png", sample=SAMPLE),
    output:
        log=temp("logs/all/all")
    singularity:
        "docker://ewels/multiqc:1.7"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        find quast -size 0 -delete -print     2>> {output.log}.err | xargs echo "empty quast files: "     >> {output.log}.err
        #find SeqSero -size 0 -delete -print   2>> {output.log}.err | xargs echo "empty seqsero files: "   >> {output.log}.err # TBD at later date
        #find blobtools -size 0 -delete -print 2>> {output.log}.err | xargs echo "empty blobtools files: " >> {output.log}.err # TBD at later date
        #find abricate -size 0 -delete -print  2>> {output.log}.err | xargs echo "empty abricate files: "  >> {output.log}.err # TBD at later date
        #find mash -size 0 -delete -print      2>> {output.log}.err | xargs echo "empty mash files: "      >> {output.log}.err # TBD at later date
        #find Prokka -size 0 -delete -print    2>> {output.log}.err | xargs echo "empty prokka files: "    >> {output.log}.err # TBD at later date
        # Coming soon? adding a pdf version to the muliqc rule
        multiqc --version >> logs/all/all.log
        multiqc -f \
            --outdir ./logs \
            --cl_config "prokka_fn_snames: True"  \
            ./abricate_results/summary \
            ./blobtools \
            ./cg-pipeline \
            ./fastqc \
            ./mash \
            ./Prokka \
            ./quast \
            ./SeqSero \
            ./Sequencing_reads/QCed/*tsv \
            2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "seqyclean version: $(seqyclean -h | grep Version)" >> {output.log}.log
        seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 {input.read1} -2 {input.read2} -o Sequencing_reads/QCed/{wildcards.sample}_clean 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        fastqc --version >> {output.log}.log
        fastqc --outdir fastqc --threads {threads} Sequencing_reads/*/*.fastq* 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        shovill --version >> {output.log}.log
        RAM=$(free -m --giga | grep "Mem:" | awk '{{ print ($2*0.8) }}' | cut -f 1 -d ".")
        echo "Using $RAM RAM and {threads} cpu for shovill" | tee -a {output.log}.log
        shovill --cpu {threads} --ram $RAM --outdir shovill_result/{wildcards.sample} --R1 {input.read1} --R2 {input.read2} --force 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        cp {output.file} {output.final}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "mash version: $(mash --version)" >> {output.log}.log
        cat {input.read1} {input.read2} | mash sketch -m 2 -o mash/{wildcards.sample} - 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "mash version: $(mash --version)" >> {output.log}.log
        mash dist -p {threads} -v 0 /db/RefSeqSketchesDefaults.msh {input} | sort -gk3 > {output.file} 2>> {output.log}.err || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        prokka -v >> {output.log}.log
        mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d "_" | cut -f 1 -d "." | sed 's/_/ /g' )) || mash_result=('none', 'none')
        prokka --cpu {threads} \
            --compliant \
            --centre URF \
            --mincontiglen 500 \
            --outdir Prokka/{wildcards.sample} \
            --locustag locus_tag \
            --prefix {wildcards.sample} \
            --genus ${{mash_result[0]}} \
            --species ${{mash_result[1]}} \
            --force {input.contig_file} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        cp {output.file} {output.final}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        quast.py --version >> {output.log}.log
        quast.py {input} --output-dir quast/{wildcards.sample} --threads {threads} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.log}.err || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        run_assembly_shuffleReads.pl -gz {input.read1} {input.read2} > {output.file} 2>> {output.log}.err || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        wget -nc https://raw.githubusercontent.com/StaPH-B/UPHL/master/URF_scripts/genome_sizes.txt -O logs/genome_sizes.txt 2>> {output.log}.err || true
        mash_result=($(head -n 1 {input.mash_file} | cut -f 1 | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d "_" | cut -f 1 -d "." )) || true
        genome_length=$(grep $mash_result logs/genome_sizes.txt | grep -v "#" | head -n 1 | cut -f 2 -d ":" | awk '{{ print $0 "e+06" }}' ) || \
            genome_length=$(grep 'Estimated genome size:' {input.mash_error}.err | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' | xargs printf "%.0f" ) || \
            genome_length="0"
        echo "The genome length for {wildcards.sample} is $genome_length" >> {output.log}.log
        run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e $genome_length 2>> {output.log}.err > {output.file} || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        rm -R SeqSero/{wildcards.sample} || true
        SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} #2>> {output.log}.err >> {output.log}.log || true
        cp {output.file} {output.final} || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        abricate --version >> {output.log}.log
        abricate --list >> {output.log}.log
        abricate --db {wildcards.database} --threads {threads} --minid 80 --mincov 80 {input} > {output.file} 2>> {output.log}.err || true
        touch {output}
        """

rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
    output:
        file="abricate_results/summary/{database}.abricate_summary.txt",
        log=temp("logs/abricate/{database}_summary")
    threads:
        1
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        abricate --version >> {output.log}.log
        abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output.file} 2>> {output.log}.err || true
        touch {output}
        """

rule bwa_index:
    input:
        rules.shovill.output.file
    output:
        index="shovill_result/{sample}/contigs.fa.sa",
        log=temp("logs/bwa/{sample}_index")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "bwa $(bwa 2>&1 | grep Version )" >> {output.log}.log
        bwa index {input} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        blastn -version >> {output.log}.log
        echo "The blastdb location is $BLASTDB" >> {output.log}.log
        blastn -query {input} \
            -out {output.tsv} \
            -num_threads {threads} \
            -db /blast/blastdb/nt \
            -outfmt '6 qseqid staxids bitscore std' \
            -max_target_seqs 10 \
            -max_hsps 1 \
            -evalue 1e-25 2>> {output.log}.log | tee -a {output.log}.err || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "bwa $(bwa 2>&1 | grep Version )" >> {output.log}.log
        samtools --version >> {output.log}.log
        if [ -s "{input.test}" ]
        then
            bwa mem -t {threads} {input.contig} {input.read1} {input.read2} 2>> {output.log}.err | \
            samtools sort -o {output.bam} 2>> {output.log}.err > {output.bam} || true
        else
            echo "BLASTDB not found" | tee -a {output.log}.log {output.log}.err
        fi
        samtools index {output.bam} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools create -o blobtools/{wildcards.sample} -i {input.contig} -b {input.bam} -t {input.blast} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools view -i {input} -o blobtools/ 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """

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
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools plot -i {input.json} -o blobtools/ -r species --format png 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
