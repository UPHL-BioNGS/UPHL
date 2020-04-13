print("UPHL ONT Pipeline v.20200406")
# conda activate ivar
# awk -F $'\t' 'BEGIN{OFS=FS;}{$5=60;print}' primer_schemes/nCoV-2019/V3/nCoV-2019.bed > primer_schemes/nCoV-2019/V3/nCoV-2019_col5_replaced.bed
# bwa index /home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta
# snakemake --snakefile ~/sandbox/COVID19/Illumina_covid_V3.smk --cores 20 -p

amplicon_bed="/home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/V3/V3_amplicons.bed"
reference="/home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
features="/home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/V3/GCF_009858895.2_ASM985889v3_genomic.gff"
# There are now two primer bed files... I'm thinking this is the one ivar will accept
primer_bed="/home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019_col5_replaced.bed"

SAMPLE=[]
with open('covid_samples.txt', 'r') as file:
    SAMPLE=file.read().splitlines()
print(SAMPLE)

rule all:
    input:
        expand("covid/bwa/{sample}.sorted.bam", sample=SAMPLE),
        expand("covid/consensus/{sample}.consensus", sample=SAMPLE),
        expand("covid/variants/{sample}.variants", sample=SAMPLE),
        expand("covid/quast/{sample}/report.tsv", sample=SAMPLE),
        expand("covid/samtools_stats/bwa/{sample}.txt", sample=SAMPLE),
        expand("covid/samtools_coverage/bwa/{sample}.txt", sample=SAMPLE),
        "covid/bedtools/multicov.txt",
        "covid/summary.txt"
    shell:
        """
        /home/linuxbrew/.linuxbrew/bin/multiqc -f --outdir covid covid --dirs
        """

rule bwa:
    input:
        read1="Sequencing_reads/QCed/{sample}_clean_PE1.fastq",
        read2="Sequencing_reads/QCed/{sample}_clean_PE2.fastq"
    threads:
        48
    output:
        bam="covid/bwa/{sample}.sorted.bam",
        bai="covid/bwa/{sample}.sorted.bam.bai",
        log=temp("logs/bwa_covid/{sample}")
    params:
        reference=reference
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "bwa $(bwa 2>&1 | grep Version )" >> {output.log}.log
        samtools --version >> {output.log}.log
        bwa mem -t {threads} {params.reference} {input.read1} {input.read2} 2>> {output.log}.err | \
            samtools sort 2>> {output.log}.err | \
            samtools view -F 4 -o {output.bam}

        samtools index {output.bam} 2>> {output.log}.err | tee -a {output.log}.log
        touch {output}
        """

# note: conda activate ivar
rule ivar_trim:
    input:
        rules.bwa.output.bam
    output:
        bam="covid/trimmed/{sample}.primertrim",
        log=temp("logs/ivar_trim/{sample}")
    params:
        primer_bed=primer_bed
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        ivar --version >> {output.log}.log
        ivar trim -e -i {input} -b {params.primer_bed} -p {output.bam}
        touch {output}
        """

rule samtools_sort:
    input:
        rules.ivar_trim.output.bam
    output:
        bam="covid/sorted/{sample}.primertrim.sorted.bam",
        log=temp("logs/samtools_sort/{sample}")
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        samtools sort {input}.bam -o {output.bam}
        samtools index {output.bam}
        touch {output}
        """

rule ivar_variants:
    input:
        rules.samtools_sort.output.bam
    output:
        tsv="covid/variants/{sample}.variants",
        log=temp("logs/ivar_variants/{sample}")
    params:
        reference=reference,
        features=features
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        samtools mpileup -A -d 600000 -B -Q 0 --reference {params.reference} {input} | \
            ivar variants -p {output.tsv} -q 20 -t 0.6 -r {params.reference} -g {params.features}
        touch {output}
        """

rule ivar_consensus:
    input:
        rules.samtools_sort.output.bam
    output:
        fasta="covid/consensus/{sample}.consensus",
        log=temp("logs/ivar_consensus/{sample}")
    params:
        reference=reference
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        samtools mpileup -A -d 6000000 -B -Q 0 --reference {params.reference} {input} | \
            ivar consensus -t 0.6 -p {output.fasta} -n N
        touch {output}
        """

rule quast:
    input:
        consensus=rules.ivar_consensus.output.fasta,
        bam=rules.bwa.output.bam
    output:
        file="covid/quast/{sample}/report.tsv",
        log=temp("logs/quast_covid/{sample}")
    threads:
        1
    params:
        reference=reference,
        features=features
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        quast {input.consensus}.fa \
            -r {params.reference} \
            --features {params.features} \
            --ref-bam {input.bam} \
            --output-dir covid/quast/{wildcards.sample}
        touch {output}
        """

rule samtools_stats:
    input:
        bwa=rules.bwa.output.bam,
        sort=rules.samtools_sort.output.bam,
    output:
        bwa="covid/samtools_stats/bwa/{sample}.txt",
        sort="covid/samtools_stats/sort/{sample}.txt",
        log=temp("logs/samtools_stats/{sample}")
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        samtools stats {input.bwa} > {output.bwa}
        samtools stats {input.sort} > {output.sort}
        touch {output}
        """

rule samtools_coverage:
    input:
        bwa=rules.bwa.output.bam,
        sort=rules.samtools_sort.output.bam,
    output:
        bwa="covid/samtools_coverage/bwa/{sample}.txt",
        sort="covid/samtools_coverage/sort/{sample}.txt",
    shell:
        """
        samtools coverage {input.bwa} -m -o {output.bwa}.hist
        samtools coverage {input.bwa} -o {output.bwa}
        samtools coverage {input.sort} -m -o {output.sort}.hist
        samtools coverage {input.sort} -o {output.sort}
        """

rule bedtools:
    input:
        bwa=expand("covid/bwa/{sample}.sorted.bam", sample=SAMPLE),
        sort=expand("covid/sorted/{sample}.primertrim.sorted.bam", sample=SAMPLE),
        #expand("covid/samtools_stats/trim/{sample}.txt", sample=SAMPLE),
    output:
        txt="covid/bedtools/multicov.txt",
        log=temp("logs/bedtools/multicov")
    params:
        amplicon_bed=amplicon_bed
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "primer" $(ls covid/bwa/*bam covid/sorted/*bam) | tr ' ' '\\t' > {output.txt}
        bedtools multicov -bams {input.bwa} {input.sort} -bed {params.amplicon_bed} | cut -f 4,6- >> {output.txt}
        touch {output}
        """

rule summary:
    input:
        expand("covid/bwa/{sample}.sorted.bam", sample=SAMPLE),
        expand("covid/consensus/{sample}.consensus", sample=SAMPLE),
        expand("covid/quast/{sample}/report.tsv", sample=SAMPLE),
        expand("covid/samtools_stats/bwa/{sample}.txt", sample=SAMPLE),
        expand("covid/samtools_coverage/bwa/{sample}.txt", sample=SAMPLE),
        "covid/bedtools/multicov.txt",
    output:
        "covid/summary.txt"
    shell:
        """
        date | tee -a {output}
        echo "sample,%_Human_reads,degenerate_check,coverage,depth,failed_amplicons,num_N" | tee -a {output}
        while read line
        do
            human_reads=$(grep "Homo" blobtools/$line*txt | cut -f 13 ) || human_reads="none"
            degenerate=$(grep -f ~/degenerate.txt covid/consensus/$line*fa | grep -v ">" | wc -l ) || degenerate="none"
            cov_and_depth=($(cut -f 6,7 covid/samtools_coverage/bwa/$line*txt | tail -n 1)) || cov_and_depth=(0 0)
            bedtools_column=$(head -n 1 covid/bedtools/multicov.txt | tr '\t' '\n' | grep -n $line | grep -v primertrim | cut -f 1 -d ":" | head -n 1 )
            amp_fail=$(cat covid/bedtools/multicov.txt | cut -f $bedtools_column | awk '{{ if ( $1 < 20 ) print $0 }}' | wc -l ) || amp_fail=0
            num_of_N=$(grep -o 'N' covid/consensus/$line*fa | wc -l ) || num_of_N=0
            echo "$line,$human_reads,$degenerate,${{cov_and_depth[0]}},${{cov_and_depth[1]}},$amp_fail,$num_of_N" | tee -a {output}
        done < covid_samples.txt
#        ~/sandbox/COVID19/files_for_submission.sh $(pwd)
        """
