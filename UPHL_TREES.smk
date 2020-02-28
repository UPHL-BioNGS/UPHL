import datetime
from pathlib import Path
import glob
import os, os.path
import fnmatch
print("OUTBREAK_120 v.0.2019.08.21")

DATE=str(datetime.date.today())
working_directory=os.getcwd()
base_directory=workflow.basedir + "/outbreak_120_scripts"
kraken_mini_db="/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208"

ANALYSIS_TYPE, GENUS, SPECIES= glob_wildcards('{ANALYSIS_TYPE}/{GENUS}/{SPECIES}/list_of_samples.txt')

def find_control(wildcards):
    if glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/*control.gff") :
        control_gff=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/*control.gff")
        control_file=os.path.splitext(os.path.basename(str(control_gff[0])))[0]
        return(str("-o " + control_file))
    else:
        return("")

rule all:
    input:
        # abricate
        expand("{analysis_type}/{genus}/{species}/ABRICATE/abricate_presence_absence.csv", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        # roary
        expand("{analysis_type}/{genus}/{species}/roary.done", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        # iqtree
        expand("{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.treefile", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        # ggtree
        expand("{analysis_type}/{genus}/{species}/GGTREE/PLOTS_IQTREE.R", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        expand("{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance.pdf", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        expand("results/{analysis_type}.{genus}.{species}.abricate_resistance.pdf", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
    params:
        working_directory=working_directory,
        base_directory=base_directory
    threads:
        48
    singularity:
        "docker://staphb/multiqc:1.8"
    shell:
        """
        which multiqc 2>> {output.log}.err | tee -a {output.log}.log
        multiqc --version 2>> {output.log}.err | tee -a {output.log}.log
        multiqc --outdir logs . 2>> {output.log}.err | tee -a {output.log}.log || true
        """

rule abricate:
    input:
        "{analysis_type}/{genus}/{species}/list_of_samples.txt",
    output:
        file="{analysis_type}/{genus}/{species}/ABRICATE/abricate_presence_absence.csv",
        log=temp("logs/abricate_summary/{analysis_type}.{genus}.{species}")
    threads:
        1
    singularity:
        "docker:/staphb/abricate:0.8.13s"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "abricate --version >> {output.log}.log ; " # version of abricate
        "gff_files=($(ls $output_directory/*gff )) ; "
        """
        for gff_file in ${gff_files[@]}
        do
          sample=$(echo $gff_file | sed 's!.*/!!' | sed 's/.gff//g')
          echo -e "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/ABRICATE/$sample.tab
          cat serotyping_results/abricate/*$sample*out.tab | grep -v '#' | awk '{{ if ($10 > 80) print $0 }}' | awk '{{ if ($9 > 80) print $0 }}' | sort | uniq >> {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/ABRICATE/$sample.tab
        done
        """
        "abricate --summary {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/ABRICATE/*tab > {output.file} 2>> {output.log}.err "
        "|| true ; touch {output}"

rule roary:
    input:
        "{analysis_type}/{genus}/{species}/list_of_samples.txt",
    output:
        file="{analysis_type}/{genus}/{species}/roary.done",
        log=temp("logs/roary/{analysis_type}.{genus}.{species}")
    singularity:
        "docker:/staphb/roary:3.12.0"
    threads:
        48
    shell:
        "date >> {output.log}.log ; " # time stamp
        "roary -w 2>> {output.log}.err | tee -a {output.log}.log ; " # version of roary
        "if [ -d \"{wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out\" ] ; then rm -R {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out ; fi  2>> {output.log}.err | tee -a {output.log}.log ; "
        "roary -p {threads} -f {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out -e -n -qc -k /kraken_db {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/*.gff --force 2>> {output.log}.err | tee -a {output.log}.log || true ; "
        "touch {output}"

rule iqtree:
    input:
        file=rules.roary.output.file,
        control_file=find_control
    output:
        file="{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.treefile",
        log=temp("logs/iqtree/{analysis_type}.{genus}.{species}")
    threads:
        48
    singularity:
        "docker:/staphb/iqtree:1.6.7"
    shell:
        "date >> {output.log}.log ; " # time stamp
        "iqtree --version 2>> {output.log}.err | tee -a {output.log}.log ; "
        "if [ -d \"{wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/IQTREE\" ] ; then rm -R {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/IQTREE/* ; fi  2>> {output.log}.err | tee -a {output.log}.log ; "
        "iqtree -s {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out/core_gene_alignment.aln -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/IQTREE/{wildcards.genus}.{wildcards.species}.iqtree -nt {threads} {input.control_file}  2>> {output.log}.err | tee -a {output.log}.log "
        "|| true ; touch {output}")

rule ggtree_create:
    input:
        abricate_result=rules.abricate.output.file,
        roary_output=rules.roary.output.file,
        treefile=rules.iqtree.output.file,
    output:
        file="{analysis_type}/{genus}/{species}/GGTREE/PLOTS_IQTREE.R",
        log=temp("logs/ggtree/{analysis_type}.{genus}.{species}_create"),
    threads:
        1
    params:
        base_directory=base_directory,
        working_directory=working_directory,
        date=DATE,
        roary_core_genome="{analysis_type}/{genus}/{species}/Roary_out/core_gene_alignment.aln",
        roary_gene_presence="{analysis_type}/{genus}/{species}/Roary_out/gene_presence_absence.Rtab",
        nucleotide_distance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.nucleotide_distance",
        roary_gene_presence_out="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.roary_gene_presence",
        abricate_resistance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance",
        tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.tree",
        bootstrap_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.bootstrap_tree",
        UFboot_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.UFboot_tree",
        distance_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.distance_tree",
    shell:
        "{params.base_directory}/ggtree_plot_organize.sh "
        "{params.base_directory}/PLOTS_IQTREE.R "
        "{input.treefile} "
        "{params.roary_gene_presence} "
        "{params.roary_core_genome} "
        "{input.abricate_result} "
        "{params.nucleotide_distance} "
        "{params.roary_gene_presence_out} "
        "{params.abricate_resistance} "
        "{wildcards.analysis_type} "
        "{wildcards.genus} "
        "{wildcards.species} "
        "{params.date} "
        "{params.tree} "
        "{params.bootstrap_tree} "
        "{params.UFboot_tree} "
        "{params.distance_tree} "
        "{output} "
        "{params.working_directory} "
        "2>> {output.log}.err | tee -a {output.log}.log"

rule ggtree:
    input:
        rules.ggtree_create.output
    output:
        nucleotide_distance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.nucleotide_distance.pdf",
        roary_gene_presence="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.roary_gene_presence.pdf",
        abricate_resistance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance.pdf",
        tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.tree.pdf",
        bootstrap_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.bootstrap_tree.pdf",
        log=temp("logs/ggtree/{analysis_type}.{genus}.{species}"),
    threads:
        1
    shell:
        "Rscript {input} 2>> {output.log}.err | tee -a {output.log}.log || true ; touch {output}"

rule ggtree_move:
    input:
        rules.ggtree.output.nucleotide_distance,
        rules.ggtree.output.roary_gene_presence,
        rules.ggtree.output.abricate_resistance,
    output:
        pdf_dist="results/{analysis_type}.{genus}.{species}.nucleotide_distance.pdf",
        pdf_gene="results/{analysis_type}.{genus}.{species}.roary_gene_presence.pdf",
        pdf_abri="results/{analysis_type}.{genus}.{species}.abricate_resistance.pdf",
        jpg_dist="logs/results/{analysis_type}.{genus}.{species}.nucleotide_distance_mqc.jpg",
        jpg_gene="logs/results/{analysis_type}.{genus}.{species}.roary_gene_presence_mqc.jpg",
        jpg_abri="logs/results/{analysis_type}.{genus}.{species}.abricate_resistance_mqc.jpg",
        log=temp("logs/ggtree/{analysis_type}.{genus}.{species}_move"),
    threads:
        1
    run:
        shell("date >> {output.log}.log") # time stamp
        shell("cp {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/GGTREE/*pdf results/. 2>> {output.log}.err | tee -a {output.log}.log || true ; touch {output.pdf_dist} {output.pdf_gene} {output.pdf_abri}")
        shell("cp {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/GGTREE/*_mqc.jpg logs/results/. 2>> {output.log}.err | tee -a {output.log}.log || true ; touch {output.jpg_dist} {output.jpg_gene} {output.jpg_abri}")
