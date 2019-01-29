import datetime
from pathlib import Path
import glob
import os, os.path
import fnmatch
print("OUTBREAK_120 v.2019.1.4")

DATE=str(datetime.date.today())
working_directory=os.getcwd()
kraken_mini_db="/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208"

ANALYSIS_TYPE, GENUS, SPECIES= glob_wildcards('{ANALYSIS_TYPE}/{GENUS}/{SPECIES}/list_of_samples.txt')
ruleorder: abricate>roary>iqtree>ggtree_create

def find_gff_files(wildcards):
    gff_files=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/*gff" )
    return(gff_files)

def find_control(wildcards):
    control_gff=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/*control.gff")
    control_file=os.path.splitext(os.path.basename(str(control_gff[0])))[0]
    if len(control_file) > 0 :
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
        expand("{analysis_type}/{genus}/{species}/GGTREE/PLOTS_OUTBREAK_120_IQTREE.R", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        expand("{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance.pdf", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        expand("results/{analysis_type}.{genus}.{species}.abricate_resistance.pdf", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
    params:
        working_directory=working_directory,
        base_directory=workflow.basedir
    threads:
        48
    run:
        # creating a table from the benchmarks
        shell("{params.base_directory}/benchmark120_multiqc.sh {params.working_directory}"),
        # roary results
        shell("{params.base_directory}/roary_multiqc.sh {params.working_directory}"),
        shell("{params.base_directory}/organism_multiqc.sh {params.working_directory}"),
        # running multiqc
        shell("which multiqc")
        shell("multiqc --version")
        shell("cp {params.base_directory}/multiqc_config_outbreak120_snakemake.yaml multiqc_config.yaml"),
        shell("multiqc --outdir {params.working_directory}/logs {params.working_directory}/logs"),

rule abricate:
    input:
        find_gff_files,
    output:
        "{analysis_type}/{genus}/{species}/ABRICATE/abricate_presence_absence.csv"
    params:
        base_directory= workflow.basedir,
        working_directory=working_directory,
    log:
        "logs/abricate_summary/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/abricate_summary/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    run:
        shell("which abricate")
        shell("abricate --version")
        shell("{params.base_directory}/abricate_organize.sh {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species} {output} {params.working_directory}")

rule roary:
    input:
        find_gff_files,
    output:
        "{analysis_type}/{genus}/{species}/roary.done"
    log:
        "logs/roary/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/roary/{analysis_type}.{genus}.{species}.log"
    params:
        kraken_mini_db
    threads:
        48
    run:
        shell("which roary")
        shell("roary -w")
        shell("if [ -d \"{wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out\" ] ; then rm -R {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out ; fi ")
        shell("roary -p {threads} -f {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out -e -n -qc -k {params} {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/*.gff --force")
        shell("touch {output}")

rule iqtree:
    input:
        touch_roary=rules.roary.output,
    output:
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.ckp.gz",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.contree",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.iqtree",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.log",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.splits.nex",
        treefile="{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.treefile",
    params:
        core_genome="{analysis_type}/{genus}/{species}/Roary_out/core_gene_alignment.aln",
        control_file=find_control
    log:
        "logs/iqtree/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/iqtree/{analysis_type}.{genus}.{species}.log"
    threads:
        48
    run:
        shell("which iqtree")
        shell("iqtree --version")
        shell("iqtree -s {params.core_genome} -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/IQTREE/{wildcards.genus}.{wildcards.species}.iqtree -nt {threads} {params.control_file} ")

rule ggtree_create:
    input:
        abricate_result=rules.abricate.output,
        roary_output=rules.roary.output,
        treefile=rules.iqtree.output.treefile,
    output:
        "{analysis_type}/{genus}/{species}/GGTREE/PLOTS_OUTBREAK_120_IQTREE.R"
    log:
        "logs/ggtree_create/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/ggtree_create/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    params:
        base_directory=workflow.basedir,
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
        "{params.base_directory}/PLOTS_OUTBREAK_120_IQTREE.R "
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
        "{params.working_directory}"

rule ggtree:
    input:
        rules.ggtree_create.output
    output:
        nucleotide_distance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.nucleotide_distance.pdf",
        roary_gene_presence="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.roary_gene_presence.pdf",
        abricate_resistance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance.pdf",
        tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.tree.pdf",
        bootstrap_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.bootstrap_tree.pdf",
        UFboot_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.UFboot_tree.pdf",
        distance_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.distance_tree.pdf",
    log:
        "logs/ggtree/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/ggtree/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    shell:
        "Rscript {input}"

rule ggtree_move:
    input:
        rules.ggtree.output.nucleotide_distance,
        rules.ggtree.output.roary_gene_presence,
        rules.ggtree.output.abricate_resistance,
    output:
        "results/{analysis_type}.{genus}.{species}.nucleotide_distance.pdf",
        "results/{analysis_type}.{genus}.{species}.roary_gene_presence.pdf",
        "results/{analysis_type}.{genus}.{species}.abricate_resistance.pdf",
        "logs/results/{analysis_type}.{genus}.{species}.nucleotide_distance_mqc.jpg",
        "logs/results/{analysis_type}.{genus}.{species}.roary_gene_presence_mqc.jpg",
        "logs/results/{analysis_type}.{genus}.{species}.abricate_resistance_mqc.jpg",
    log:
        "logs/ggtree_move/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/ggtree_move/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    run:
        shell("cp {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/GGTREE/*pdf results/.")
        shell("cp {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/GGTREE/*_mqc.jpg logs/results/.")
