import datetime
from pathlib import Path
import glob
import os, os.path
import fnmatch
#snakemake --snakefile /home/eyoung/sandbox/UPHL/UPHL_reference_free.smk --directory /home/IDGenomics_NAS/WGS_Serotyping/UT-M03999-181210 --cores 24 --rerun-incomplete
DATE=str(datetime.date.today())
working_directory=os.getcwd()
kraken_mini_db="/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208"

ANALYSIS_TYPE, GENUS, SPECIES= glob_wildcards('{ANALYSIS_TYPE}/{GENUS}/{SPECIES}/list_of_samples.txt')

ruleorder: abricate>roary>iqtree

def find_gff_files(wildcards):
    gff_files=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/*gff" )
    return(gff_files)

def find_core_genome(wildcards):
    core_genome=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/Roary_out*/core_gene_alignment.aln" )
    return(core_genome[0])

def find_gene_presence(wildcards):
    core_genome=glob.glob(wildcards.analysis_type + "/" + wildcards.genus + "/" + wildcards.species + "/Roary_out*/gene_presence_absence.Rtab" )
    return(core_genome[0])

def control_files(directory):
    control_file=glob.glob(directory + "/*control.gff")
    if len(control_file) > 0 :
        return(str("-o " + control_file[0]))
    else:
        return("")

rule all:
    input:
        # iqtree
        expand("{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.iqtree", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        # abricate
        expand("{analysis_type}/{genus}/{species}/ABRICATE/abricate_presence_absence.csv", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
        # roary
        expand("{analysis_type}/{genus}/{species}/roary.done", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES),
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
        for directory in expand("{analysis_type}.{genus}.{species}", zip, analysis_type=ANALYSIS_TYPE, genus=GENUS, species=SPECIES):
            shell("{params.base_directory}/benchmark_multiqc.sh {params.working_directory} " + directory),
        # roary results
        shell("{params.base_directory}/roary_multiqc.sh {params.working_directory}"),
        shell("{params.base_directory}/organism_multiqc.sh {params.working_directory}"),
        # running multiqc
        shell("cp {params.base_directory}/multiqc_config_outbreak120_snakemake.yaml multiqc_config.yaml"),
        shell("multiqc --outdir {params.working_directory}/logs {params.working_directory}/logs"),

rule abricate:
    input:
        find_gff_files,
    output:
        "{analysis_type}/{genus}/{species}/ABRICATE/abricate_presence_absence.csv"
    params:
        base_directory= workflow.basedir
    log:
        "logs/abricate_summary/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/abricate_summary/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    shell:
        "{params.base_directory}/abricate_organize.sh {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species} {output}"

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
        shell("roary -p {threads} -f {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/Roary_out -e -n -qc -k {params} {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/*.gff --force")
        shell("touch {output}")

rule iqtree:
    input:
        touch_roary=rules.roary.output,
        core_genome= find_core_genome,
    output:
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.ckp.gz",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.contree",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.iqtree",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.log",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.splits.nex",
        "{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.treefile",
    log:
        "logs/iqtree/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/iqtree/{analysis_type}.{genus}.{species}.log"
    threads:
        48
    run:
        control_file= control_files("{wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}")
        shell("iqtree -s {input.core_genome} -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/IQTREE/{wildcards.genus}.{wildcards.species}.iqtree -nt AUTO " + control_file)

rule ggtree_create:
    input:
        abricate_result=rules.abricate.output,
        roary_core_genome= find_core_genome,
        roary_gene_presence= find_gene_presence,
        treefile="{analysis_type}/{genus}/{species}/IQTREE/{genus}.{species}.iqtree.treefile",
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
        nucleotide_distance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.nucleotide_distance",
        roary_gene_presence="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.roary_gene_presence",
        abricate_resistance="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.abricate_resistance",
        tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.tree",
        bootstrap_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.bootstrap_tree",
        UFboot_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.UFboot_tree",
        distance_tree="{analysis_type}/{genus}/{species}/GGTREE/{analysis_type}.{genus}.{species}.distance_tree",
    shell:
        "{params.base_directory}/ggtree_plot_organize.sh "
        "{params.base_directory}/PLOTS_OUTBREAK_120_IQTREE.R "
        "{input.treefile} "
        "{input.roary_gene_presence} "
        "{input.roary_core_genome} "
        "{input.abricate_result} "
        "{params.nucleotide_distance} "
        "{params.roary_gene_presence} "
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
    log:
        "logs/ggtree_move/{analysis_type}.{genus}.{species}.log"
    benchmark:
        "logs/benchmark/ggtree_move/{analysis_type}.{genus}.{species}.log"
    threads:
        1
    shell:
        "cp {wildcards.analysis_type}/{wildcards.genus}/{wildcards.species}/GGTREE/*pdf results/."