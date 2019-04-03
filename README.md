# UPHL-Reference-Free
##### The UPHL-Reference-Free pipeline takes paired-end fastq files to contigs for microbial WGS. The pipeline utilizes the following programs, which must be included in PATH.
- [seqyclean](https://github.com/ibest/seqyclean)
- [shovill](https://github.com/tseemann/shovill)
- [prokka](https://github.com/tseemann/prokka)
- [fastqc](https://github.com/s-andrews/FastQC)
- [cg-pipeline](https://github.com/lskatz/CG-Pipeline)
- [quast](https://github.com/ablab/quast)
- [multiqc](https://github.com/ewels/MultiQC)
- [mash](https://github.com/marbl/Mash)
- [seqsero](https://github.com/denglab/SeqSero)
- [abricate](https://github.com/tseemann/abricate)

##### To turn this data into trees, the gff files generated with prokka are put through roary and iqtree with visualization done with ggtree (in R)
- [roary](https://github.com/sanger-pathogens/Roary)
- [iqtree](http://www.iqtree.org/)
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [ggtree](http://bioconductor.org/packages/release/bioc/html/ggtree.html)

## A. The SNAKEMAKE files that connects everything at the Utah Public Health Laboratory

### 1. [UPHL_reference_free.smk](UPHL_reference_free.smk): UPHL's reference-free pipeline as a set of snakemake rules 
- (requires [benchmark_multiqc.sh](URF_scripts/benchmark_multiqc.sh), [cgpipeline_multiqc.sh](URF_scripts/cgpipeline_multiqc.sh), [check_multiqc.sh](URF_scripts/check_multiqc.sh), [genome_length_cg.sh](URF_scripts/genome_length_cg.sh), [mash_multiqc.sh](URF_scripts/mash_multiqc.sh), [multiqc_config_URF_snakemake.yaml](URF_scripts/multiqc_config_URF_snakemake.yaml), [seqsero_multiqc.sh](URF_scripts/seqsero_multiqc.sh), and [seqyclean_multiqc.sh](URF_scripts/seqyclean_multiqc.sh))

Usage:
```
snakemake --snakefile UPHL_reference_free.smk --directory <path to *directory> --cores 24
```
* indicates that directory contains paired-end fastq files in *directory/Sequencing_reads/Raw

There are two variables that will need to be adjusted: 
- line 10: seqyclean_adaptors="/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta"
- line 11: mash_sketches="/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh"

### 2. [OUTBREAK_120.smk](OUTBREAK_120.smk): UPHL's method of looking through files from the last 120 days.
- requires [PLOTS_IQTREE.R](outbreak_120_scripts/PLOTS_IQTREE.R), [abricate_organize.sh](outbreak_120_scripts/abricate_organize.sh), [benchmark120_multiqc.sh](outbreak_120_scripts/benchmark120_multiqc.sh), [ggtree_plot_organize.sh](outbreak_120_scripts/ggtree_plot_organize.sh), [multiqc_config_outbreak120_snakemake.yaml](outbreak_120_scripts/multiqc_config_outbreak120_snakemake.yaml), [organism_multiqc.sh](outbreak_120_scripts/organism_multiqc.sh), and [roary_multiqc.sh](outbreak_120_scripts/roary_multiqc.sh)

Usage:
```
./outbreak_120_scripts/outbreak_120_organize.sh -i <location of serotyping and gff files> -o <path to final directory>

snakemake --snakefile OUTBREAK_120_core_trees.smk --directory <path to final directory>

```
Full list of options for outbreak_120_organize.sh
```
OUTBREAK_120 looks through a directory (or directories) for contig files.

REQUIRED OPTIONS:

-i <path to look through for input files>

OUTBREAK_120 will search this directory for files matching the following patterns:
Prokka files:     <path>/*/ALL_gff/NAME.gff
Mash files:       <path>/*/mash/NAME*.sorted.txt
SeqSero files:    <path>/*/SeqSero/NAME.Seqsero_result.txt
Abricate results: <path>/*/abricate_results/database/database.NAME.out.txt

-o <path for outfolder>

OPTIONAL OPTIONS:

-e <file>                 List of sample names to be removed from analysis if found.
-t mash/seqsero/abricate  Restrict results to specific analysis.
-p <species>              Restrict results to species in mash results.
-s <yyyy-mm-dd>           Specify end date.
-d <number of days>       Specify number of days prior to end date (default is 120)
-m <number of samples>    Adjust maximum number of samples (default is 100)
-a                        Evaluates all subfolders, not just those with samples less than 10 days old.
-P <PATH>                 Add path to PATH. Can be used to specify multiple paths.
```

There is one variable in the snakemake file that will need to be adjusted:
- line 9: kraken_mini_db="/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208"

Additionally, the controls with paths in [outbreak_120_organize.sh](outbreak_120_scripts/outbreak_120_organize.sh) will have to be adjusted. These are reference genomes downloaded from ncbi and run through prokka.


## B. Other ideas and scripts

1. Turning ABRICATE results into something multiqc-parsable with the following [BASH commands](other_ideas/abricate_to_multiqc_heatmap.sh) and associated [multiqc_config](other_ideas/abricate_multiqc_config_options.yaml) options
