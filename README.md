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

## A. The shell scripts and SNAKEMAKE files that connects everything at the Utah Public Health Laboratory
### 1. [UPHL_reference_free.sh](UPHL_reference_free.sh): UPHL's reference-free pipeline

Usage:
```
./UPHL_reference_free.sh -i /path/to/fastq -o /out/directory
```
Full list of options:
```
REQUIRED OPTIONS:
-i /path/to/fastq/files          Full path to fastq files. All fastq files must be in the same folder.
-o /output/directory             The full path of the directory where all files will be moved or created.

OPTIONAL FLAGS:
-C /path/to/seqyclean/adaptors   Path to SEQYCLEAN adaptors
-M /path/to/refseq/mash/sketches Path to RefSeqSketchesDefaults.msh
-v                               Display versions and paths of dependencies and exit
-b                               Display script commands and exit
-P <PATH>                        Add path(s) to PATH variable

PARTIAL PIPELINE FLAGS (Does not run full pipeline):
-c                               Run SEQYCLEAN
-l                               Run SHOVILL     (Requires SEQYCLEAN files)
-p                               Run PROKKA      (Requires SHOVILL contigs)
-f                               Run FASTQC      (Requires fastq files, but will also work on SEQYCLEAN files)
-q                               Run Quast       (Requires SHOVILL contigs)
-e                               Run CG-pipeline (Requires SEQYCLEAN files and QUAST results)
-m                               Run Mash        (Requires SEQYCLEAN files)
-g                               Locate Salmonella and E. coli samples from MASH results
-s                               Run SEQSERO     (Requires SEQYCLEAN files); NOTE: Salmonella samples can be listed in out/directory/salmonella.txt
-a                               Run ABRICATE    (Requires SHOVILL contigs); NOTE: E. coli samples can be listed in out/directory/ecoli.txt
-x                               Check if all files have been generated

./UPHL_reference_free.sh -o /home/NGS -m -g -s -a
```


### 2. [UPHL_reference_free.smk](UPHL_reference_free_smk): UPHL's reference-free pipeline as a set of snakemake rules 
- (requires [benchmark_multiqc.sh](benchmark_multiqc.sh), [cgpipeline_multiqc.sh](cgpipeline_multiqc.sh), [check_multiqc.sh](check_multiqc.sh), [genome_length_cg.sh](genome_length_cg.sh), [mash_multiqc.sh](mash_multiqc.sh), [multiqc_config_URF_snakemake.yaml](multiqc_config_URF_snakemake.yaml), [seqsero_multiqc.sh](seqsero_multiqc.sh), and [seqyclean_multiqc.sh](seqyclean_multiqc.sh) to be in the same folder as UPHL_reference_free.smk)

Usage:
```
snakemake --snakefile UPHL_reference_free.smk --directory <path to directory*> --cores 24 --keep-going
```
* indicates that directory contains paired-end fastq files in directory/Sequencing_reads/Raw

There are two variables that will need to be adjusted: 
- line 10: seqyclean_adaptors="/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta"
- line 11: mash_sketches="/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh"


### 3 [outbreak_120_organize.sh](outbreak_120_organize.sh) and [OUTBREAK_120_core_trees.smk](OUTBREAK_120_core_trees.smk): UPHL's method of looking through files from the last 120 days.
- requires [PLOTS_OUTBREAK_120_IQTREE.R](PLOTS_OUTBREAK_120_IQTREE.R), [abricate_organize.sh](abricate_organize.sh), [benchmark120_multiqc.sh](benchmark120_multiqc.sh), [ggtree_plot_organize.sh](ggtree_plot_organize.sh), [multiqc_config_outbreak120_snakemake.yaml](multiqc_config_outbreak120_snakemake.yaml), [organism_multiqc.sh](organism_multiqc.sh), and [roary_multiqc.sh](roary_multiqc.sh)

Usage:
```
./outbreak_120_organize.sh -i <location of serotyping and gff files> -o <path to final directory>

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

Additionally, the controls with paths in [outbreak_120_organize.sh](outbreak_120_organize.sh) will have to be adjusted. These are reference genomes downloaded from ncbi and run through prokka.


## B. Other ideas and scripts

1. Turning ABRICATE results into something multiqc-parsable with the following [BASH commands](abricate_to_multiqc_heatmap.sh) and associated [multiqc_config](abricate_multiqc_config_options.yaml) options
