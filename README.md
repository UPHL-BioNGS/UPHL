![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)

# Bioinformatics at the Utah Public Health Laboratory (UPHL)

Bioinformatic analysis combines biology, computer science, mathematics, statistics, and engineering in order to analyze, interpret and make sense of the sequencing data.

This analysis allows
* For a better understanding of the genetic basis of disease
* Aids in bacterial pathogen identification
* Can ultimately assist with outbreak detection

(More information can be found at [UPHL's website](https://uphl.utah.gov/infectious-diseases/next-generation-sequencing/))

## How to download this repository:
```
git clone https://github.com/StaPH-B/UPHL.git
```
* note: to make your life easier type the following commands 
```
cd UPHL
git init
```
and then you can `git pull` any time you need to sync your version of the workflow with the one posted here, on github.

### UPHL Reference-Free Workflow
[UPHL_reference_free_docker.smk](UPHL_reference_free_docker.smk) is a snakemake workflow utilizing docker containers pulled via singularity to take Illumina NGS fastq files to annotated contigs. In order to utilize this workflow, both snakemake and singularity must be installed. Installation instructions can be found [here](READMEs/installation/README.md). 

The snakemake workflow _should_ work for short-read pair-end fastq files that either end in formats like the following: `_R1_001.fastq.gz` (Illumina) or `_1.fastq` (SRA) and all fastq files must be located in `<path to directory>/Sequencing_reads/Raw/`

Usage:
```
snakemake --snakefile UPHL_reference_free_docker.smk \          
  --directory <path to directory> \
  --use-singularity \
  --singularity-args "--bind <path to directory>:/data,<path to blast databases>:/blast/blastdb" \
  --singularity-prefix <path to where singularity images should be stored> \
  --cores 20
```
Using the `-n` option is highly encouraged before starting a run.

The majority of singularity containers are actually converted docker containers maintained by [STAPHB](https://github.com/StaPH-B/docker-builds)

A full list of commands can be found [here](READMEs/URF_commands/README.md).

### Comparing isolates and looking for outbreaks
[OUTBREAK_120.smk](OUTBREAK_120.smk) is a follow-up workflow that combines gff files from accross 120 days of sequencing runs into phylogenetic trees. Currently does not use containers, but will _soon_ in the future.
- requires [PLOTS_IQTREE.R](outbreak_120_scripts/PLOTS_IQTREE.R), [abricate_organize.sh](outbreak_120_scripts/abricate_organize.sh), [benchmark120_multiqc.sh](outbreak_120_scripts/benchmark120_multiqc.sh), [ggtree_plot_organize.sh](outbreak_120_scripts/ggtree_plot_organize.sh), [multiqc_config_outbreak120_snakemake.yaml](outbreak_120_scripts/multiqc_config_outbreak120_snakemake.yaml), [organism_multiqc.sh](outbreak_120_scripts/organism_multiqc.sh), and [roary_multiqc.sh](outbreak_120_scripts/roary_multiqc.sh)

Usage:
```
./outbreak_120_scripts/outbreak_120_organize.sh \
  -i <location of serotyping and gff files> \
  -o <path to final directory>
```
Followed by
```
snakemake --snakefile OUTBREAK_120_core_trees.smk \
  --directory <path the final directory specified above> \
  --cores 20
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
```

Admittedly, the snakefile is still fairly specific for the files and structure at UPHL. That should be changing _soon_. For others, a few lines of code will need to be adjusted:

There is one variable in the snakemake file that will need to be adjusted:
- line 9: kraken_mini_db="/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208"

Additionally, the controls with paths in [outbreak_120_organize.sh](outbreak_120_scripts/outbreak_120_organize.sh) will have to be adjusted. These are reference genomes downloaded from ncbi and run through prokka.

A full list of commands can be found [here](READMEs/OUTBREAK_120_commands/README.md).
