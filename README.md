![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)

# Bioinformatics at the Utah Public Health Laboratory (UPHL)

Bioinformatic analysis combines biology, computer science, mathematics, statistics, and engineering in order to analyze, interpret and make sense of the sequencing data.

This analysis allows
* For a better understanding of the genetic basis of disease
* Aids in bacterial pathogen identification
* Can ultimately assist with outbreak detection

(More information can be found at [UPHL's website](https://uphl.utah.gov/infectious-diseases/next-generation-sequencing/))

# The UPHL-Reference-Free pipeline takes paired-end fastq files to contigs for microbial WGS. The pipeline utilizes the following programs.
- [seqyclean](https://github.com/ibest/seqyclean)  
```
seqyclean -minlen 25 -qual -c Adapters.fasta -1 sample_raw1.fastq -2 sample_raw2.fastq -o sample_clean
```
- [shovill](https://github.com/tseemann/shovill)
```
shovill --cpu 48 --ram 200 --outdir shovill_result/sample --R1 sample_cleanPE1.fastq --R2 sample_cleanPE2.fastq
```
- [prokka](https://github.com/tseemann/prokka)
```
prokka --cpu 48 --compliant --centre --UPHL --mincontiglen 500 --outdir Prokka/sample --locustag locus_tag --prefix sample --genus $genus --species $species --force shovill_result/sample/contigs.fa 
```
- [fastqc](https://github.com/s-andrews/FastQC)
```
fastqc *fastq*
```
- [cg-pipeline](https://github.com/lskatz/CG-Pipeline)
```
run_assembly_shuffleReads.pl -gz sample_raw1.fastq sample_raw2.fastq > sample_shuffled.fastq
run_assembly_readMetrics.pl sample_shuffled.fastq --fast --numcpus 48 -e $genome_size
```
- [quast](https://github.com/ablab/quast)
```
quast shovill_result/sample/contigs.fa --output-dir quast/sample
```
- [multiqc](https://github.com/ewels/MultiQC)
```
multiqc .
```
- [mash](https://github.com/marbl/Mash)
```
cat sample_cleanPE1.fastq sample_cleanPE2.fastq sample_cleanSE.fastq > sample_all.fastq
mash sketch -m 2 sample_all.fastq
mash dist RefSeqSketchesDefaults.msh sample_all.fastq.msh > sample_all.fastq.msh.distance.txt
```
- [seqsero](https://github.com/denglab/SeqSero)
```
SeqSero.py -m 2 -d SeqSero/sample -i sample_cleanPE1.fastq sample_cleanPE2.fastq
```
- [abricate](https://github.com/tseemann/abricate)
```
abricate --db database --threads 48 shovill_result/sample/contigs.fa > database.sample.out.tab
abricate --summary database*out.tab > database.summary.txt
```

##### To turn this data into trees, the gff files generated with prokka are put through roary and iqtree with visualization done with ggtree (in R)
- [roary](https://github.com/sanger-pathogens/Roary)
```
roary -p 48 -f Roary_out -e -n -qc -k kraken_mini_db/minikraken_20141208 Prokka/*/*.gff
```
- [iqtree](http://www.iqtree.org/)
```
iqtree -s Roary_out/core_gene_alignment.aln -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre sample.iqtree -nt 48 -o $control
```
- [ape](https://cran.r-project.org/web/packages/ape/index.html) & [ggtree](http://bioconductor.org/packages/release/bioc/html/ggtree.html) (For more information, see [PLOTS_IQTREE.R](outbreak_120_scripts/PLOTS_IQTREE.R))
```
library(ggplot2)
library(ggtree)
library(treeio)
library(seqinr)
library(ape)
library(ade4)
library(gplots)
library(ggstance)
library(phytools)
library(viridis)
```

# Automatic detection of new sequencing run completion and pipeline initiation. 

Most of our sequencing is done on the Illumina MiSeq. Fastq generation is done through BaseSpace and seemless is transfer is made possible via [basemount](https://basemount.basespace.illumina.com/). 

Mounting BaseSpace to our linux workstation:

```
basemount /home/BaseSpace
```  

Ideally, we look for new runs as they are completed and begin our pipeline. Basemount as two main folders for our purposes:
`BaseSpace/Runs`
and
`BaseSpace/Projects`

Although our script to do this is specific for our file structure, the essential principles can be used for any bash script:
#### Part 1 : Finding a new run on basemount

```
while [ -d "/home/BaseSpace" ] # This directory _should_ always exist, so this keeps the while loop open
do
  basespace_runs=($(ls /home/BaseSpace/Runs/*/ -d | rev | cut -f 1 -d "/" | rev )) # lists all the runs in basespace
  for sequencing_run in ${basespace_runs[@]}
  do
    if [ ! -d "/WGS_DIRECTORY/$sequencing_run" ]
    then
      echo "New sequencing run detected: $sequencing_run"
      script_to_move_files_and_start_pipeline.sh "$sequencing_run"
    else 
      echo "No new sequencing run detected"
    fi
  done

  sleep 10m # waits for 10 minutes
  if [ -n "$(basemount-cmd --path /home/BaseSpace/Runs refresh | grep Error)" ] # refreshes basemount and if there's an error basemount is restarted
  then
    yes | basemount --unmount /home/BaseSpace
    basemount /home/BaseSpace
  fi

done
```
#### Part 2 : Copying files from basemount and starting the pipeline

_script_to_move_files_and_start_pipeline.sh_ probably contains a chunk of code like the following:

```
sequencing_run=$1

test_file=$(timeout -k 2m 1m find /home/BaseSpace/Projects/$sequencing_run/Samples/ -iname *fastq.gz | head -n 1 )
while [ -z "$test_file" ]
do
  echo "Fastq files are not located in /home/BaseSpace/Projects/$sequencing_run/Samples/"
  if [ -n "$(basemount-cmd --path /home/BaseSpace/Projects refresh | grep Error)" ]
  then
    yes | basemount --unmount /home/BaseSpace
    basemount /home/BaseSpace
  fi
  sleep 20m
  test_file=$(timeout -k 2m 1m find /home/BaseSpace/Projects/$run/Samples/ -iname *$run*fastq.gz | head -n 1)
done

mkdir -p /WGS_DIRECTORY/$sequencing_run/Sequencing_reads/Raw

echo "Copying files from /home/BaseSpace/Projects/$sequencing_run/Samples/*/Files/*fastq.gz to /WGS_DIRECTORY/$sequencing_run/Sequencing_reads/Raw/"
cp /home/BaseSpace/Projects/$sequencing_run/Samples/*/Files/*$sequencing_run*fastq.gz /WGS_DIRECTORY/$sequencing_run/Sequencing_reads/Raw/.

echo "Now starting snakemake"
snakemake --snakefile /home/Bioinformatics/UPHL/UPHL_reference_free.smk --directory /WGS_DIRECTORY/$sequencing_run --cores 22 

```

Now all of the fastq files are in `/WGS_DIRECTORY/$sequencing_run/Sequencing_reads/Raw/`

The snakefile will work for short-read pair-end fastq files that either end in formats like the following: `_R1_001.fastq.gz` (Illumina) or `_1.fastq` (SRA)

## SNAKEMAKE WORKFLOWS for UPHL's pipelines

### [UPHL_reference_free.smk](UPHL_reference_free.smk): UPHL's reference-free pipeline as a set of snakemake rules 
- requires internal scripts [benchmark_multiqc.sh](URF_scripts/benchmark_multiqc.sh), [cgpipeline_multiqc.sh](URF_scripts/cgpipeline_multiqc.sh), [check_multiqc.sh](URF_scripts/check_multiqc.sh), [genome_length_cg.sh](URF_scripts/genome_length_cg.sh), [mash_multiqc.sh](URF_scripts/mash_multiqc.sh), [multiqc_config_URF_snakemake.yaml](URF_scripts/multiqc_config_URF_snakemake.yaml), [seqsero_multiqc.sh](URF_scripts/seqsero_multiqc.sh), and [seqyclean_multiqc.sh](URF_scripts/seqyclean_multiqc.sh)
- requires external commands grep, touch, sed, awk, cat, find, parallel, multiqc, seqyclean, abricate, fastqc, shovill, mash, prokka, quast, cg-pipeline's run_assemply_shuffleReads.pl and run_assembly_readMetrics.pl, SeqSero.py

Usage:
```
snakemake --snakefile UPHL_reference_free.smk --directory <path to *directory> --cores 24
```
* indicates that directory contains paired-end fastq files in *directory/Sequencing_reads/Raw

There are two variables that will need to be adjusted: 
- line 10: seqyclean_adaptors="/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta"
- line 11: mash_sketches="/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh"

### __OR__
(since we are aware that using cloud computing is a popular option)

### [UPHL_reference_free_docker.smk](UPHL_reference_free_docker.smk): UPHL's reference-free pipeline as a set of snakemake rules that are compatible with singularity images
- (requires internal scripts [benchmark_multiqc.sh](URF_scripts/benchmark_multiqc.sh), [check_multiqc.sh](URF_scripts/check_multiqc.sh), [genome_length_cg.sh](URF_scripts/genome_length_cg.sh), [mash_multiqc.sh](URF_scripts/mash_multiqc.sh), [multiqc_config_URF_snakemake.yaml](URF_scripts/multiqc_config_URF_snakemake.yaml), and [seqsero_multiqc.sh](URF_scripts/seqsero_multiqc.sh))
- requires external commands singularity, multiqc, grep, touch, sed, awk, cat, find, and parallel

Usage:
```
snakemake --snakefile UPHL_reference_free_docker.smk --directory <path to *directory> --use-singularity --singularity-args "--bind <path to *directory>:/data" --cores 20
```

The majority of singularity containers are actually converted docker containers maintained by [STAPHB](https://github.com/StaPH-B/docker-builds)

## Organizing recent organisms for outbreak and cluster detection:

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
