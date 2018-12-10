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


## B. Helpful ideas and scripts

1. Turning ABRICATE results into something multiqc-parsable with the following [BASH commands](abricate_to_multiqc_heatmap.sh) and associated [multiqc_config](abricate_multiqc_config_options.yaml) options
