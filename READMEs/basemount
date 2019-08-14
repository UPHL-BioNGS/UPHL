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
