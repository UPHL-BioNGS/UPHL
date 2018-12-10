#!/bin/bash

output_directory=$1
SAMPLES=($(ls $output_directory/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq ))
ANALYSES=($(ls $output_directory/logs/benchmark/*/ -d | awk -F "/" '{ print $(NF-1) }' | grep -v "combine" | grep -v "multiqc" ))

header="sample"
for analysis in ${ANALYSES[@]}
do
  header="$header,$analysis"
done
echo $header > $output_directory/logs/benchmark_summary.csv

for sample in ${SAMPLES[@]}
do
  sample_benchmark="$sample"
  for analysis in ${ANALYSES[@]}
  do
    running_time=$(tail -n 1 $output_directory/logs/benchmark/$analysis/$sample*.log | cut -f 1 | sort -n -r | head -n 1 )
    if [ -z "$running_time" ]
    then
      running_time=0
    fi
    sample_benchmark="$sample_benchmark,$running_time"
  done
  echo "$sample_benchmark" >> $output_directory/logs/benchmark_summary.csv
done
