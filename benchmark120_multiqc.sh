#!/bin/bash

output_directory=$1
analysis_type_genus_species=$2

ANALYSES=($(ls $output_directory/logs/benchmark/*/ -d | awk -F "/" '{ print $(NF-1) }' | grep -v "combine" | grep -v "multiqc" ))

if [ ! -f "$output_directory/logs/benchmark_summary.csv" ]
then
  header="sample"
  for analysis in ${ANALYSES[@]}
  do
    header="$header,$analysis"
  done
  echo $header > $output_directory/logs/benchmark_summary.csv
fi

sample_benchmark="$analysis_type_genus_species"
for analysis in ${ANALYSES[@]}
do
  running_time=$(tail -n 1 $output_directory/logs/benchmark/$analysis/$analysis_type_genus_species*.log | cut -f 1 | sort -n -r | head -n 1 )
  if [ -z "$running_time" ] ; then running_time=0 ; fi
  sample_benchmark="$sample_benchmark,$running_time"
done
echo "$sample_benchmark" >> $output_directory/logs/benchmark_summary.csv
