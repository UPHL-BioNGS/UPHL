#!/bin/bash

output_directory=$1
ANALYSES=($(ls $output_directory/logs/benchmark/*/ -d | awk -F "/" '{ print $(NF-1) }' | grep -v "combine" | grep -v "multiqc" ))
echo ${ANALYSES[@]}
header="sample"
for analysis in ${ANALYSES[@]}
do
  header="$header,$analysis"
done
echo $header > $output_directory/logs/benchmark_summary.csv
log_files=($(ls $output_directory/logs/benchmark/*/*log | rev | cut -f 1 -d "/" | rev | grep -v "home" | sort | uniq ))
echo ${log_files[@]}
for log_file in ${log_files[@]}
do
  sample_benchmark=$log_file
  for analysis in ${ANALYSES[@]}
  do
    running_time=$(tail -n 1 $output_directory/logs/benchmark/$analysis/$log_file | cut -f 1 | sort -n -r | head -n 1 )
    if [ -z "$running_time" ] ; then running_time=0 ; fi
    sample_benchmark="$sample_benchmark,$running_time"
  done
echo "$sample_benchmark" >> $output_directory/logs/benchmark_summary.csv
done
