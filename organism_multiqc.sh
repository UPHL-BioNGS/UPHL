#!/bin/bash
output_directory=$1
roary_summaries=($(ls $output_directory/*/*/*/Roary_out*/qc_report.csv ))
organisms=("$(cat $output_directory/*/*/*/Roary_out*/qc_report.csv | grep -v "Sample,Genus,Species" | cut -f 3 -d "," | sort | uniq | tr ' ' '_' )")
if [ ! -f "$output_directory/logs/Roary_qc_report.csv" ]
then
  header="tree"
  for organism in ${organisms[@]}
  do
    header="$header,$organism"
  done
  echo "$header" > $output_directory/logs/Roary_qc_report.csv
fi

for summary in ${roary_summaries[@]}
do
  kraken_results="$summary"
  for organism in ${organisms[@]}
  do
    hit=$(cat $summary | cut -f 3 -d "," | tr ' ' '_' | grep $organism | cut -f 3 -d "," | sort | uniq -c | awk '{ print $1 }' )
    if [ -z "$hit" ] ; then hit=0 ; fi
    kraken_results="$kraken_results,$hit"
  done
  echo "$kraken_results" >> $output_directory/logs/Roary_qc_report.csv
done
