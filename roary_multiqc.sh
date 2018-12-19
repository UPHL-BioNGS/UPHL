#!/bin/bash
output_directory=$1
roary_summaries=($(ls $output_directory/*/*/*/Roary_out/summary_statistics.txt ))
if [ ! -f "$output_directory/logs/Roary_summary_statistics.txt" ]
then
  echo -e "tree\tcore_genes\tsoft_core_genes\tshell_genes\tcloud_genes" > $output_directory/logs/Roary_summary_statistics.txt
fi

for summary in ${roary_summaries[@]}
do
  core_genes=$(grep "Core genes" $summary | cut -f 3 )
  soft_core_genes=$(grep "Soft core genes" $summary | cut -f 3 )
  shell_genes=$(grep "Shell genes" $summary | cut -f 3 )
  cloud_genes=$(grep "Cloud genes" $summary | cut -f 3 )
  echo -e "$summary\t$core_genes\t$soft_core_genes\t$shell_genes\t$cloud_genes" >> $output_directory/logs/Roary_summary_statistics.txt
done
