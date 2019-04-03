#!/bin/bash
output_directory=$1
abricate_output=$2
working_directory=$3

gff_files=($(ls $output_directory/*gff ))
for gff_file in ${gff_files[@]}
do
  sample=$(echo $gff_file | sed 's!.*/!!' | sed 's/.gff//g')
  echo -e "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > $output_directory/ABRICATE/$sample.tab
  cat $working_directory/serotyping_results/abricate/*$sample*out.tab | grep -v '#' | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | sort | uniq >> $output_directory/ABRICATE/$sample.tab
done

abricate --summary $output_directory/ABRICATE/*tab > $abricate_output.orig
cat $abricate_output.orig | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }' | sed 's/\t/,/g' | sed 's/ /,/g' | sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' > $abricate_output

if [ ! -f "$abricate_output" ] ; then touch $abricate_output ; fi
