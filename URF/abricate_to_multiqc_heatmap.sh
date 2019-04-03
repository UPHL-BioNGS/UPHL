#!/bin/bash

# optional steps: creating one file for all samples, and then separating it into individual files.
abricate -db $database $out/ALL_assembled/*_contigs.fa > $out/abricate_results/$database.out.txt
ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | parallel "echo -e '#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT' > $out/abricate_results/$database/$database.{}.tab "
ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | parallel "grep {} $out/abricate_results/$database.out.txt >> $out/abricate_results/$database/$database.{}.tab"

# abricate can also be done on individual samples:
ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | parallel "abricate -db $database $out/ALL_assembled/{}_contigs.fa > $out/abricate_results/$database.{}.tab"

# create the abricate heatmaps for multiqc
abricate --summary $out/abricate_results/$database/$database*tab > $out/abricate_results/$database/$database.summary.txt

# making the file parsible for multiqc
cat $out/abricate_results/$database/$database.summary.txt \
| sed 's/#//g' \ # removing the '#' at the beginning of the header line
| sed 's/.tab//g' \ # removing the '.tab' at the end of each file name
| sed "s/$database.//g" \ # fixing the sample name
| awk '{ sub("^.*/", "", $1); print}' \ # fixing the sample name
| awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' \ # some abricate results are split in two values. Only the first is used in the heatmap
| awk '{ $2="" ; print $0 }' \ # the second column is the total number of hits, so this column is removed
| sed 's/\t/,/g' \ # getting rid of tabs and spaces and converting to a csv format
| sed 's/ /,/g' \
| sed 's/[.],/0,/g' \ # changing empty values to 0
| sed 's/,[.]/,0/g' \
| sed 's/,,/,/g' > $out/abricate_results/$database/$database.summary.csv
