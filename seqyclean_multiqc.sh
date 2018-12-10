#!/bin/bash

out=$1

grep "Version" $out/Sequencing_reads/Logs/*clean_SummaryStatistics.tsv | cut -f 2- -d ':' | head -n 1 | cut -f 2- | awk '{ print "Sample," $0 }' | sed 's/PE1PE2/PE1,PE2/g' | sed 's/ /,/g' | sed 's/\t/,/g' | sed 's/,,/,/g'  > $out/Sequencing_reads/Logs/seqyclean_summary.txt
grep "Sequencing_reads" $out/Sequencing_reads/Logs/*_clean_SummaryStatistics.tsv | cut -f 2- -d ':' | sort | uniq | cut -f 2- | awk '{ print $1 " " $0 }' | awk '{ sub("^.*/", "", $1); print $0 }' | awk '{ sub("_.*fastq.*,", "", $1 ); print $0 }' | sed 's/ /,/g' | sed 's/\t/,/g' | sed 's/,,/,/g' >> $out/Sequencing_reads/Logs/seqyclean_summary.txt
