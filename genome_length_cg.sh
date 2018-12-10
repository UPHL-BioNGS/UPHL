#!/bin/bash

shuffled_fastq=$1
quast_file=$2
threads=$3
output=$4

genome_length=$(grep 'Total length (>= 0 bp)' $quast_file | grep -v "All statistics" | sed 's/Total length (>= 0 bp)//g' | sed 's/ //g' )
run_assembly_readMetrics.pl $shuffled_fastq --fast --numcpus $threads -e $genome_length > $output
