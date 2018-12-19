#!/bin/bash

PLOTS_OUTBREAK_120_IQTREE=$1
treefile=$2
roary_gene_presence=$3
roary_core_genome=$4
abricate_result=$5
nucleotide_distance=$6
roary_gene_presence_out=$7
abricate_resistance=$8
analysis_type=$9
genus=${10}
species=${11}
date=${12}
tree=${13}
bootstrap_tree=${14}
UFboot_tree=${15}
distance_tree=${16}
output=${17}
working_directory=${18}

treefile=$(echo "$working_directory/$treefile" | sed 's/\//\\\//g')
roary_gene_presence=$(echo "$working_directory/$roary_gene_presence" | sed 's/\//\\\//g')
roary_core_genome=$(echo "$working_directory/$roary_core_genome" | sed 's/\//\\\//g')
abricate_result=$(echo "$working_directory/$abricate_result" | sed 's/\//\\\//g')
nucleotide_distance=$(echo "$working_directory/$nucleotide_distance" | sed 's/\//\\\//g')
roary_gene_presence_out=$(echo "$working_directory/$roary_gene_presence_out" | sed 's/\//\\\//g')
abricate_resistance=$(echo "$working_directory/$abricate_resistance" | sed 's/\//\\\//g')
analysis_type=$(echo $analysis_type | sed 's/\//\\\//g')
genus=$(echo $genus | sed 's/\//\\\//g')
species=$(echo $species | sed 's/\//\\\//g')
date=$(echo $date | sed 's/\//\\\//g')
tree=$(echo "$working_directory/$tree" | sed 's/\//\\\//g')
bootstrap_tree=$(echo "$working_directory/$bootstrap_tree" | sed 's/\//\\\//g')
UFboot_tree=$(echo "$working_directory/$UFboot_tree" | sed 's/\//\\\//g')
distance_tree=$(echo "$working_directory/$distance_tree" | sed 's/\//\\\//g')

cat $PLOTS_OUTBREAK_120_IQTREE |
sed "s/FULLPATHTOIQTREETREEFILE/$treefile/g" |
sed "s/FULLPATHTOROARYGENETABLE/$roary_gene_presence/g" |
sed "s/FULLPATHTOROARYALIGNMENT/$roary_core_genome/g" |
sed "s/FULLPATHTOABRICATE_TABLE/$abricate_result/g" |
sed "s/FULLPATHTODISTANCEIMAGE/$nucleotide_distance/g" |
sed "s/FULLPATHTOGENETABLIMAGE/$roary_gene_presence_out/g" |
sed "s/FULLPATHTORESISTENIMAGE/$abricate_resistance/g" |
sed "s/TITLEOFTREE/$analysis_type $genus $species/g" |
sed "s/DATEOFRUN/$date/g" |
sed "s/FULLPATHTOTREEIMAGE/$tree/g" |
sed "s/FULLPATHTOBOOTSTRAPTREEIMAGE/$bootstrap_tree/g" |
sed "s/FULLPATHTO_UFBOOT_TREE_IMAGE/$UFboot_tree/g" |
sed "s/FULLPATHTODISTANCETREEIMAGE/$distance_tree/g" > $output
