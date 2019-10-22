#!/bin/bash
out=$1
SAMPLES=($(ls $out/Sequencing_reads/Raw/*fastq* | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | cut -f 1 -d "." | sort | uniq ))

#find $out -type f -empty | grep -v final.txt | grep -v fastqc.complete
#wc -l abricate_results*/*/*out.tab | awk '{ if ( $1 == 1 ) print $2 }' | parallel " rm $out/{} ; touch $out/{} "
#grep -P "Predicted antigenic profile:\t-:-:-" SeqSero/*Seqsero_result.txt | cut -f 1 -d ":" | parallel " rm $out/{} ; touch $out/{} "

# fastqc results
FASTQC_RESULTS=()
echo "sample,x,y" > $out/logs/raw_clean_scatter.csv
for sample in ${SAMPLES[@]}
do
  if [ -n "$(find fastqc -iname $sample*zip | grep -v shuffled | grep -v clean )" ]
  then
    file_raw1=$(ls fastqc/$sample*zip | grep -v "clean" | rev | cut -f 1 -d "/" | rev | head -n 1)
    file_raw2=$(ls fastqc/$sample*zip | grep -v "clean" | rev | cut -f 1 -d "/" | rev | tail -n 1)
  else
    file_raw1="not_found" ; file_raw2="not_found"
  fi

  if [ -n "$(find fastqc -iname $sample*zip | grep clean_PE1 )" ]
  then
    file_cln1=$(ls fastqc/$sample*zip | grep "clean_PE1" | rev | cut -f 1 -d "/" | rev | head -n 1)
    file_cln2=$(ls fastqc/$sample*zip | grep "clean_PE2" | rev | cut -f 1 -d "/" | rev | tail -n 1)
  else
    file_cln1="not_found" ; file_cln2="not_found"
  fi

  if [ -s "fastqc/$file_raw1" ]
  then
    fastqc_summry=$(unzip -l fastqc/$file_raw1 | grep fastqc_data.txt | awk '{ print $4 }' )
    result_raw1=$(unzip -p fastqc/$file_raw1 $fastqc_summry | grep "Total Sequences" | awk '{ print $3 }' )
  else
    file_raw1="not_found"
    result_raw1="not_found"
  fi
  if [ -s "fastqc/$file_raw2" ]
  then
    fastqc_summry=$(unzip -l fastqc/$file_raw2 | grep fastqc_data.txt | awk '{ print $4 }' )
    result_raw2=$(unzip -p fastqc/$file_raw2 $fastqc_summry | grep "Total Sequences" | awk '{ print $3 }' )
  else
    file_raw2="not_found"
    result_raw2="not_found"
  fi
  if [ -s "fastqc/$file_cln1" ]
  then
    fastqc_summry=$(unzip -l fastqc/$file_cln1 | grep fastqc_data.txt | awk '{ print $4 }' )
    result_cln1=$(unzip -p fastqc/$file_cln1 $fastqc_summry | grep "Total Sequences" | awk '{ print $3 }' )
  else
    file_cln1="not_found"
    result_cln1="not_found"
  fi
  if [ -s "fastqc/$file_cln2" ]
  then
    fastqc_summry=$(unzip -l fastqc/$file_cln2 | grep fastqc_data.txt | awk '{ print $4 }' )
    result_cln2=$(unzip -p fastqc/$file_cln2 $fastqc_summry | grep "Total Sequences" | awk '{ print $3 }' )
  else
    file_cln2="not_found"
    result_cln2="not_found"
  fi
  echo "$sample,$result_raw2,$result_cln2" >> $out/logs/raw_clean_scatter.csv
  FASTQC_RESULTS=(${FASTQC_RESULTS[@]} "$sample:$result_raw2:$result_cln2")
done
#echo ${FASTQC_RESULTS[@]}

# cg-pipeline results
CGPIPELINE_RESULTS=()
echo "sample x y" > $out/logs/raw_clean_coverage.txt
for sample in ${SAMPLES[@]}
do
  if [ -f "cg-pipeline/$sample.raw.out.txt" ]
  then
    cg_raw_coverage=$(grep "raw_shuffled" cg-pipeline/$sample.raw.out.txt | head -n 1 | awk '{print $9 }' )
    if [ -z "$cg_raw_coverage" ]; then cg_raw_coverage="no_result"; fi
  else
    cg_raw_coverage="no_result"
  fi
  if [ -f "cg-pipeline/$sample.clean.out.txt" ]
  then
    cg_cln_coverage=$(grep "clean_shuffled" cg-pipeline/$sample.clean.out.txt | head -n 1 | awk '{print $9 }' )
    if [ -z "$cg_cln_coverage" ]; then cg_cln_coverage="no_result"; fi
  else
    cg_cln_coverage="no_result"
  fi
  echo "$sample $cg_raw_coverage $cg_cln_coverage" >> $out/logs/raw_clean_coverage.txt
  CGPIPELINE_RESULTS=(${CGPIPELINE_RESULTS[@]} "$sample:$cg_raw_coverage:$cg_cln_coverage")
done
#echo ${CGPIPELINE_RESULTS[@]}

# mash results, seqsero results, and abricate serotype results
MASH_RESULTS=()
SEQSERO_RESULTS=()
SERO_ABRICATE_RESULTS=()
for sample in ${SAMPLES[@]}
do
  # mash_results
  if [ -n "$(find $out/mash -iname $sample*mashdist.txt )" ]
  then
    mash_results=($(cat $out/mash/$sample*mashdist.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2,3,4 -d "_" | cut -f 1 -d "." | tr "_" " " ))
    mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}""_""${mash_results[2]}""_""${mash_results[3]}")
    if [ -z "$mash_result" ]; then mash_result="no_result"; fi
    simple_mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}")
    if [ -z "$simple_mash_result" ]; then simple_mash_result="no_result"; fi
  else
    mash_result="no_result"
    simple_mash_result="no_result"
  fi
  MASH_RESULTS=(${MASH_RESULTS[@]} "$sample:$mash_result:$simple_mash_result")
  # seqsero results
  if [ -n "$(echo $mash_result | grep "Salmonella_enterica" )" ]
  then
    if [ -n "$(find $out/SeqSero -iname $sample*Seqsero_result.txt )" ]
    then
      seqsero_serotype=$(grep "Predicted serotype" $out/SeqSero/$sample.Seqsero_result.txt | cut -f 2 | tr ' ' '_' )
      if [ -z "$seqsero_serotype" ]; then seqsero_serotype="no_result"; fi
      seqsero_profile=$(grep "Predicted antigenic profile" $out/SeqSero/$sample.Seqsero_result.txt | cut -f 2 | tr ' ' '_' )
      if [ -z "$seqsero_profile" ]; then seqsero_profile="no_result"; fi
      simple_seqsero_result=$(echo $seqsero_serotype | sed 's/potentialmonophasicvariantof//g' | sed 's/potential_monophasic_variant_of_//g' | sed 's/O5-//g' | sed 's/See_comments_below/Enteritidis/g' )
      if [ -z "$simple_seqsero_result" ]; then simple_seqsero_result="no_result"; fi
      seqsero_serotype=$(echo $seqsero_serotype | sed 's/Enteritidis/Enteritidis(sdf+)/g' | sed 's/See_comments_below/Enteritidis(sdf-)/g' )
    else
      seqsero_serotype="no_result"
      seqsero_profile="no_result"
      simple_seqsero_result="no_result"
    fi
  else
    seqsero_serotype="not_salmonella"
    seqsero_profile="not_salmonella"
    simple_seqsero_result="not_salmonella"
  fi
  SEQSERO_RESULTS=(${SEQSERO_RESULTS[@]} "$sample;$seqsero_serotype;$seqsero_profile;$simple_seqsero_result")

  # abricate SerotypeFinder
  if [ -n "$(echo $mash_result | grep "Escherichia_coli" )" ]
  then
    if [ -f "$out/abricate_results/serotypefinder/serotypefinder.$sample.out.tab" ]
    then
      O_group_sero=($(grep $sample $out/abricate_results/serotypefinder/serotypefinder.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "O" | sed 's/\///g' ))
      H_group_sero=($(grep $sample $out/abricate_results/serotypefinder/serotypefinder.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "H" | sed 's/\///g' ))
      abricate_serotype_O=$(echo ${O_group_sero[@]} | tr ' ' '_' )
      if [ -z "$abricate_serotype_O" ]; then abricate_serotype_O="none"; fi
      abricate_serotype_H=$(echo ${H_group_sero[@]} | tr ' ' '_' )
      if [ -z "$abricate_serotype_H" ]; then abricate_serotype_H="none"; fi
    else
      abricate_serotype_H="not_ecoli"
      abricate_serotype_O="not_ecoli"
    fi
  else
    abricate_serotype_H="not_ecoli"
    abricate_serotype_O="not_ecoli"
  fi
  SERO_ABRICATE_RESULTS=(${SERO_ABRICATE_RESULTS[@]} "$sample:$abricate_serotype_O:$abricate_serotype_H")
done
#echo ${MASH_RESULTS[@]}
#echo ${SEQSERO_RESULTS[@]}
#echo ${SERO_ABRICATE_RESULTS[@]}

#abricate results : ncbi database
NCBI_ABRICATE_RESULTS=()
for sample in ${SAMPLES[@]}
do
  if [ -f "$out/abricate_results/ncbi/ncbi.$sample.out.tab" ]
  then
    abricate_results=($(grep $sample $out/abricate_results/ncbi/ncbi.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq | sed 's/[0]//g' ))
    abricate_result=$(echo ${abricate_results[@]} | tr ' ' '_' )
    if [ -z "$abricate_result" ]; then abricate_result="not_found"; fi
  else
    abricate_result="no_result"
  fi
  NCBI_ABRICATE_RESULTS=(${NCBI_ABRICATE_RESULTS[@]} "$sample:$abricate_result")
done
#echo ${NCBI_ABRICATE_RESULTS[@]}

STX_ABRICATE_RESULTS=()
for sample in ${SAMPLES[@]}
do
  if [ -f "$out/abricate_results/vfdb/vfdb.$sample.out.tab" ]
  then
    stxeae_results=($(grep $sample $out/abricate_results/vfdb/vfdb.$sample.out.tab | grep -e "stx" -e "eae" | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq | sed 's/[0]//g' ))
    stxeae_result=$(echo ${stxeae_results[@]} | tr ' ' '_' )
    if [ -z "$stxeae_result" ]; then stxeae_result="not_found"; fi
  else
    stxeae_result="not_found"
  fi
  STX_ABRICATE_RESULTS=(${STX_ABRICATE_RESULTS[@]} "$sample:$stxeae_result")
done
#echo ${STX_ABRICATE_RESULTS[@]}

BLOBTOOLS_RESULTS=()
for sample in ${SAMPLES[@]}
do
  if [ -f "$out/blobtools/$sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt" ]
  then
    blobtools_result=($(grep -v ^"#" $out/blobtools/$sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt | tr ' ' '_' | sort -k 3,3 -h | grep -v all | tail -n 1 | cut -f 1,13 ))
    if [ -z "$blobtools_result" ]; then blobtools_result="not_found"; fi
  else
    blobtools_result="not_found"
  fi
  BLOBTOOLS_RESULTS=(${BLOBTOOLS_RESULTS[@]} "$sample:${blobtools_result[0]}(${blobtools_result[1]})")
done
#echo ${BLOBTOOLS_RESULTS[@]}

# creating a heatmap for the files for easy visualization
echo "sample,fastqc,seqyclean,cg-pipeline,mash,shovill,prokka,quast,seqsero,abricate:serotypefinder,abricate:ncbi,blobtools" > $out/logs/File_heatmap.csv
for sample in ${SAMPLES[@]}
do
  sample=$sample # woot! A useless line
  fastqc_check=($(history -p ${FASTQC_RESULTS[@]} | sort | uniq | grep $sample | cut -f 2,3 -d ":" | tr ":" " " ))
  if [ "$fastq_check[0]" == "not_found" ]
  then
    fastqc_file="0"
  else
    fastqc_file="1"
  fi
  if [ "$fastq_check[1]" == "not_found" ]
  then
    seqyclean_file="0"
  else
    seqyclean_file="1"
  fi
  cg_check=$(history -p ${CGPIPELINE_RESULTS[@]} | sort | uniq | grep $sample | grep -v "not_found" | head -n 1 )
  if [ -z "$cg_check" ]
  then
    cg_file="0"
  else
    cg_file="1"
  fi
  mash_check=$(history -p ${MASH_RESULTS[@]} | sort | uniq | grep $sample | grep -v "not_found" | head -n 1 )
  if [ -z "$mash_check" ]
  then
    mash_file="0"
  else
    mash_file="1"
  fi
  if [ -s "shovill_result/$sample/contigs.fa" ]
  then
    shovill_file="1"
  else
    shovill_file="0"
  fi
  if [ -s "Prokka/$sample/$sample.gff" ]
  then
    prokka_file="1"
  else
    prokka_file="0"
  fi
  if [ -s "quast/$sample/report.tsv" ]
  then
    quast_file="1"
  else
    quast_file="0"
  fi
  if [ -s "SeqSero/$sample/Seqsero_result.txt" ]
  then
    seqsero_file="1"
  else
    seqsero_file="0"
  fi
  sero_check=$(history -p ${SERO_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | grep -v "not_found" | grep -v "not_ecoli" | head -n 1 )
  if [ -z "$sero_check" ]
  then
    serotypefinder_file="0"
  else
    serotypefinder_file="1"
  fi
  ncbi_check=$(history -p ${NCBI_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | grep -v "not_found" | head -n 1 )
  if [ -z "$ncbi_check" ]
  then
    ncbi_file="0"
  else
    ncbi_file="1"
  fi
  stx_check=$(history -p ${STX_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | grep -v "not_found" | head -n 1 )
  if [ -z "$stx_check" ]
  then
    stx_file="0"
  else
    stx_file="1"
  fi
  if [ -s "blobtools/$sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt" ]
  then
    blobtools_file="1"
  else
    blobtools_file="0"
  fi
  echo "$sample,$fastqc_file,$seqyclean_file,$cg_file,$mash_file,$shovill_file,$prokka_file,$quast_file,$seqsero_file,$serotypefinder_file,$ncbi_file,$blobtools_file" >> $out/logs/File_heatmap.csv
done
# Getting all the results in one file
echo -e "sample_id\tsample\tmash_result\tsimple_mash_result\tseqsero_serotype\tseqsero_profile\tsimple_seqsero_result\tcg_cln_coverage\tcg_raw_coverage\tfastqc_raw_reads_1\tfastqc_raw_reads_2\tfastqc_clean_reads_PE1\tfastqc_clean_reads_PE2\tabricate_ecoh_O\tabricate_ecoh_H\tabricate_serotype_O\tabricate_serotype_H\tstxeae_result\targannot\tresfinder\tcard\tplasmidfinder\tvfdb\tecoli_vf\tncbi\tblobtools_result" > run_results.txt
for sample in ${SAMPLES[@]}
do
  sample_id=$(echo $sample | rev | cut -f 4- -d "-" | rev )
  sample=$sample # woot! A useless line
  mash_result_split=($(history -p ${MASH_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  mash_result=${mash_result_split[1]}
  simple_mash_result=${mash_result_split[2]}
  seqsero_result_split=($(history -p ${SEQSERO_RESULTS[@]} | sort | uniq | grep $sample | tr ';' ' ' ))
  seqsero_serotype=${seqsero_result_split[1]}
  seqsero_profile=$(echo ${seqsero_result_split[2]} | tr ',' ';' )
  simple_seqsero_result=${seqsero_result_split[3]}
  cg_result_split=($(history -p ${CGPIPELINE_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  cg_cln_coverage=${cg_result_split[1]}
  cg_raw_coverage=${cg_result_split[2]}
  fastqc_result_split=($(history -p ${FASTQC_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  fastqc_raw_reads_1=${fastqc_result_split[1]}
  fastqc_raw_reads_2=${fastqc_result_split[1]}
  fastqc_clean_reads_PE1=${fastqc_result_split[2]}
  fastqc_clean_reads_PE2=${fastqc_result_split[2]}
  abricate_ecoh_O="X"
  abricate_ecoh_H="X"
  sero_result_split=($(history -p ${SERO_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  abricate_serotype_O=${sero_result_split[1]}
  abricate_serotype_H=${sero_result_split[2]}
  stxeae_result="X"
  argannot="X"
  resfinder="X"
  card="X"
  plasmidfinder="X"
  vfdb="X"
  stx_result_split=($(history -p ${STX_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  ecoli_vf=${stx_result_split[1]}
  ncbi_result_split=($(history -p ${NCBI_ABRICATE_RESULTS[@]} | sort | uniq | grep $sample | tr ':' ' ' ))
  ncbi=${ncbi_result_split[1]}
  blobtools=($(history -p ${BLOBTOOLS_RESULTS[@]} | sort | uniq | grep $sample | cut -f 2 -d ":" ))
  echo -e "$sample_id\t$sample\t$mash_result\t$simple_mash_result\t$seqsero_serotype\t$seqsero_profile\t$simple_seqsero_result\t$cg_cln_coverage\t$cg_raw_coverage\t$fastqc_raw_reads_1\t$fastqc_raw_reads_2\t$fastqc_clean_reads_PE1\t$fastqc_clean_reads_PE2\t$abricate_ecoh_O\t$abricate_ecoh_H\t$abricate_serotype_O\t$abricate_serotype_H\t$stxeae_result\t$argannot\t$resfinder\t$card\t$plasmidfinder\t$vfdb\t$ecoli_vf\t$ncbi\t$blobtools" >> run_results.txt
done

echo "Mash results count"
mash_column=$(head -n 1 $out/run_results.txt | tr "\t" "\n" | grep -n ^"simple_mash_result" | cut -f 1 -d ":" )
cut -f $mash_column $out/run_results.txt | awk '{if(NR>1)print}' | sed 's/.-//g' | sort | uniq -c | sort -k 1 -n | grep -v "no_result"

echo "Blobtools results count"
blobtools_column=$(head -n 1 $out/run_results.txt | tr "\t" "\n" | grep -n ^"blobtools_result" | cut -f 1 -d ":" )
cut -f $blobtools_column $out/run_results.txt | awk '{if(NR>1)print}' | sed 's/.-//g' | sed 's/(.*)//g' | sort | uniq -c | sort -k 1 -n | grep -v "no_result"

echo "Seqsero results count"
seqsero_column=$(head -n 1 $out/run_results.txt | tr "\t" "\n" | grep -n ^"simple_seqsero_result" | cut -f 1 -d ":" )
cut -f $seqsero_column $out/run_results.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" | grep -v "not_salmonella"

echo "Abricate serotype results count"
O_serotype_column=$(head -n 1 $out/run_results.txt | tr "\t" "\n" | grep -n ^"abricate_serotype_O" | cut -f 1 -d ":" )
H_serotype_column=$(head -n 1 $out/run_results.txt | tr "\t" "\n" | grep -n ^"abricate_serotype_H" | cut -f 1 -d ":" )
cut -f $O_serotype_column,$H_serotype_column $out/run_results.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"

date
echo "Finding each file for each sample complete!"
