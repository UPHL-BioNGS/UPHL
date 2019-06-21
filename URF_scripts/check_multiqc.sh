#!/bin/bash
out=$1

echo "sample,x,y" > $out/logs/raw_clean_scatter.csv
echo "sample x y" > $out/logs/raw_clean_coverage.txt

FILE_TEMPLATE=(
"raw_fastq_files:Sequencing_reads/Raw:.f"
"clean_fastq_files:Sequencing_reads/QCed:clean_*.fastq"
"fastqc_files:fastqc:_fastqc."
"mash_file:mash:.fastq.msh.distance.sorted.txt"
"shovill_file:ALL_assembled:_contigs.fa"
"plasmid_shovill_file:ALL_assembled_plasmids:_plasmidcontigs.fa"
"prokka_file:ALL_gff:.gff"
"plasmid_prokka_file:ALL_gff_plasmids:_plasmid.gff"
"seqsero_file:SeqSero:.Seqsero_result.txt"
"cg-pipeline_files:cg-pipeline:.out.txt"
"quast_file:quast:report.txt"
"plasmid_quast:quast_plasmids:report.txt"
"abricate_all:abricate_results:.out.tab"
"plasmid_abricate_all:abricate_results_plasmids:plasmids.out.tab"
)

wc -l abricate_results*/*/*out.tab | awk '{ if ( $1 == 1 ) print $2 }' | parallel " rm $out/{} ; touch $out/{} "
grep -P "Predicted antigenic profile:\t-:-:-" SeqSero/*Seqsero_result.txt | cut -f 1 -d ":" | parallel " rm $out/{} ; touch $out/{} "
analysis_types=($(history -p ${FILE_TEMPLATE[@]} | cut -f 1 -d ":" ))
echo "sample\t${analysis_types[@]}" | tr ' ' '\t' | parallel "echo -e {} > $out/run_file_summary.txt "
file_summary_header=($(head -n 1 $out/run_file_summary.txt ))
echo ${file_summary_header[@]} | tr ' ' ',' > $out/logs/File_heatmap.csv

RESULTS=(
"pnusa"
"sample"
"mash_result"
"simple_mash_result"
"seqsero_serotype"
"seqsero_profile"
"simple_seqsero_result"
"cg_cln_coverage"
"cg_raw_coverage"
"fastqc_raw_reads_1"
"fastqc_raw_reads_2"
"fastqc_clean_reads_PE1"
"fastqc_clean_reads_PE2"
"abricate_ecoh_O"
"abricate_ecoh_H"
"abricate_serotype_O"
"abricate_serotype_H"
"stxeae_result"
)
#ABRICATE_DATABASES=("argannot" "resfinder" "card" "plasmidfinder" "vfdb" "ecoli_vf" "ncbi" "ecoh" "serotypefinder")
for database in argannot resfinder card plasmidfinder vfdb ecoli_vf ncbi
do
  RESULTS=($(echo -e "${RESULTS[@]}\t$database" ))
done
echo ${RESULTS[@]} | tr ' ' '\t' | parallel " echo -e {} > $out/run_results_summary.txt "
SAMPLES=($(ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq ))
for sample in ${SAMPLES[@]}
do
  pnusa=$(echo $sample | sed 's/-UT.*//g' )
  FILES=("sample:$sample")
  for template in ${FILE_TEMPLATE[@]}
  do
    template_information=($(echo $template | awk -F ":" '{print $1 " " $2 " " $3}' ))
    files=($(find $out/${template_information[1]} -wholename *$pnusa*${template_information[2]}* -size +0 -type f ))
    if [ -n "${files[0]}" ]
    then
      file=$(echo ${files[@]} | tr ' ' ',' | tr '\t' ',' | tr '\n' ',' )
    else
      file="not_found,"
    fi
    FILES=("$(echo -e "${FILES[@]}\t${template_information[0]}:$file" )")
  done

  #########################################heatmap_of_files and summary_of_files

  file_summary=""
  file_heatmap=""
  for column in ${file_summary_header[@]}
  do
    file=$(history -p ${FILES[@]} | sort | uniq | awk -F ":" -v analysis=$column '{ if ($1==analysis) print $2}' | tr '\n' ' ' )
    if [ -z "$file_summary" ]
    then
      file_summary="$file"
    else
      file_summary="$(echo -e "$file_summary\t$file" )"
    fi

    if [ -z "$file_heatmap" ]
    then
      file_heatmap="$file"
    else
      if [ -z "$file" ]
      then
        file_heatmap="$file_heatmap,0"
      elif [ -n "$(echo $file | grep not_ecoli )" ]
      then
        file_heatmap="$file_heatmap,0.5"
      elif [ -n "$(echo $file | grep not_salmonella )" ]
      then
        file_heatmap="$file_heatmap,0.5"
      elif [ -n "$(echo $file | grep not_found )" ]
      then
        file_heatmap="$file_heatmap,0"
      elif [ -n "$(echo $file | grep no_result )" ]
      then
        file_heatmap="$file_heatmap,0.5"
      else
        file_heatmap="$file_heatmap,1"
      fi
    fi
  done

  echo "$file_heatmap" | sed 's/,,/,/g' >> $out/logs/File_heatmap.csv
  echo "$file_summary" | sed 's/,,/,/g' >> $out/run_file_summary.txt

  #########################################results_from_files

  # mash_results
  if [ -n "$(find $out/mash -iname $sample*sorted.txt )" ]
  then
#    mash_results=($(cat $out/mash/$sample*sorted.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2,3,4 -d "_" | cut -f 1 -d "." | awk -F "_" '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' ))
    mash_results=($(cat $out/mash/$sample*sorted.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2,3,4 -d "_" | cut -f 1 -d "." | tr "_" " " ))
    mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}""_""${mash_results[2]}""_""${mash_results[3]}")
    if [ -z "$mash_result" ]; then mash_result="no_result"; fi
    simple_mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}")
    if [ -z "$simple_mash_result" ]; then simple_mash_result="no_result"; fi
  else
    mash_result="no_result"
    simple_mash_result="no_result"
  fi

  # seqsero_results
  if [ -n "$(echo $mash_result | grep "Salmonella_enterica" )" ]
  then
    if [ -n "$(find $out/SeqSero -iname $sample*Seqsero_result.txt )" ]
    then
      seqsero_serotype=$(grep "Predicted serotype" $out/SeqSero/$sample.Seqsero_result.txt | cut -f 2 | tr ' ' '_' )
      if [ -z "$seqsero_serotype" ]; then seqsero_serotype="no_result"; fi
      seqsero_profile=$(grep "Predicted antigenic profile" $out/SeqSero/$sample.Seqsero_result.txt | cut -f 2 | tr ' ' '_' )
      if [ -z "$seqsero_profile" ]; then seqsero_profile="no_result"; fi
      simple_seqsero_result=$(echo $seqsero_serotype | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/potential_monophasic_variant_of_//g' | sed 's/O5-//g' )
      if [ -z "$simple_seqsero_result" ]; then simple_seqsero_result="no_result"; fi
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

  # cg-pipeline_results
  if [ -f "$out/cg-pipeline/$sample.raw.out.txt" ]
  then
    cg_raw_coverage=$(grep "raw_shuffled" $out/cg-pipeline/$sample.raw.out.txt | grep "$pnusa" | head -n 1 | awk '{print $9 }' )
    if [ -z "$cg_raw_coverage" ]; then cg_raw_coverage="no_result"; fi
  else
    cg_raw_coverage="no_result"
  fi
  if [ -f "$out/cg-pipeline/$sample.clean.out.txt" ]
  then
    cg_cln_coverage=$(grep "clean_shuffled" $out/cg-pipeline/$sample.clean.out.txt | grep "$pnusa" | head -n 1 | awk '{print $9 }' )
    if [ -z "$cg_cln_coverage" ]; then cg_cln_coverage="no_result"; fi
  else
    cg_cln_coverage="no_result"
  fi

  # fastqc_results
  FASTQC_RESULTS=()
  if [ -n "$(find $out/fastqc -iname $sample*zip )" ]
  then
    fastqc_files=(ls $out/fastqc/$sample*zip)
    for fastqc_file in ${fastqc_files[@]}
    do
      if [ -f "$fastqc_file" ]
      then
        fastqc_summry=$(unzip -l $fastqc_file | grep fastqc_data.txt | awk '{ print $4 }' )
        fastqc_result=$(unzip -p $fastqc_file $fastqc_summry | grep "Total Sequences" | awk '{ print $3 }' )
        if [ -z "$fastqc_result" ]; then fastqc_result="not_found"; fi
        FASTQC_RESULTS=($(echo -e "${FASTQC_RESULTS[@]}\t$fastqc_file:$fastqc_result" ))
      else
        FASTQC_RESULTS=($(echo -e "${FASTQC_RESULTS[@]}\t$fastqc_file:no_result" ))
      fi
    done
    fastqc_raw_reads_1=$(history -p ${FASTQC_RESULTS[@]} | grep $sample | grep -v "shuffled" | grep -v "clean" | grep -v "ls:no_result" | cut -f 2 -d ":" | head -n 1 )
    fastqc_raw_reads_2=$(history -p ${FASTQC_RESULTS[@]} | grep $sample | grep -v "shuffled" | grep -v "clean" | grep -v "ls:no_result" | cut -f 2 -d ":" | tail -n 1 )
    fastqc_clean_reads_PE1=$(history -p ${FASTQC_RESULTS[@]} | grep $sample | grep -v "shuffled" | grep "clean_PE1.fastq" | grep -v "ls:no_result" | cut -f 2 -d ":" )
    fastqc_clean_reads_PE2=$(history -p ${FASTQC_RESULTS[@]} | grep $sample | grep -v "shuffled" | grep "clean_PE2.fastq" | grep -v "ls:no_result" | cut -f 2 -d ":" )
  else
    fastqc_raw_reads_1="no_result"
    fastqc_raw_reads_2="no_result"
    fastqc_clean_reads_PE1="no_result"
    fastqc_clean_reads_PE2="no_result"
  fi

  if [ -n "$(echo $mash_result | grep "Escherichia_coli" )" ]
  then
    O_group_ecoh=($(grep $pnusa $out/abricate_results/ecoh/ecoh.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "O" | sed 's/\///g' ))
    H_group_ecoh=($(grep $pnusa $out/abricate_results/ecoh/ecoh.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "H" | sed 's/\///g' ))
    abricate_ecoh_O=$(echo ${O_group_ecoh[@]} | tr ' ' '_' )
    if [ -z "$abricate_ecoh_O" ]; then abricate_ecoh_O="none"; fi
    abricate_ecoh_H=$(echo ${H_group_ecoh[@]} | tr ' ' '_' )
    if [ -z "$abricate_ecoh_H" ]; then abricate_ecoh_H="none"; fi

    O_group_sero=($(grep $pnusa $out/abricate_results/serotypefinder/serotypefinder.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "O" | sed 's/\///g' ))
    H_group_sero=($(grep $pnusa $out/abricate_results/serotypefinder/serotypefinder.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | awk -F "_" '{print $NF}' | awk -F "-" '{print $NF}' | sort | uniq | grep "H" | sed 's/\///g' ))
    abricate_serotype_O=$(echo ${O_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_O" ]; then abricate_serotype_O="none"; fi
    abricate_serotype_H=$(echo ${H_group_sero[@]} | tr ' ' '_' )
    if [ -z "$abricate_serotype_H" ]; then abricate_serotype_H="none"; fi
  else
    abricate_serotype_H="not_ecoli"
    abricate_serotype_O="not_ecoli"
    abricate_ecoh_H="not_ecoli"
    abricate_ecoh_O="not_ecoli"
  fi

  if [ -f "$out/abricate_results/vfdb/vfdb.$sample.out.tab" ]
  then
    stxeae_results=($(grep $pnusa $out/abricate_results/vfdb/vfdb.$sample.out.tab | grep -e "stx" -e "eae" | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq | sed 's/[0]//g' ))
    stxeae_result=$(echo ${stxeae_results[@]} | tr ' ' '_' )
    if [ -z "$stxeae_result" ]; then stxeae_result="not_found"; fi
  fi

  results="$pnusa\t$sample\t$mash_result\t$simple_mash_result\t$seqsero_serotype\t$seqsero_profile\t$simple_seqsero_result\t$cg_cln_coverage\t$cg_raw_coverage\t$fastqc_raw_reads_1\t$fastqc_raw_reads_2\t$fastqc_clean_reads_PE1\t$fastqc_clean_reads_PE2\t$abricate_ecoh_O\t$abricate_ecoh_H\t$abricate_serotype_O\t$abricate_serotype_H\t$stxeae_result"
  # abricate_results
  for database in argannot resfinder card plasmidfinder vfdb ecoli_vf ncbi
  do
    if [ -f "$out/abricate_results/$database/$database.$sample.out.tab" ]
    then
      abricate_results=($(grep $sample $out/abricate_results/$database/$database.$sample.out.tab | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq | sed 's/[0]//g' ))
      abricate_result=$(echo ${abricate_results[@]} | tr ' ' '_' )
      if [ -z "$abricate_result" ]; then abricate_result="not_found"; fi
      results="$results\t$abricate_result"
    else
      results="$results\tno_result"
    fi
  done
  echo -e $results >> $out/run_results_summary.txt
  echo "$sample,$fastqc_raw_reads_2,$fastqc_clean_reads_PE2" | sed 's/no_result/0/g' >> $out/logs/raw_clean_scatter.csv
  echo "$sample $cg_raw_coverage $cg_cln_coverage" | sed 's/no_result/0/g' >> $out/logs/raw_clean_coverage.txt
done
echo "Results for $out"

echo "Mash results count"
mash_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"simple_mash_result" | cut -f 1 -d ":" )
cut -f $mash_column $out/run_results_summary.txt | awk '{if(NR>1)print}' | sed 's/.-//g' | sort | uniq -c | sort -k 1 -n | grep -v "no_result"

echo "Seqsero results count"
seqsero_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"simple_seqsero_result" | cut -f 1 -d ":" )
cut -f $seqsero_column $out/run_results_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" | grep -v "not_salmonella"

echo "Abricate serotype results count"
O_serotype_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"abricate_serotype_O" | cut -f 1 -d ":" )
H_serotype_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"abricate_serotype_H" | cut -f 1 -d ":" )
cut -f $O_serotype_column,$H_serotype_column $out/run_results_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"

echo "Abricate ecoh results count"
O_ecoh_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"abricate_ecoh_O" | cut -f 1 -d ":" )
H_ecoh_column=$(head -n 1 $out/run_results_summary.txt | tr "\t" "\n" | grep -n ^"abricate_ecoh_H" | cut -f 1 -d ":" )
cut -f $O_ecoh_column,$H_ecoh_column $out/run_results_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"

date
echo "Finding each file for each sample complete!"
