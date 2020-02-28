#!/bin/bash

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

USAGE="
This is a script that takes the prokka, abricate, mash, and/or SeqSero output
from the Utah Reference Free Pipeline (aka URF of the Oakeson method) and
organizes it for roary and iqtree in order to create a phylogenetic tree of
selected organisms.

Version: 0.20200225

Usage: ./organize_URF_for_TREES.sh -i runs -o outbreak -f samples.txt

A full list of options is available with
./organize_URF_for_TREES.sh -h

"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

HELP="
The Utah Reference Free Pipeline (aka URF or the Oakeson method) generates files
with the following structure:

Version: 0.20200225

-RUN
  |-run_results.txt
  |-abricate_results
  | |-ncbi
  | | |-ncbi.{sample}.out.tab
  | |-serotypefinder
  |   |-serotypefinder.{sample}.out.tab
  |-ALL_gff
    |-{sample}.gff

This script organizes files into the following structure:

-DATE
  |-ecoli
  | |-{O group}
  |   |-{H group}
  |     |-gff
  |     |-abricate
  |-mash_results
  | |-{genus}
  |   |-{species}
  |     |-gff
  |     |-abricate
  |-Salmonella
    |-enterica
      |-{SeqSero Serotype}
        |-gff
        |-abricate

Example usage:

./organize_URF_for_TREES.sh -i /home/IDGenomics_NAS/WGS_Serotyping -o /home/IDGenomics_NAS/TREES -f outbreak_samples.txt

./organize_URF_for_TREES.sh -i /home/IDGenomics_NAS/WGS_Serotyping -o /home/IDGenomics_NAS/TREES -f outbreak_samples.txt -e failed_runs.txt -a 5

./organize_URF_for_TREES.sh -i /home/IDGenomics_NAS/WGS_Serotyping -o /home/IDGenomics_NAS/TREES -f outbreak_samples.txt -e failed_runs.txt -c /home/IDGenomics_NAS/WGS_Serotyping/tree_outgroups

REQUIRED OPTIONS:
  -i <URF path>             Path to parent directory of URF results
  -o <output path>          Path to directory where files will get organized to.
  -f <sample file>          Simple text file with one sample ID or regex match per line.

OPTIONAL OPTIONS:
  -e <file>                 Simple text file with one sample ID or regex match per line to be removed from analysis if found.
  -A                        Flag to include all available samples of genus/species
  -a <optional: number>     Specify maximum number of recent samples to include
  -c <path>                 Path to gff files to be included whenever directory structure matches organized samples

"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

# default variable values, mainly for Erin. Sorry!

input_files=()
exclude_file=""
run_date=$(date +%Y-%m-%d)
search_path="/home/IDGenomics_NAS/WGS_Serotyping"
out="/home/IDGenomics_NAS/UPHL_TREES"
min_cov=20 # On the todo list

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

while getopts 'Aa:c:e:f:hi:o:' OPTION
do
  case "$OPTION" in
    A)
    echo "Adding all genus/species matches to specified samples"
    add_flag=1
    ;;
    a)
    echo "Adding $OPTARG random genus/species matches to specified samples"
    random_samples=$OPTARG
    if [ -n "${random_samples//[0-9]}" ] && [ -n "$random_samples" ]
    then
      echo "$random_samples is not a valid number of samples"
      echo $USAGE
      exit
    fi
    ;;
    c)
    echo "Samples that must be included are found in $OPTARG"
    if [ -d "$OPTARG" ]
    then
      include_path=$OPTARG
    else
      echo "Directory containing samples could not be located"
      echo $USAGE
      exit
    fi
    ;;
    e)
    echo "The names of samples or regex to exclude are in $OPTARG"
    exclude_file=$OPTARG
    if [ ! -f "$exclude_file" ]
    then
      echo "$exclude_file could not be located"
      echo $USAGE
      exit
    fi
    ;;
    f)
    echo "Finding samples listed in $OPTARG"
    list_of_samples=$OPTARG
    if [ ! -f "$list_of_samples" ]
    then
      echo "$list_of_samples could not be located"
      echo "$USAGE"
      exit
    fi
    ;;
    h)
    echo "$HELP"
    exit 1
    ;;
    i)
    echo "Looking for files in $OPTARG"
    search_path=$OPTARG
    if [ ! -d "$search_path" ]
    then
      echo "Could not locate directory containing relevant files."
      echo $USAGE
      exit
    fi
    ;;
    o)
    echo "Outfiles and subfolders will be created in $OPTARG"
    out=$OPTARG
    ;;
    :)
    echo "Invalid option: $OPTARG requires an argument"
    echo "$USAGE"
    exit 1
    ;;
    \?)
    echo "$USAGE"
    exit 1
    ;;
  esac
done
shift "$(($OPTIND -1))"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

date
echo "Looking through $search_path for run_results.txt"

if [ -d "$out/$run_date" ]
then
  echo "$out/$run_date has already been created! Making new directory"
  i=1
  while [ -d "$out/$run_date" ]
  do
    run_date="$(date +%Y-%m-%d)""_$i"
    i=$(( $i + 1 ))
  done
  mkdir -p "$out/$run_date"
  if [ -d "$out/$run_date" ]
  then
    echo "The final directory is now $out/$run_date"
  else
    echo "Could not create a new directory"
    exit
  fi
fi

mkdir -p $out/$run_date/logs
echo -e "sample\treason" > $out/$run_date/logs/samples.rm

# creating shortcuts to gff, abricate, and seqsero results for samples
if [ -f "$exclude_file" ]
then
  grep -v "sample_id" $search_path/*/run_results.txt | grep -v "Undetermined" | grep -v -f $exclude_file | grep -f $list_of_samples | awk '{ if ( $9>20 ) print $0}' | sed 's/run_results.txt:/run_results.txt\t/g' > $out/$run_date/logs/sample_list.txt
else
  grep -v "sample_id" $search_path/*/run_results.txt | grep -v "Undetermined" |                            grep -f $list_of_samples | awk '{ if ( $9>20 ) print $0}' | sed 's/run_results.txt:/run_results.txt\t/g' > $out/$run_date/logs/sample_list.txt
fi

# As of February 25, 2020
# cg_cln_coverage = 9
mash_column=5
seqs_column=8
abrO_column=18
abrH_column=19

organize_samples () {
  while read line
  do
    sample=$(echo $line | cut -f 2 -d " " )
    run=$(echo $line | cut -f 1 -d " " | rev | cut -f 2 -d "/" | rev )
    if [ -n "$sample" ]
    then
      gff_file=$(ls -t $search_path/$run/ALL_gff/$sample*.gff | head -n 1) 2>> $out/$run_date/logs/errors.txt
      abr_file=$(ls -t $search_path/$run/abricate_results/ncbi/ncbi.$sample*.out.tab | head -n 1) 2>> $out/$run_date/logs/errors.txt
      if [ -n "$gff_file" ]
      then
        mash_sero=($(echo $line | sed 's/\t/ /g' | cut -f $mash_column -d " " | sed 's/_/\t/g'))
        seqs_sero=$(echo  $line | sed 's/\t/ /g' | cut -f $seqs_column -d " " | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/O5-//g')
        abrO_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrO_column -d " " | sed 's/O/ O/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
        abrH_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrH_column -d " " | sed 's/H/ H/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))

        if [ -z "${mash_sero[0]}" ] ; then mash_sero[0]="none" ; fi
        if [ -z "${mash_sero[1]}" ] ; then mash_sero[1]="none" ; fi
        if [ -z "$seqs_sero" ]      ; then seqs_sero="not_salmonella"    ; fi
        if [ -z "${abrO_sero[0]}" ] ; then abrO_sero=("notecoli")  ; fi
        if [ -z "${abrH_sero[0]}" ] ; then abrH_sero=("notecoli")  ; fi

        if [ "$seqs_sero" = "not_salmonella" ] && [ "$abrO_sero" = "notecoli" ]
        then
          mkdir -p $out/$run_date/mash_results/${mash_sero[0]}/{all,${mash_sero[1]}}/{gff,abricate}
          cp $gff_file $out/$run_date/mash_results/${mash_sero[0]}/all/gff/. 2>> $out/$run_date/logs/errors.txt
          cp $abr_file $out/$run_date/mash_results/${mash_sero[0]}/all/abricate/. 2>> $out/$run_date/logs/errors.txt
          cp $gff_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/gff/. 2>> $out/$run_date/logs/errors.txt
          cp $abr_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/abricate/. 2>> $out/$run_date/logs/errors.txt

        elif [ "$seqs_sero" = "not_salmonella" ] && [ "$abrO_sero" != "notecoli" ]
        then
          for O_group in ${abrO_sero[@]}
          do
            mkdir -p $out/$run_date/ecoli/$O_group/all/{gff,abricate}
            cp $gff_file $out/$run_date/ecoli/$O_group/all/gff/. 2>> $out/$run_date/logs/errors.txt
            cp $abr_file $out/$run_date/ecoli/$O_group/all/abricate/. 2>> $out/$run_date/logs/errors.txt
            for H_group in ${abrH_sero[@]}
            do
              mkdir -p $out/$run_date/ecoli/$O_group/$H_group/{gff,abricate}
              cp $gff_file $out/$run_date/ecoli/$O_group/$H_group/gff/. 2>> $out/$run_date/logs/errors.txt
              cp $abr_file $out/$run_date/ecoli/$O_group/$H_group/abricate/. 2>> $out/$run_date/logs/errors.txt
            done
          done
        elif [ "$seqs_sero" != "not_salmonella" ] && [ "$abrO_sero" = "notecoli" ]
        then
          mkdir -p $out/$run_date/Salmonella/enterica/$seqs_sero/{gff,abricate}
          cp $gff_file $out/$run_date/Salmonella/enterica/$seqs_sero/gff/. 2>> $out/$run_date/logs/errors.txt
          cp $abr_file $out/$run_date/Salmonella/enterica/$seqs_sero/abricate/. 2>> $out/$run_date/logs/errors.txt
        else
          echo -e "$sample\tundetermined_serotype_in_$run" >> $out/$run_date/logs/samples.rm
        fi
      else
        echo -e "$sample\tgff_not_found_in_$run" >> $out/$run_date/logs/samples.rm
      fi
    fi
  done < $1
}

organize_samples $out/$run_date/logs/sample_list.txt

date
echo "Removing gff files that are too small"
ls $out/$run_date/*/*/*/gff -d | parallel wc -l {}/*gff | awk '{ if ( $1 < 10000) print $2 "\tsize_too_small"}' >> $out/$run_date/logs/samples.rm
ls $out/$run_date/*/*/*/gff -d | parallel wc -l {}/*gff | grep -v "total" | awk '{ if ( $1 < 10000) print $2 }' | parallel rm {}

if [ -n "$add_flag" ] || [ -n "$random_samples" ]
then
  date
  echo "Adding additional samples"
  directories=($(ls -d $out/$run_date/*/*/*/gff | sed 's/gff//g' | grep -v "ecoli" | grep -v "all" ))
  for directory in ${directories[@]}
  do
    if [ -z "$exclude_file" ] ; then exclude_file=$list_of_samples ; fi
    species=($(echo $directory | rev | cut -f 1-4 -d "/" | rev | sed 's/\// /g' | sed 's/mash_results//g' ))
    grep_command="$(echo ${species[@]:1} | sed 's/ / | grep /g')"
    if [ -z "$grep_command" ] ; then grep_command="." ; fi
    if [ -n "$random_samples" ]
    then
      eval "grep ${species[0]} $search_path/*/run_results.txt | grep $grep_command " | sed 's/\/run_results.txt/-run_results.txt/g' | sort -r -t "-" -k3 | grep -v -f $list_of_samples | grep -v -f $exclude_file | awk '{ if ( $9>20 ) print $0}' | head -n $random_samples | sed 's/-run_results.txt:/\/run_results.txt\t/g' >> $out/$run_date/logs/sample_addition_list.txt
    else
      eval "grep ${species[0]} $search_path/*/run_results.txt | grep $grep_command " |                                                                     grep -v -f $list_of_samples | grep -v -f $exclude_file | awk '{ if ( $9>20 ) print $0}' |                           sed 's/run_results.txt:/run_results.txt\t/g'    >> $out/$run_date/logs/sample_addition_list.txt
    fi
  done

  ecoli_directories=($(ls -d $out/$run_date/ecoli/*/*/gff | sed 's/gff//g'))
  for directory in ${ecoli_directories[@]}
  do
    if [ -z "$exclude_file" ] ; then exclude_file=$list_of_samples ; fi
    species=($(echo $directory | rev | cut -f 1-4 -d "/" | rev | sed 's/\// /g' | sed 's/mash_results//g' | sed 's/all//g'))
    grep_command="$(echo ${species[@]:1} | sed 's/ / | grep /g')"
    if [ -z "$grep_command" ] ; then grep_command="." ; fi
    if [ -n "$random_samples" ]
    then
      grep Escherichia_coli $search_path/*/run_results.txt | eval "grep $grep_command " | sed 's/\/run_results.txt/-run_results.txt/g' | sort -r -t "-" -k3 | grep -v -f $list_of_samples | grep -v -f $exclude_file | awk '{ if ( $9>20 ) print $0}' | head -n $random_samples | sed 's/-run_results.txt:/\/run_results.txt\t/g' >> $out/$run_date/logs/sample_addition_list.txt
    else
      grep Escherichia_coli $search_path/*/run_results.txt | eval "grep $grep_command " |                                                                     grep -v -f $list_of_samples | grep -v -f $exclude_file | awk '{ if ( $9>20 ) print $0}' |                           sed 's/run_results.txt:/run_results.txt\t/g'    >> $out/$run_date/logs/sample_addition_list.txt
    fi
  done
  cat $out/$run_date/logs/sample_addition_list.txt | sort | uniq > $out/$run_date/logs/sample_addition_list.txt_temp
  mv $out/$run_date/logs/sample_addition_list.txt_temp $out/$run_date/logs/sample_addition_list.txt
  organize_samples $out/$run_date/logs/sample_addition_list.txt
fi

date
echo "Removing redundant directories"
directories=($(ls -d $out/$run_date/*/*/*/gff | grep -v "all/gff" | grep -v "Salmonella/enterica" | sort | uniq ))
for specific_directory in ${directories[@]}
do
  cut_directory=$(echo $specific_directory | rev | cut -f 3- -d "/" | rev)
  diff_check=$(diff -q $specific_directory $cut_directory/all/gff | grep "Only" | head -n 1 )
  if [ -z "$diff_check" ]
  then
    ls $cut_directory/all/gff/*gff | awk '{ print $0 "\tdirectory_redundant" }' >> $out/$run_date/logs/samples.rm
    rm -R $cut_directory/all
  fi
done

date
echo "Creating list_of_samples.txt for each directory and removing empty directories"
directories=($(ls -d $out/$run_date/*/*/*/gff ))
for directory in ${directories[@]}
do
  if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
  then
    rm -R $directory
  else
    newdirectory=$(echo $directory | sed 's/gff//g')
    ls $directory/*gff > $newdirectory/list_of_samples.txt
  fi
done

date
echo "Removing directories with too few samples (Roary & IQTREE require 4)"
small_directories=($(wc -l $out/$run_date/*/*/*/list_of_samples.txt | awk '{if ( $1<4 ) print $2}' | sed 's/list_of_samples.txt//g'))
for directory in ${small_directories[@]}
do
  ls $directory/gff/*gff | awk '{ print $0 "\ttoo_few_samples"}' >> $out/$run_date/logs/samples.rm
  rm -R $directory
done

date
echo "Checking samples"
missing_samples=($(cat $list_of_samples | parallel ls $out/$run_date/*/*/*/gff/*{}* 2>&1 | grep "cannot access" | cut -f 4 -d " " | rev | cut -f 1 -d "/" | rev | sed 's/\*//g' | sed 's/://g'))
echo "Missing samples:"
history -p ${missing_samples[@]}
#echo -e "\nSample\tReason"
#for sample in ${missing_samples[@]}
#do
#  grep $sample $out/$run_date/logs/samples.rm)
#done

find $out/$run_date -type d -empty -delete

if [ -d "$include_path" ]
then
  date
  echo "Including samples found in $include_path"
  ecli_directories=($(ls -d $out/$run_date/ecoli/*/*/ | rev | cut -f 1-3 -d "/" | rev ))
  salm_directories=($(ls -d $out/$run_date/Salmonella/enterica/*/ | rev | cut -f 2 -d "/" | rev ))
  mash_directories=($(ls -d $out/$run_date/mash_results/*/*/ | rev | cut -f 1-3 -d "/" | rev ))
  for ecli_directory in ${ecli_directories[@]}
  do
    O_H=($(echo $ecli_directory | sed 's/\// /g'))
    if [ -d "$include_path/ecoli/${O_H[0]}/${O_H[1]}" ] && [ -d "$out/$run_date/ecoli/${O_H[0]}/${O_H[1]}" ]
    then
      cp $include_path/ecoli/${O_H[0]}/${O_H[1]}/abricate/*out.tab $out/$run_date/ecoli/${O_H[0]}/${O_H[1]}/abricate/.
      cp $include_path/ecoli/${O_H[0]}/${O_H[1]}/gff/*gff          $out/$run_date/ecoli/${O_H[0]}/${O_H[1]}/gff/.
    fi
    if [ -d "$include_path/ecoli/${O_H[0]}" ] && [ -d "$out/$run_date/ecoli/${O_H[0]}/all" ]
    then
      cp $include_path/ecoli/${O_H[0]}/*/abricate/*out.tab $out/$run_date/ecoli/${O_H[0]}/all/abricate/.
      cp $include_path/ecoli/${O_H[0]}/*/gff/*gff          $out/$run_date/ecoli/${O_H[0]}/all/gff/.
    fi
  done
  for salm_serotype in ${salm_directories[@]}
  do
    if [ -d "$include_path/Salmonella/enterica/$salm_serotype" ] && [ -d "$out/$run_date/Salmonella/enterica/$salm_serotype" ]
    then
      cp $include_path/Salmonella/enterica/$salm_serotype/abricate/*out.tab $out/$run_date/Salmonella/enterica/$salm_serotype/abricate/.
      cp $include_path/Salmonella/enterica/$salm_serotype/gff/*gff          $out/$run_date/Salmonella/enterica/$salm_serotype/gff/.
    fi
  done
  for mash_directory in ${mash_directories[@]}
  do
    genus_species=($(echo $mash_directory | sed 's/\// /g'))
    if [ -d "$include_path/mash_results/${genus_species[0]}/${genus_species[1]}" ] && [ -d "$out/$run_date/mash_results/${genus_species[0]}/${genus_species[1]}" ]
    then
      cp $include_path/mash_results/${genus_species[0]}/${genus_species[1]}/abricate/*out.tab $out/$run_date/mash_results/${genus_species[0]}/${genus_species[1]}/abricate/.
      cp $include_path/mash_results/${genus_species[0]}/${genus_species[1]}/gff/*gff          $out/$run_date/mash_results/${genus_species[0]}/${genus_species[1]}/gff/.
      if [ -d "$out/$run_date/mash_results/${genus_species[0]}/all" ]
      then
        cp $include_path/mash_results/${genus_species[0]}/${genus_species[1]}/abricate/*out.tab $out/$run_date/mash_results/${genus_species[0]}/all/abricate/.
        cp $include_path/mash_results/${genus_species[0]}/${genus_species[1]}/gff/*gff          $out/$run_date/mash_results/${genus_species[0]}/all/gff/.
      fi
    fi
  done
fi

date
echo "log files can be found at $out/$run_date/logs"
echo "Organizing gff files for UPHL_TREES is finished!"

exit 0
