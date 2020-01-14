#!/bin/bash

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

USAGE="
This is a script that takes the output from prokka, abricate, mash, and/or
SeqSero through roary, and iqtree in order to create a phylogenetic tree of
organisms with files created less than 120 days ago.

This current script is dated September 19, 2019 by Erin Young.

Usage: ./organize_URF_for_OUTBREAK_120.sh -i /files/location -o /output/directory/

A full list of options is available with
./organize_URF_for_OUTBREAK_120.sh -h

"
#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
HELP="
OUTBREAK_120 looks through a directory (or directories) for contig files.

REQUIRED OPTIONS:

-i <path to look through for */run_results_summary.txt>

OUTBREAK_120 will search this directory for files matching the following patterns:
Summary files:    <path>/*/run_results_summary.txt
Prokka files:     <path>/*/ALL_gff/NAME.gff
Abricate results: <path>/*/abricate_results/ncbi/ncbi.*.out.txt

-o <path for outfolder>

OPTIONAL OPTIONS:

-e <file>                         List of sample names to be removed from analysis if found.
-p <species>                      Restrict to specific species.
-s <yyyy-mm-dd>                   Specify end date.
-f <path to file>                 Specify samples via a file.
-d <number of days>               Specify number of days prior to end date (default is 120)
-m <number of samples>            Adjust maximum number of samples (default is 100)
-c <path to control info>         Use file similar to results file for controls
-a                                Evaluates all subfolders, not just those with samples less than 10 days old.
-r                                Expands maximum number of samples for directories with too few samples.

"

control_file_help="
The control file needs the following columns:
* sample                with full path of gff file of control
* simple_mash_result    with Species_genus of control as determined by mash
* simple_seqsero_result with Salmonella_enteria_predicted_serotype
* abricate_serotype_O   O group as determined via abricate with the serotypefinder database
* abricate_serotype_H   H group as determined via abricate with the serotypefinder database
"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

# default variable values, mainly for Erin. Sorry!

input_files=()
exclude_file=""
d=$(date +%Y-%m-%d)
run_date=$(date +%Y-%m-%d)
age=120
pd=$(date -d "$d - $age days" +%Y-%m-%d)
rd=$(date -d "$d - 10 days" +%Y-%m-%d)
control_file=/home/Bioinformatics/Data/NCBI_refs/OUTBREAK_120_controls/controls_file.txt
control_path=/home/Bioinformatics/Data/NCBI_refs/OUTBREAK_120_controls/ASSEMBLED_GENOMES # On the todo list
max_samples=100
search_path="/home/IDGenomics_NAS/WGS_Serotyping"
out="/home/IDGenomics_NAS/WGS_Serotyping/OUTBREAK_120"
species="."
min_cov=20 # On the todo list

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

while getopts 'ac:C:d:e:f:hi:m:o:p:rs:t:' OPTION
do
  case "$OPTION" in
    a)
    echo "Creating directories for all samples, not just those with samples less than 10 days old"
    recent_flag=1
    ;;
    C)
    echo "Looking for controls in $OPTARG"
    control_path=$OPTARG
    echo "This is still in development. Sorry if you were hoping this would be functional!"
    echo "email eriny@utah.gov to get working on this"
    exit
    ;;
    c)
    echo "Control information is found in $OPTARG"
    if [ ! -f "$OPTARG" ]
    then
      echo "$OPTARG not found!"
      echo "$USAGE"
      echo "$control_file_help"
      exit
    else
      control_file=$OPTARG
    fi
    ;;
    s)
    echo "Ending date is now $OPTARG instead of $d"
    d=$OPTARG
    end_flag=1
    ;;
    d)
    echo "Looking for files that have been modified $OPTARG days prior"
    age=$OPTARG
    pd=$(date -d "$d - $age days" +%Y-%m-%d)
    if [ -n "${age//[0-9]}" ] && [ -n "$age" ]
    then
      echo "$age is not a valid number of days"
      echo $USAGE
      exit
    fi
    ;;
    e)
    echo "The names of samples to exclude are in $OPTARG"
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
    m)
    max_samples=$OPTARG
    if [ -n "${max_samples//[0-9]}" ] && [ -n "$max_samples" ]
    then
      echo "$max_samples is not a valid number of samples"
      echo $USAGE
      exit
    fi
    echo "The maximum number of samples has changed from 100 to $max_samples"
    ;;
    o)
    echo "Outfiles and subfolders will be created in $OPTARG"
    out=$OPTARG
    ;;
    p)
    species=$OPTARG
    echo "Will look through results for samples that are $species"
    ;;
    r)
    expand_flag=1
    echo "For directories with less than 4 samples, will look for files older than $age days"
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
echo "Looking through $search_path for gff files between $d and $pd"

if [ -d "$out/$run_date" ]
then
  echo "$out/$run_date has already been created! Making new directory"
  i=1
  while [ -d "$out/$run_date" ]
  do
    run_date="$(date +%Y-%m-%d)""_$i"
    i=$(( $i + 1 ))
  done
  mkdir "$out/$run_date"
  if [ -d "$out/$run_date" ]
  then
    echo "The final directory is now $out/$run_date"
  else
    echo "Could not create a new directory"
    exit
  fi
fi

mkdir -p $out/$run_date/logs
mkdir -p $out/$run_date/serotyping_results/abricate
echo -e "sample\treason" > $out/$run_date/logs/samples.rm

# finding the run_result summaries
if [ -n "$end_flag" ]
then
  list_of_run_result_files=($(find $search_path -maxdepth 2 -newermt $pd ! -newermt $d -name run_results*))
else
  list_of_run_result_files=($(find $search_path -maxdepth 2 -newermt $pd -name run_results*))
fi

# getting the exclude file in place
if [ -z "$exclude_file" ]
then
  touch $out/$run_date/logs/exclude_file.txt
else
  cp $exclude_file $out/$run_date/logs/exclude_file.txt
fi

if [ -z "$list_of_samples" ]
then
  echo "." > $out/$run_date/logs/sample_list.txt
else
  cp $list_of_samples $out/$run_date/logs/sample_list.txt
fi

# creating shortcuts to gff, abricate, and seqsero results for samples
for result_file in ${list_of_run_result_files[@]}
do
  time_stamp=$(date)
  echo "$time_stamp: Looking through samples in $result_file"
  mash_column=$(head -n 1 $result_file | tr "\t" "\n" | grep -n "simple_mash_result"    | cut -f 1 -d ":" )
  seqs_column=$(head -n 1 $result_file | tr "\t" "\n" | grep -n "simple_seqsero_result" | cut -f 1 -d ":" )
  abrO_column=$(head -n 1 $result_file | tr "\t" "\n" | grep -n "abricate_serotype_O"   | cut -f 1 -d ":" )
  abrH_column=$(head -n 1 $result_file | tr "\t" "\n" | grep -n "abricate_serotype_H"   | cut -f 1 -d ":" )
  cvrg_column=$(head -n 1 $result_file | tr "\t" "\n" | grep -n "cg_cln_coverage"       | cut -f 1 -d ":" )

  run=$(echo $result_file | rev | cut -f 2 -d "/" | rev )

  ln -s $search_path/$run/abricate_results/ncbi/ncbi*out.tab $out/$run_date/serotyping_results/abricate/. 2>> $out/$run_date/logs/shortcut_overlap.txt

  while read line
  do
    sample=$(echo $line | awk -v col="$cvrg_column" '{ if ( $col>20 ) print $1}' )
    if [ -n "$sample" ]
    then
      if [ -n "$end_flag" ]
      then
        gff_file=$(find $search_path/$run/ALL_gff -maxdepth 1 -newermt $pd ! -newermt $d -iname $sample*.gff | head -n 1 )
      else
        gff_file=$(find $search_path/$run/ALL_gff -maxdepth 1 -newermt $pd -iname $sample*.gff | head -n 1 )
      fi

      if [ -n "$gff_file" ]
      then
        mash_sero=($(echo $line | sed 's/\t/ /g' | cut -f $mash_column -d " " | sed 's/_/\t/g'))
        seqs_sero=$(echo  $line | sed 's/\t/ /g' | cut -f $seqs_column -d " " | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/O5-//g')
        abrO_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrO_column -d " " | sed 's/O/ O/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
        abrH_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrH_column -d " " | sed 's/H/ H/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))

        if [ -z "${mash_sero[0]}" ] ; then mash_sero[0]="none" ; fi
        if [ -z "${mash_sero[1]}" ] ; then mash_sero[1]="none" ; fi
        if [ -z "$seqs_sero" ]      ; then seqs_sero="none"    ; fi
        if [ -z "${abrO_sero[0]}" ] ; then abrO_sero=("none")  ; fi
        if [ -z "${abrH_sero[0]}" ] ; then abrH_sero=("none")  ; fi

        mkdir -p $out/$run_date/mash_results/${mash_sero[0]}/all
        mkdir -p $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}
        ln -s $gff_file $out/$run_date/mash_results/${mash_sero[0]}/all/. 2>> $out/$run_date/logs/shortcut_overlap.txt
        ln -s $gff_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/. 2>> $out/$run_date/logs/shortcut_overlap.txt

        mkdir -p $out/$run_date/Salmonella/enterica/all
        mkdir -p $out/$run_date/Salmonella/enterica/$seqs_sero
        ln -s $gff_file $out/$run_date/Salmonella/enterica/all/. 2>> $out/$run_date/logs/shortcut_overlap.txt
        ln -s $gff_file $out/$run_date/Salmonella/enterica/$seqs_sero/. 2>> $out/$run_date/logs/shortcut_overlap.txt

        for O_group in ${abrO_sero[@]}
        do
          mkdir -p $out/$run_date/ecoli/$O_group/all
          ln -s $gff_file $out/$run_date/ecoli/$O_group/all/. 2>> $out/$run_date/logs/shortcut_overlap.txt
          for H_group in ${abrH_sero[@]}
          do
            mkdir -p $out/$run_date/ecoli/$O_group/$H_group
            ln -s $gff_file $out/$run_date/ecoli/$O_group/$H_group/. 2>> $out/$run_date/logs/shortcut_overlap.txt
          done
        done
      else
        echo -e "$sample-$run\tno_gff_in_timespan" >> $out/$run_date/logs/samples.rm
      fi
    fi
  done < <(grep -v pnusa "$result_file" | grep -v "sample_id" | grep -v "Undetermined" | grep "$species" | grep -v -f $out/$run_date/logs/exclude_file.txt | grep -f $out/$run_date/logs/sample_list.txt )
done

date
echo "Now removing redundant serotypes"
# first: remove the notecoli, notsalmonella, notmash stuff
if [ -d "$out/$run_date/Salmonella/enterica/no_result" ]
then
  ls $out/$run_date/Salmonella/enterica/no_result/*gff      2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_seqsero" }'                 >> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/Salmonella/enterica/no_result
fi
if [ -d "$out/$run_date/Salmonella/enterica/not_salmonella" ]
then
  ls $out/$run_date/Salmonella/enterica/not_salmonella/*gff 2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_salmonella" }'             >> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/Salmonella/enterica/not_salmonella
fi
if [ -d "$out/$run_date/mash_results/none" ]
then
  ls $out/$run_date/mash_results/none/*/*gff                2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_mash" }'                    >> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/mash_results/none
fi
if [ -d "$(find $out/$run_date/mash_results/ -name none | head -n 1 )" ]
then
  ls $out/$run_date/mash_results/*/none/*gff                2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_mash" }'                    >> $out/$run_date/logs/samples.rm
  ls $out/$run_date/mash_results/*/none/ -d | parallel rm -R {}
fi
if [ -d "$out/$run_date/ecoli/none" ]
then
  ls $out/$run_date/ecoli/none/*/*gff                       2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_abricate_serotypefinder" }' >> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/ecoli/none
fi
if [ -d "$(find $out/$run_date/ecoli/ -name none | head -n 1 )" ]
then
  ls $out/$run_date/ecoli/*/none/*gff                       2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_abricate_serotypefinder" }' >> $out/$run_date/logs/samples.rm
  ls $out/$run_date/ecoli/*/none/ -d | parallel rm -R {}
fi
if [ -d "$out/$run_date/ecoli/notecoli" ]
then
  ls $out/$run_date/ecoli/notecoli/*/*gff                   2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_ecoli" }'                  >> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/ecoli/notecoli
fi
if [ -d "$(find $out/$run_date/ecoli/ -name notecoli | head -n 1 )" ]
then
  ls $out/$run_date/ecoli/*/notecoli/*gff                   2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_ecoli" }'                  >> $out/$run_date/logs/samples.rm
  ls $out/$run_date/ecoli/*/notecoli/ -d | parallel rm -R {}
fi
if [ -d "$out/$run_date/ecoli/no" ]
then
  ls $out/$run_date/ecoli/no/*/*gff 2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tunknown_serotype" }'>> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/ecoli/no
fi
if [ -d "$out/$run_date/ecoli/result" ]
then
  ls $out/$run_date/ecoli/result/*/*gff 2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tunknown_serotype" }'>> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/ecoli/result
fi
if [ -d "$out/$run_date/Salmonella/enterica/all" ]
then
  ls $out/$run_date/Salmonella/enterica/all/*gff 2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tlisted_as_salmonella" }'>> $out/$run_date/logs/samples.rm
  rm -R $out/$run_date/Salmonella/enterica/all
fi

# to avoid duplicates
if [ -d "$out/$run_date/Salmonella/enterica" ] & [ -d "$out/$run_date/mash_results/Salmonella/enterica" ] ; then rm -R $out/$run_date/mash_results/Salmonella/enterica ; fi
if [ -d "$out/$run_date/ecoli" ]               & [ -d "$out/$run_date/mash_results/Escherichia/coli" ]    ; then rm -R $out/$run_date/mash_results/Escherichia/coli    ; fi

# to remove files that are "too small"
ls $out/$run_date/*/*/*/ -d | parallel wc -l {}*gff | awk '{ if ( $1 < 10000) print $2 "\tsize_too_small"}' >> $out/$run_date/logs/samples.rm
ls $out/$run_date/*/*/*/ -d | parallel wc -l {}*gff | grep -v "total" | awk '{ if ( $1 < 10000) print $2 }' | parallel rm {}

date
echo "Creating list_of_samples.txt for each directory and removing empty directories"
directories=($(ls -d $out/$run_date/*/*/*/ | sed s/.$// ))
for directory in ${directories[@]}
do
  if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
  then
    rm -R $directory
  else
    ls $directory/*gff > $directory/list_of_samples.txt
    if [ -z "$recent_flag" ]
    then
      while read line
      do
        file_age=$(stat -L $line | grep "Modify" | awk '{print $2}' )
        if [[ "$rd" < "$file_age" ]] ; then echo $line >> $directory/recent_samples.txt ; fi
      done < $directory/list_of_samples.txt
      if [ ! -f $directory/recent_samples.txt ] ; then ls $directory/*gff | awk '{ print $0 "\tno_recent" }' >> $out/$run_date/logs/samples.rm ; rm -R $directory ; fi
    fi
  fi
done

date
echo "Removing redundant directories"
directories=($(ls -d $out/$run_date/*/*/*/ | sed s/.$// | grep -v "/all" | grep -v "Salmonella/enterica"))
for specific_directory in ${directories[@]}
do
  cut_directory=$(echo $specific_directory | rev | cut -f 2- -d "/" | rev)
  diff_check=$(diff -q $specific_directory $cut_directory/all | grep "Only" | head -n 1 )
  if [ -z "$diff_check" ]
  then
    ls $cut_directory/all/*gff | awk '{ print $0 "\tdirectory_redundant" }' >> $out/$run_date/logs/samples.rm
    rm -R $cut_directory/all
#    rm -R $specific_directory
  fi
done

date
echo "Removing directories with more than $max_samples samples or too few samples"
directories=($(ls -d $out/$run_date/*/*/*/ | sed s/.$// ))
for directory in ${directories[@]}
do
  number_of_files=$(ls $directory/*gff | wc -l )
  if (( $number_of_files > $max_samples ))
  then
    ls $directory/*gff | awk '{ print $0 "\ttoo_many_samples" }' >> $out/$run_date/logs/samples.rm
    rm -R $directory
  elif (( $number_of_files < 4 ))
  then
    if [ -z "$expand_flag" ]
    then
      echo "removing $directory due to only containing $number_of_files files"
      ls $directory/*gff | awk '{ print $0 "\ttoo_few_samples" }' >> $out/$run_date/logs/samples.rm
      rm -R $directory
    else
      date
      echo "Looking for additional samples for $directory"
      grep_options=($(echo $directory | rev | cut -f 1-3 -d "/" | rev | tr '/' ' ' | sed 's/notecoli//g' | sed 's/ecoli/Escherichia coli/g' | sed 's/all/./g' | sed 's/mash_results//g' | sed 's/^ //g' | sed 's/ / | grep /g' ))
      additives=($(ls -t $search_path/*/run_results* | parallel --jobs 1 grep ${grep_options[@]} {} | grep "UT-" | awk -v col="$cvrg_column" '{ if ( $col>20 ) print $2}' | head -n 20 | sort | uniq ))
      for addition in ${additives[@]}
      do
        gff_file=$(ls -t $search_path/*/ALL_gff/$addition*.gff | head -n 1 )
        if [ -n "$gff_file" ]
        then
          wc -l $gff_file | awk '{ if ( $1 > 10000) print $2 }' | parallel ln -s {} $directory/. 2>> $out/$run_date/logs/shortcut_overlap.txt
        fi
      done
      number_of_files=$(ls $directory/*gff | wc -l )
      if (( $number_of_files < 4 ))
      then
        echo "removing $directory due to only containing $number_of_files sample(s) (need 4)"
        ls $directory/*gff | awk '{ print $0 "\ttoo_few_samples" }' >> $out/$run_date/logs/samples.rm
        rm -R $directory
      fi
    fi
  fi
done

find $out/$run_date -type d -empty -delete

date
echo "Adding controls to samples"
if [ -f "$control_file" ]
then
  date
  echo "Adding controls for outgroups listed in $control_file"
  mash_column=$(head -n 1 $control_file | tr "\t" "\n" | grep -n "simple_mash_result"    | cut -f 1 -d ":" )
  seqs_column=$(head -n 1 $control_file | tr "\t" "\n" | grep -n "simple_seqsero_result" | cut -f 1 -d ":" )
  abrO_column=$(head -n 1 $control_file | tr "\t" "\n" | grep -n "abricate_serotype_O"   | cut -f 1 -d ":" )
  abrH_column=$(head -n 1 $control_file | tr "\t" "\n" | grep -n "abricate_serotype_H"   | cut -f 1 -d ":" )
  gfff_column=$(head -n 1 $control_file | tr "\t" "\n" | grep -n "sample"                | cut -f 1 -d ":" )

  while read line
  do
    mash_sero=($(echo $line | sed 's/\t/ /g' | cut -f $mash_column -d " " | sed 's/_/\t/g'))
    seqs_sero=$(echo  $line | sed 's/\t/ /g' | cut -f $seqs_column -d " " | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/O5-//g')
    abrO_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrO_column -d " " | sed 's/O/ O/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
    abrH_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrH_column -d " " | sed 's/H/ H/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
    gfff_file=$(echo  $line | sed 's/\t/ /g' | cut -f $gfff_column -d " " )

    if [ -f "$gfff_file" ]
    then
      # mash_results
      if [ -d "$out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}" ]
      then
        echo "Adding $gfff_file"
        ln -s $gfff_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/. 2>> $out/$run_date/logs/shortcut_overlap.txt
        if [ -d "$out/$run_date/mash_results/${mash_sero[0]}/all" ]
        then
          ln -s $gfff_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/. 2>> $out/$run_date/logs/shortcut_overlap.txt
        fi
      fi
      # seqsero_results
      if [ -d "$out/$run_date/Salmonella/enterica/$seqs_sero" ]
      then
        echo "Adding $gfff_file"
        ln -s $gfff_file $out/$run_date/Salmonella/enterica/$seqs_sero/. 2>> $out/$run_date/logs/shortcut_overlap.txt
      fi
      #abricate_results
      for O_group in ${abrO_sero[@]}
      do
        if [ -d "$out/$run_date/ecoli/$O_group" ]
        then
                  echo "Adding $gfff_file"
          if [ -d "$out/$run_date/ecoli/$O_group/all" ] ; then ln -s $gfff_file $out/$run_date/ecoli/$O_group/all/. 2>> $out/$run_date/logs/shortcut_overlap.txt ; fi
          for H_group in ${abrH_sero[@]}
          do
            if [ -d "$out/$run_date/ecoli/$O_group/$H_group" ]
            then
              ln -s $gfff_file $out/$run_date/ecoli/$O_group/$H_group/. 2>> $out/$run_date/logs/shortcut_overlap.txt
            fi
          done
        fi
      done
    else
      echo "gff file for $gfff_file not found. Did you include the full path?"
    fi
  done < <(grep -v "#" $control_file)
elif [ -n $(find $control_path -iname "*gff" ) ]
then
  date
  echo "Adding controls for outgropus found in $control_path"
  echo "This is still in development. Sorry if you were hoping this would be functional!"
  echo "email eriny@utah.gov to get working on this"
else
  date
  echo "Not adding control files."
fi

find $out/$run_date -type d -empty -delete

date
echo "log files can be found at $out/$run_date/logs"
echo "Organizing gff files for OUTBREAK_120 is finished!"

exit
