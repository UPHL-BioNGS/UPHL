#!/bin/bash

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

USAGE="
This is a script that takes the output from prokka, abricate, mash, and/or
SeqSero through roary, and iqtree in order to create a phylogenetic tree of
organisms with files created less than 120 days ago.

This current script is dated September 24, 2018 by Erin Young.

Usage: ./outbreak_120_organize.sh -i /files/location -o /output/directory/

A full list of options is available with
./outbreak_120_organize_2.sh -h
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
-p <species>                      Restrict specific species.
-s <yyyy-mm-dd>                   Specify end date.
-d <number of days>               Specify number of days prior to end date (default is 120)
-m <number of samples>            Adjust maximum number of samples (default is 100)
-a                                Evaluates all subfolders, not just those with samples less than 10 days old.
-P <PATH>                         Add path to PATH. Can be used to specify multiple paths.

"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

# default variable values

input_files=()
exclude_file=""
d=$(date +%Y-%m-%d)
run_date=$(date +%Y-%m-%d)
age=120
pd=$(date -d "$d - $age days" +%Y-%m-%d)
rd=$(date -d "$d - 10 days" +%Y-%m-%d)
control_path=/home/Bioinformatics/Data/NCBI_refs/OUTBREAK_120_controls/ASSEMBLED_GENOMES
max_samples=100
search_path="/home/IDGenomics_NAS/WGS_Serotyping"
out="/home/IDGenomics_NAS/WGS_Serotyping/OUTBREAK_120"
species="."

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

while getopts 'ac:d:e:f:hi:m:o:p:s:t:' OPTION
do
  case "$OPTION" in
    a)
    echo "Creating directories for all samples, not just those with samples less than 10 days old"
    recent_flag=1
    ;;
    c)
    echo "The controls are located in $OPTARG"
    control_path=$OPTARG
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

run_control ()
{
  out=$1
  d=$2
  control_path=$3

  date
  echo "Adding NCBI genome controls for outgroups"
  control_list=(
  Acinetobacter_baumannii_FA_control.gff:mash_results:Acinetobacter:baumannii
  Acinetobacter_nosocomialis_control.gff:mash_results:Acinetobacter:nosocomialis
  Acinetobacter_pittii_FA_control.gff:mash_results:Acinetobacter:pittii
  Bacteroides_fragilis_control.gff:mash_results:Bacteroides:fragilis
  Campylobacter_coli_FA_control.gff:mash_results:Campylobacter:coli
  Campylobacter_fetus_FA_control.gff:mash_results:Campylobacter:fetus
  Campylobacter_hyointestinalis_FA_control.gff:mash_results:Campylobacter:hyointestinalis
  Campylobacter_jejuni_FA_control.gff:mash_results:Campylobacter:jejuni
  Campylobacter_lari_FA_control.gff:mash_results:Campylobacter:lari
  Campylobacter_peloridis_control.gff:mash_results:Campylobacter:peloridis
  Campylobacter_upsaliensis_FA_control.gff:mash_results:Campylobacter:upsaliensis
  Clostridium_perfringens_FA_control.gff:mash_results:Clostridium:perfringens
  Cronobacter_sakazakii_control.gff:mash_results:Cronobacter:sakazakii
  Cronobacter_turicensis_control.gff:mash_results:Cronobacter:turicensis
  Enterobacter_aerogenes_control.gff:mash_results:Enterobacter:aerogenes
  Escherichia_albertii_control.gff:mash_results:Escherichia:albertii
  Escherichia_coli_K12_control.gff:ecoli:none:none
  Escherichia_coli_O103H11_control.gff:ecoli:O103:H11
  Escherichia_coli_O103H25_control.gff:ecoli:O103:H25
  Escherichia_coli_O103H2_control.gff:ecoli:O103:H2
  Escherichia_coli_O111H8_FA_control.gff:ecoli:O111:H8
  Escherichia_coli_O113H21_control.gff:ecoli:O113:H21
  Escherichia_coli_O118H16_control.gff:ecoli:O118:H16
  Escherichia_coli_O119H4_control.gff:ecoli:O119:H4
  Escherichia_coli_O121H19_FA_control.gff:ecoli:O121:H19
  Escherichia_coli_O123H11_control.gff:ecoli:O123:H11
  Escherichia_coli_O128H27_control.gff:ecoli:O128:H27
  Escherichia_coli_O128H2_control.gff:ecoli:O128:H2
  Escherichia_coli_O145H25_control.gff:ecoli:O145:H25
  Escherichia_coli_O145H28_control.gff:ecoli:O145:H28
  Escherichia_coli_O157H7_FA_control.gff:ecoli:O157:H7
  Escherichia_coli_O15H11_control.gff:ecoli:O15:H11
  Escherichia_coli_O165H25_control.gff:ecoli:O165:H25
  Escherichia_coli_O174H21_control.gff:ecoli:O174:H21
  Escherichia_coli_O177_control.gff:ecoli:O177:none
  Escherichia_coli_O25bH4_control.gff:ecoli:O25:H4
  Escherichia_coli_O26H11_FA_control.gff:ecoli:O26:H11
  Escherichia_coli_O45H2_control.gff:ecoli:O45:H2
  Escherichia_coli_O5H4_control.gff:ecoli:O5:H4
  Escherichia_coli_O69H11_control.gff:ecoli:O69:H11
  Escherichia_coli_O81_control.gff:ecoli:O81:none
  Helicobacter_pullorum_control.gff:mash_results:Helicobacter:pullorum
  Klebsiella_aerogenes_control.gff:mash_results:Klebsiella:aerogenes
  Klebsiella_oxytoca_control.gff:mash_results:Klebsiella:oxytoca
  Klebsiella_pneumoniae_control.gff:mash_results:Klebsiella:pneumoniae
  Legionella_pneumophila_control.gff:mash_results:Legionella:pneumophila
  Listeria_monocytogenes_FA_control.gff:mash_results:Listeria:monocytogenes
  Pantoea_ananatis_control.gff:mash_results:Pantoea:ananatis
  Pseudomonas_aeruginosa_control.gff:mash_results:Pseudomonas:aeruginosa
  Salmonella_bongori_control.gff:mash_results:Salmonella:bongori
  Salmonella_enterica_Aberdeen_control.gff:salmonella:enterica:Aberdeen
  Salmonella_enterica_Abony_control.gff:salmonella:enterica:Abony
  Salmonella_enterica_Adelaide_control.gff:salmonella:enterica:Adelaide
  Salmonella_enterica_Agbeni_contigs_control.gff:salmonella:enterica:BronorAgbeni
  Salmonella_enterica_Agona_control.gff:salmonella:enterica:Agona
  Salmonella_enterica_Alachua_contigs_control.gff:salmonella:enterica:Alachua
  Salmonella_enterica_Albany_control.gff:salmonella:enterica:Albany
  Salmonella_enterica_Altona_contigs_control.gff:salmonella:enterica:Altona
  Salmonella_enterica_Anatum_control.gff:salmonella:enterica:Anatum
  Salmonella_enterica_Apapa_control.gff:salmonella:enterica:Apapa
  Salmonella_enterica_Baildon_contigs_control.gff:salmonella:enterica:Baildon
  Salmonella_enterica_Bareilly_control.gff:salmonella:enterica:Bareilly
  Salmonella_enterica_Berta_contigs_control.gff:salmonella:enterica:Berta
  Salmonella_enterica_Blockley_contigs_control.gff:salmonella:enterica:HaardtorBlockley
  Salmonella_enterica_Bovismorbificans_control.gff:salmonella:enterica:HindmarshorBovismorbificans
  Salmonella_enterica_Braenderup_control.gff:salmonella:enterica:Braenderup
  Salmonella_enterica_Brandenburg_contigs_control.gff:salmonella:enterica:Brandenburg
  Salmonella_enterica_Bredeney_control.gff:salmonella:enterica:Bredeney
  Salmonella_enterica_Cerro_control.gff:salmonella:enterica:Cerro
  Salmonella_enterica_Chester_control.gff:salmonella:enterica:Chester
  Salmonella_enterica_Choleraesuis_control.gff:salmonella:enterica:Choleraesuis
  Salmonella_enterica_Concord_control.gff:salmonella:enterica:Concord
  Salmonella_enterica_Corvallis_control.gff:salmonella:enterica:Corvallis
  Salmonella_enterica_Cubana_control.gff:salmonella:enterica:Cubana
  Salmonella_enterica_Daytona_contigs_control.gff:salmonella:enterica:Daytona
  Salmonella_enterica_Derby_control.gff:salmonella:enterica:Derby
  Salmonella_enterica_Dublin_control.gff:salmonella:enterica:Dublin
  Salmonella_enterica_Ealing_contigs_control.gff:salmonella:enterica:Ealing
  Salmonella_enterica_Eastbourne_contigs_control.gff:salmonella:enterica:Eastbourne
  Salmonella_enterica_Edinburgh_contigs_control.gff:salmonella:enterica:Edinburgh
  Salmonella_enterica_Enteritidis_control.gff:salmonella:enterica:Enteritidis
  Salmonella_enterica_Gaminara_control.gff:salmonella:enterica:Gaminara
  Salmonella_enterica_Give_control.gff:salmonella:enterica:Give
  Salmonella_enterica_Glostrup_contigs_control.gff:salmonella:enterica:GlostruporChomedey
  Salmonella_enterica_Grumpensis_contigs_control.gff:salmonella:enterica:Grumpensis
  Salmonella_enterica_Hadar_contigs_control.gff:salmonella:enterica:HadarorIstanbul
  Salmonella_enterica_Hartford_contigs_control.gff:salmonella:enterica:Hartford
  Salmonella_enterica_Havana_contigs_control.gff:salmonella:enterica:Havana
  Salmonella_enterica_Heidelberg_control.gff:salmonella:enterica:Heidelberg
  Salmonella_enterica_Hvittingfoss_control.gff:salmonella:enterica:II16benxorHvittingfoss
  Salmonella_enterica_Ibadan_contigs_control.gff:salmonella:enterica:IbadanorMississippiorII11323b15
  Salmonella_enterica_Indiana_control.gff:salmonella:enterica:Indiana
  Salmonella_enterica_Infantis_control.gff:salmonella:enterica:Infantis
  Salmonella_enterica_Istanbul_contigs_control.gff:salmonella:enterica:HadarorIstanbul
  Salmonella_enterica_Javiana_control.gff:salmonella:enterica:II912lz2815orJaviana
  Salmonella_enterica_Johannesburg_control.gff:salmonella:enterica:Johannesburg
  Salmonella_enterica_Kentucky_control.gff:salmonella:enterica:Kentucky
  Salmonella_enterica_Kiambu_contigs_control.gff:salmonella:enterica:KiambuorII141227z15
  Salmonella_enterica_Kintambo_contigs_control.gff:salmonella:enterica:KintamboorWashington
  Salmonella_enterica_Kottbus_contigs_control.gff:salmonella:enterica:Kottbus
  Salmonella_enterica_Litchfield_contigs_control.gff:salmonella:enterica:PakistanorLitchfield
  Salmonella_enterica_Livingstone_contigs_control.gff:salmonella:enterica:Livingstone
  Salmonella_enterica_Lomalinda_contigs_control.gff:salmonella:enterica:LomalindaorII1912aenx
  Salmonella_enterica_London_contigs_control.gff:salmonella:enterica:London
  Salmonella_enterica_Manhattan_control.gff:salmonella:enterica:YovokomeorManhattan
  Salmonella_enterica_Mbandaka_control.gff:salmonella:enterica:Mbandaka
  Salmonella_enterica_Meleagridis_contigs_control.gff:salmonella:enterica:Meleagridis
  Salmonella_enterica_Miami_contigs_control.gff:salmonella:enterica:Miami
  Salmonella_enterica_Michigan_contigs_control.gff:salmonella:enterica:Michigan
  Salmonella_enterica_Minnesota_control.gff:salmonella:enterica:Minnesota
  Salmonella_enterica_Mississippi_contigs_control.gff:salmonella:enterica:IbadanorMississippiorII11323b15
  Salmonella_enterica_Monschaui_contigs_control.gff:salmonella:enterica:Monschaui
  Salmonella_enterica_Montevideo_control.gff:salmonella:enterica:Montevideo
  Salmonella_enterica_Muenchen_contigs_control.gff:salmonella:enterica:VirginiaorMuenchen
  Salmonella_enterica_Muenster_control.gff:salmonella:enterica:Muenster
  Salmonella_enterica_Napoli_contigs_control.gff:salmonella:enterica:Napoli
  Salmonella_enterica_Newport_control.gff:salmonella:enterica:Newport
  Salmonella_enterica_Norwich_contigs_control.gff:salmonella:enterica:Norwich
  Salmonella_enterica_Ohio_FA_control.gff:salmonella:enterica:Ohio
  Salmonella_enterica_Onderstepoort_contigs_control.gff:salmonella:enterica:BahrenfeldorOnderstepoort
  Salmonella_enterica_Oranienburg_control.gff:salmonella:enterica:OranienburgorII67mt-
  Salmonella_enterica_Oslo_contigs_control.gff:salmonella:enterica:Oslo
	Salmonella_enterica_Panama_control.gff:salmonella:enterica:PanamaorHouston
	Salmonella_enterica_ParatyphiA_control.gff:salmonella:enterica:ParatyphiA
	Salmonella_enterica_ParatyphiB_control.gff:salmonella:enterica:ParatyphiB
	Salmonella_enterica_ParatyphiC_control.gff:salmonella:enterica:ParatyphiC
	Salmonella_enterica_Pomona_control.gff:salmonella:enterica:Pomona
	Salmonella_enterica_Poona_control.gff:salmonella:enterica:FarmsenorPoona
	Salmonella_enterica_Potsdam_contigs_control.gff:salmonella:enterica:Potsdam
	Salmonella_enterica_Reading_contigs_control.gff:salmonella:enterica:Reading
	Salmonella_enterica_Rissen_contigs_control.gff:salmonella:enterica:Rissen
	Salmonella_enterica_Rubislaw_control.gff:salmonella:enterica:Rubislaw
	Salmonella_enterica_Saintpaul_control.gff:salmonella:enterica:Saintpaul
	Salmonella_enterica_Sandiego_contigs_control.gff:salmonella:enterica:Sandiego
	Salmonella_enterica_Schwarzengrund_control.gff:salmonella:enterica:Schwarzengrund
	Salmonella_enterica_Senftenberg_control.gff:salmonella:enterica:SenftenbergorDessau
	Salmonella_enterica_Singapore_contigs_control.gff:salmonella:enterica:Singapore
	Salmonella_enterica_Stanley_control.gff:salmonella:enterica:Stanley
	Salmonella_enterica_Stanleyville_control.gff:salmonella:enterica:Stanleyville
	Salmonella_enterica_Sundsvall_contigs_control.gff:salmonella:enterica:Sundsvall
	Salmonella_enterica_Tallahassee_contigs_control.gff:salmonella:enterica:Tallahassee
	Salmonella_enterica_Telelkebir_contigs_control.gff:salmonella:enterica:Telelkebir
	Salmonella_enterica_Tennessee_control.gff:salmonella:enterica:II67z29z42orTennessee
	Salmonella_enterica_Thompson_control.gff:salmonella:enterica:Thompson
	Salmonella_enterica_Typhi_FA_control.gff:salmonella:enterica:Typhi
	Salmonella_enterica_Typhimurium_FA_control.gff:salmonella:enterica:Typhimurium
	Salmonella_enterica_Uganda_contigs_control.gff:salmonella:enterica:Uganda
	Salmonella_enterica_Urbana_contigs_control.gff:salmonella:enterica:Urbana
	Salmonella_enterica_Virchow_contigs_control.gff:salmonella:enterica:Virchow
	Salmonella_enterica_Wandsworth_control.gff:salmonella:enterica:Wandsworth
	Salmonella_enterica_Weltevreden_control.gff:salmonella:enterica:Weltevreden
	Salmonella_enterica_Worthington_contigs_control.gff:salmonella:enterica:Worthington
	Serratia_marcescens_control.gff:mash_results:Serratia:marcescens
	Shigella_flexneri_control.gff:mash_results:Shigella:flexneri
	Shigella_sonnei_FA_control.gff:mash_results:Shigella:sonnei
	Variovorax_paradoxus_control.gff:mash_results:Variovorax:paradoxus
	Vibrio_navarrensis_control.gff:mash_results:Vibrio:navarrensis
	Vibrio_parahaemolyticus_control.gff:mash_results:Vibrio:parahaemolyticus
  )

  for control in ${control_list[@]}
  do
    control_parts=($(echo "$control_path/$control" | awk -F ":" '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' ))
    if [ -d "$out/$d/${control_parts[1]}/all/all" ] ; then ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/all/all/. 2>> $out/$d/logs/shortcut_overlap.txt ; fi

    if [ -d "$out/$d/${control_parts[1]}/${control_parts[2]}/all" ]
    then
      echo "Adding ${control_parts[0]}"
      ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/${control_parts[2]}/all/. 2>> $out/$d/logs/shortcut_overlap.txt

      if [ -d "$out/$d/${control_parts[1]}/${control_parts[2]}/${control_parts[3]}" ]
      then
        ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/${control_parts[2]}/${control_parts[3]}/. 2>> $out/$d/logs/shortcut_overlap.txt
      fi
    fi
  done

  date
  echo "Adding controls has completed."
}

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
date
echo "Looking through $search_path for gff files between $d and $pd"

mkdir -p $out/$run_date/logs
mkdir -p $out/$run_date/serotyping_results/abricate
echo -e "sample\treason" > $out/$run_date/logs/samples.rm
if [ -n "$end_flag" ]
then
  list_of_run_result_files=($(find $search_path -maxdepth 2 -newermt $pd ! -newermt $d -name run_results_summary.txt))
else
  list_of_run_result_files=($(find $search_path -maxdepth 2 -newermt $pd -name run_results_summary.txt))
fi

if [ -z "$exclude_file" ]
then
  touch $out/$run_date/logs/exclude_file.txt
else
  cp $exclude_file $out/$run_date/logs/exclude_file.txt
fi

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
        seqs_sero=$(echo  $line | sed 's/\t/ /g' | cut -f $seqs_column -d " " )
        abrO_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrO_column -d " " | sed 's/O/ O/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
        abrH_sero=($(echo $line | sed 's/\t/ /g' | cut -f $abrH_column -d " " | sed 's/H/ H/g' | sed 's/not_ecoli/notecoli/g' | sed 's/_/ /g'))
        linkfile=$(echo "$sample-$run.gff")
        if [ -z "${mash_sero[0]}" ] ; then mash_sero[0]="none" ; fi
        if [ -z "${mash_sero[1]}" ] ; then mash_sero[1]="none" ; fi
        if [ -z "$seqs_sero" ]      ; then seqs_sero="none"    ; fi
        if [ -z "${abrO_sero[0]}" ] ; then abrO_sero=("none")  ; fi
        if [ -z "${abrH_sero[0]}" ] ; then abrH_sero=("none")  ; fi

        mkdir -p $out/$run_date/mash_results/${mash_sero[0]}/all
        mkdir -p $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}
        ln -s $gff_file $out/$run_date/mash_results/${mash_sero[0]}/all/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt
        ln -s $gff_file $out/$run_date/mash_results/${mash_sero[0]}/${mash_sero[1]}/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt

        mkdir -p $out/$run_date/salmonella/enterica/all
        mkdir -p $out/$run_date/salmonella/enterica/$seqs_sero
        ln -s $gff_file $out/$run_date/salmonella/enterica/all/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt
        ln -s $gff_file $out/$run_date/salmonella/enterica/$seqs_sero/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt

        for O_group in ${abrO_sero[@]}
        do
          mkdir -p $out/$run_date/ecoli/$O_group/all
          ln -s $gff_file $out/$run_date/ecoli/$O_group/all/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt
          for H_group in ${abrH_sero[@]}
          do
            mkdir -p $out/$run_date/ecoli/$O_group/$H_group
            ln -s $gff_file $out/$run_date/ecoli/$O_group/$H_group/$linkfile 2>> $out/$run_date/logs/shortcut_overlap.txt
          done
        done
      else
        echo -e "$sample-$run\tno_gff_found" >> $out/$run_date/logs/samples.rm
      fi
    fi
  done < <(grep -v pnusa "$result_file" | grep -v "Undetermined" | grep "$species" | grep -v -f $out/$run_date/logs/exclude_file.txt )
done

date
echo "Now removing redundant serotypes"
# first: remove the notecoli, notsalmonella, notmash stuff
ls $out/$run_date/salmonella/enterica/no_result/*gff      2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_seqsero" }'                 >> $out/$run_date/logs/samples.rm
ls $out/$run_date/salmonella/enterica/not_salmonella/*gff 2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_salmonella" }'             >> $out/$run_date/logs/samples.rm
ls $out/$run_date/mash_results/none/*/*gff                2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_mash" }'                    >> $out/$run_date/logs/samples.rm
ls $out/$run_date/mash_results/*/none/*gff                2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_mash" }'                    >> $out/$run_date/logs/samples.rm
ls $out/$run_date/ecoli/none/*/*gff                       2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_abricate_serotypefinder" }' >> $out/$run_date/logs/samples.rm
ls $out/$run_date/ecoli/*/none/*gff                       2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tno_abricate_serotypefinder" }' >> $out/$run_date/logs/samples.rm
ls $out/$run_date/ecoli/notecoli/*/*gff                   2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_ecoli" }'                  >> $out/$run_date/logs/samples.rm
ls $out/$run_date/ecoli/*/notecoli/*gff                   2>> $out/$run_date/logs/shortcut_overlap.txt | awk '{ print $0 "\tnot_ecoli" }'                  >> $out/$run_date/logs/samples.rm
ls $out/$run_date/mash_results/*/none/ -d                 2>> $out/$run_date/logs/shortcut_overlap.txt | parallel rm -R {}
ls $out/$run_date/ecoli/*/none/        -d                 2>> $out/$run_date/logs/shortcut_overlap.txt | parallel rm -R {}
ls $out/$run_date/ecoli/*/notecoli/    -d                 2>> $out/$run_date/logs/shortcut_overlap.txt | parallel rm -R {}
if [ -d "$out/$run_date/salmonella/enterica/no_result" ]      ; then rm -R $out/$run_date/salmonella/enterica/no_result      ; fi
if [ -d "$out/$run_date/salmonella/enterica/not_salmonella" ] ; then rm -R $out/$run_date/salmonella/enterica/not_salmonella ; fi
if [ -d "$out/$run_date/mash_results/none" ]                  ; then rm -R $out/$run_date/mash_results/none                  ; fi
if [ -d "$out/$run_date/ecoli/notecoli" ]                     ; then rm -R $out/$run_date/ecoli/notecoli                     ; fi
if [ -d "$out/$run_date/ecoli/none" ]                         ; then rm -R $out/$run_date/ecoli/none                         ; fi

# to avoid duplicates - should these be put into samples.rm?
if [ -d "$out/$run_date/salmonella/enterica" ] ; then rm -R $out/$run_date/mash_results/Salmonella/enterica ; fi
if [ -d "$out/$run_date/ecoli" ]               ; then rm -R $out/$run_date/mash_results/Escherichia/coli ; fi

# to remove files that are "too small"
ls $out/$run_date/*/*/*/ -d | parallel wc -l {}*gff | awk '{ if ( $1 < 10000) print $2 "\tsize_too_small"}' >> $out/$run_date/logs/samples.rm
ls $out/$run_date/*/*/*/ -d | parallel wc -l {}*gff | awk '{ if ( $1 < 10000) print $2 }' | parallel rm {}

date
echo "Removing empty directories, directories with less than 4 samples, and directories with more than $max_samples"
directories=($(ls -d $out/$run_date/*/*/*/ | sed s/.$// ))
for directory in ${directories[@]}
do
  if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
  then
    rm -R $directory
  else
    number_of_files=$(ls $directory/*gff | wc -l )
    if (( $number_of_files > $max_samples ))
    then
      ls $directory/*gff | awk '{ print $0 "\ttoo_many_samples" }' >> $out/$run_date/logs/samples.rm
      rm -R $directory
    elif (( $number_of_files < 4 ))
    then
      ls $directory/*gff | awk '{ print $0 "\ttoo_few_samples" }' >> $out/$run_date/logs/samples.rm
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
        if [ ! -f $directory/recent_samples.txt ]
        then
          ls $directory/*gff | awk '{ print $0 "\tno_recent" }' >> $out/$run_date/logs/samples.rm
          rm -R $directory
        fi
      fi
    fi
  fi
done
find $out/$run_date -type d -empty -delete

date
echo "Adding controls to samples"
mkdir -p $out/$run_date/salmonella/enterica/all
run_control $out $run_date $control_path
rm -R $out/$run_date/salmonella/enterica/all
find $out/$run_date -type d -empty -delete

date
echo "log files can be found at $out/$run_date/logs"
echo "Organizing gff files for OUTBREAK_120 is finished!"

exit
