#!/bin/bash

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

USAGE="
This is a script that takes the output from prokka, abricate, mash, and/or
SeqSero through roary, and iqtree in order to create a phylogenetic tree of
organisms with files created less than 120 days ago.

This current script is dated September 24, 2018 by Erin Young.

Usage: ./outbreak_120_organize.sh -i /files/location -o /output/directory/

A full list of options is available with
./outbreak_120_organize.sh -h
"
#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
HELP="
OUTBREAK_120 looks through a directory (or directories) for contig files.

REQUIRED OPTIONS:

-i <path to look through for input files>

OUTBREAK_120 will search this directory for files matching the following patterns:
Prokka files:     <path>/*/ALL_gff/NAME.gff
Mash files:       <path>/*/mash/NAME*.sorted.txt
SeqSero files:    <path>/*/SeqSero/NAME.Seqsero_result.txt
Abricate results: <path>/*/abricate_results/database/database.NAME.out.txt

-o <path for outfolder>

OPTIONAL OPTIONS:

-e <file>                 List of sample names to be removed from analysis if found.
-t mash/seqsero/abricate  Restrict results to specific analysis.
-p <species>              Restrict results to species in mash results.
-s <yyyy-mm-dd>           Specify end date.
-d <number of days>       Specify number of days prior to end date (default is 120)
-m <number of samples>    Adjust maximum number of samples (default is 100)
-a                        Evaluates all subfolders, not just those with samples less than 10 days old.
-P <PATH>                 Add path to PATH. Can be used to specify multiple paths.

"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

input_files=()
exclude_file=""
d=$(date +%Y-%m-%d)
run_date=$(date +%Y-%m-%d)
age=120
pd=$(date -d "$d - $age days" +%Y-%m-%d)
control_path=/home/Bioinformatics/Data/NCBI_refs/OUTBREAK_120_controls/ASSEMBLED_GENOMES
max_samples=100
serotyping=("mash" "seqsero" "abricate")
search_path="/home/IDGenomics_NAS/WGS_Serotyping"
out="/home/IDGenomics_NAS/WGS_Serotyping/OUTBREAK_120"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

run_clean.files ()
{
  out=$1
  d=$2
  max_samples=$3

  echo "Removing empty directories"
  directories=($(ls -d $out/$d/*/*/*/ ))
  for directory in ${directories[@]}
  do
    if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
    then
      rm -R $directory
      echo "$directory had no remaining gff files"
    else
      number_of_files=$(ls $directory/*gff | wc -l )
      if (( $number_of_files > $max_samples ))
      then
        ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\ttoo_many_samples" }' >> $out/$d/logs/samples.rm
        rm -R $directory
      elif (( $number_of_files < 4 ))
      then
        ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\ttoo_few_samples" }' >> $out/$d/logs/samples.rm
        rm -R $directory
      fi
    fi
  done
  salmo_directories=($(ls -d $out/$d/salmonella/*/ ))
  ecoli_directories=($(ls -d $out/$d/ecoli/*/ ))
  other_directories=($(ls -d $out/$d/other/*/ ))
  echo "These are the directories:"
  echo ${salmo_directories[@]} ${ecoli_directories[@]} ${other_directories[@]}
  for directory in ${salmo_directories[@]} ${ecoli_directories[@]} ${other_directories[@]} other salmonella ecoli
  do
    echo $directory
    if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
    then
      rm -R $directory
    fi
  done
}
run_clean.prune ()
{
  out=$1
  d=$2

  date
  echo "Removing E coli samples missing O or H groups"
  ls $out/$d/ecoli/none/*/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tmissing_O_group" }' >> $out/$d/logs/samples.rm
  ls $out/$d/ecoli/*/none/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tmissing_H_group" }' >> $out/$d/logs/samples.rm
  rm -R $out/$d/ecoli/none
  ls $out/$d/ecoli/*/none | parallel "rm -R {}"

  date
  echo "Removing samples in the salmonella/salmonella/NAThepredictedantigenicprofiledoesnotexistintheWhite-Kauffmann-LeMinorscheme folder"
  ls $out/$d/salmonella/salmonella/NAThepredictedantigenicprofiledoesnotexistintheWhite-Kauffmann-LeMinorscheme/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tNA_salmonella" }' >> $out/$d/logs/samples.rm
  rm -R $out/$d/salmonella/salmonella/NAThepredictedantigenicprofiledoesnotexistintheWhite-Kauffmann-LeMinorscheme

  if [ -d "$out/$d/other/Salmonella/enterica" ] && [ -d "$out/$d/salmonella/salmonella/all" ]
  then
    date
    echo "To avoid duplicates, Salmonella samples in $out/$d/other will be removed."
    ls $out/$d/other/Salmonella/enterica/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/logs/samples.rm
    ls $out/$d/other/Salmonella/all/*gff      | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/logs/samples.rm
    rm -R $out/$d/other/Salmonella/enterica
    rm -R $out/$d/other/Salmonella/all
  fi

  if [ -d "$out/$d/other/Escherichia/coli" ] && [ -d "$out/$d/ecoli/all/all" ]
  then
    date
    echo "To avoid duplicates, E. coli samples in $out/$d/other will be removed."
    ls $out/$d/other/Escherichia/coli/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/logs/samples.rm
    ls $out/$d/other/Escherichia/all/*gff  | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/logs/samples.rm
    rm -R $out/$d/other/Escherichia/coli
    rm -R $out/$d/other/Escherichia/all
  fi
}
run_clean.recent ()
{
  out=$1
  d=$2
  recent_flag=$3

  directories=($(ls -d $out/$d/*/*/*/ ))
  for directory in ${directories[@]}
  do
    ls $directory/*gff > $directory/list_of_samples.txt
    if [ -z "$recent_flag" ]
    then
      recent_check=$(grep -f $out/$d/files/recent_files.txt $directory/list_of_samples.txt | head -n 1 )
      if [ -z "$recent_check" ]
      then
        ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tno_recent" }' >> $out/$d/logs/samples.rm
        rm -R $directory
      fi
    fi
  done
}
run_clean.size ()
{
  out=$1
  d=$2

  date
  wc -l $out/$d/*/*/*/*gff | awk '{ if ( $1 < 10000) print $2 "\tsize_too_small"}' >> $out/$d/logs/samples.rm
  wc -l $out/$d/*/*/*/*gff | awk '{ if ( $1 < 10000) print $2 }' | xargs rm
}
run_control ()
{
  out=$1
  d=$2
  control_path=$3

  date
  echo "Adding NCBI genome controls for outgroups"
  control_list=(
  Acinetobacter_baumannii_FA_control.gff:other:Acinetobacter:baumannii
  Acinetobacter_nosocomialis_control.gff:other:Acinetobacter:nosocomialis
  Acinetobacter_pittii_FA_control.gff:other:Acinetobacter:pittii
  Bacteroides_fragilis_control.gff:other:Bacteroides:fragilis
  Campylobacter_coli_FA_control.gff:other:Campylobacter:coli
  Campylobacter_fetus_FA_control.gff:other:Campylobacter:fetus
  Campylobacter_hyointestinalis_FA_control.gff:other:Campylobacter:hyointestinalis
  Campylobacter_jejuni_FA_control.gff:other:Campylobacter:jejuni
  Campylobacter_lari_FA_control.gff:other:Campylobacter:lari
  Campylobacter_peloridis_control.gff:other:Campylobacter:peloridis
  Campylobacter_upsaliensis_FA_control.gff:other:Campylobacter:upsaliensis
  Clostridium_perfringens_FA_control.gff:other:Clostridium:perfringens
  Cronobacter_sakazakii_control.gff:other:Cronobacter:sakazakii
  Cronobacter_turicensis_control.gff:other:Cronobacter:turicensis
  Enterobacter_aerogenes_control.gff:other:Enterobacter:aerogenes
  Escherichia_albertii_control.gff:other:Escherichia:albertii
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
  Helicobacter_pullorum_control.gff:other:Helicobacter:pullorum
  Klebsiella_aerogenes_control.gff:other:Klebsiella:aerogenes
  Klebsiella_oxytoca_control.gff:other:Klebsiella:oxytoca
  Klebsiella_pneumoniae_control.gff:other:Klebsiella:pneumoniae
  Legionella_pneumophila_control.gff:other:Legionella:pneumophila
  Listeria_monocytogenes_FA_control.gff:other:Listeria:monocytogenes
  Pantoea_ananatis_control.gff:other:Pantoea:ananatis
  Pseudomonas_aeruginosa_control.gff:other:Pseudomonas:aeruginosa
  Salmonella_bongori_control.gff:other:Salmonella:bongori
  Salmonella_enterica_Aberdeen_control.gff:salmonella:salmonella:Aberdeen
  Salmonella_enterica_Abony_control.gff:salmonella:salmonella:Abony
  Salmonella_enterica_Adelaide_control.gff:salmonella:salmonella:Adelaide
  Salmonella_enterica_Agbeni_contigs_control.gff:salmonella:salmonella:BronorAgbeni
  Salmonella_enterica_Agona_control.gff:salmonella:salmonella:Agona
  Salmonella_enterica_Alachua_contigs_control.gff:salmonella:salmonella:Alachua
  Salmonella_enterica_Albany_control.gff:salmonella:salmonella:Albany
  Salmonella_enterica_Altona_contigs_control.gff:salmonella:salmonella:Altona
  Salmonella_enterica_Anatum_control.gff:salmonella:salmonella:Anatum
  Salmonella_enterica_Apapa_control.gff:salmonella:salmonella:Apapa
  Salmonella_enterica_Baildon_contigs_control.gff:salmonella:salmonella:Baildon
  Salmonella_enterica_Bareilly_control.gff:salmonella:salmonella:Bareilly
  Salmonella_enterica_Berta_contigs_control.gff:salmonella:salmonella:Berta
  Salmonella_enterica_Blockley_contigs_control.gff:salmonella:salmonella:HaardtorBlockley
  Salmonella_enterica_Bovismorbificans_control.gff:salmonella:salmonella:HindmarshorBovismorbificans
  Salmonella_enterica_Braenderup_control.gff:salmonella:salmonella:Braenderup
  Salmonella_enterica_Brandenburg_contigs_control.gff:salmonella:salmonella:Brandenburg
  Salmonella_enterica_Bredeney_control.gff:salmonella:salmonella:Bredeney
  Salmonella_enterica_Cerro_control.gff:salmonella:salmonella:Cerro
  Salmonella_enterica_Chester_control.gff:salmonella:salmonella:Chester
  Salmonella_enterica_Choleraesuis_control.gff:salmonella:salmonella:Choleraesuis
  Salmonella_enterica_Concord_control.gff:salmonella:salmonella:Concord
  Salmonella_enterica_Corvallis_control.gff:salmonella:salmonella:Corvallis
  Salmonella_enterica_Cubana_control.gff:salmonella:salmonella:Cubana
  Salmonella_enterica_Daytona_contigs_control.gff:salmonella:salmonella:Daytona
  Salmonella_enterica_Derby_control.gff:salmonella:salmonella:Derby
  Salmonella_enterica_Dublin_control.gff:salmonella:salmonella:Dublin
  Salmonella_enterica_Ealing_contigs_control.gff:salmonella:salmonella:Ealing
  Salmonella_enterica_Eastbourne_contigs_control.gff:salmonella:salmonella:Eastbourne
  Salmonella_enterica_Edinburgh_contigs_control.gff:salmonella:salmonella:Edinburgh
  Salmonella_enterica_Enteritidis_control.gff:salmonella:salmonella:Enteritidis
  Salmonella_enterica_Gaminara_control.gff:salmonella:salmonella:Gaminara
  Salmonella_enterica_Give_control.gff:salmonella:salmonella:Give
  Salmonella_enterica_Glostrup_contigs_control.gff:salmonella:salmonella:GlostruporChomedey
  Salmonella_enterica_Grumpensis_contigs_control.gff:salmonella:salmonella:Grumpensis
  Salmonella_enterica_Hadar_contigs_control.gff:salmonella:salmonella:HadarorIstanbul
  Salmonella_enterica_Hartford_contigs_control.gff:salmonella:salmonella:Hartford
  Salmonella_enterica_Havana_contigs_control.gff:salmonella:salmonella:Havana
  Salmonella_enterica_Heidelberg_control.gff:salmonella:salmonella:Heidelberg
  Salmonella_enterica_Hvittingfoss_control.gff:salmonella:salmonella:II16benxorHvittingfoss
  Salmonella_enterica_Ibadan_contigs_control.gff:salmonella:salmonella:IbadanorMississippiorII11323b15
  Salmonella_enterica_Indiana_control.gff:salmonella:salmonella:Indiana
  Salmonella_enterica_Infantis_control.gff:salmonella:salmonella:Infantis
  Salmonella_enterica_Istanbul_contigs_control.gff:salmonella:salmonella:HadarorIstanbul
  Salmonella_enterica_Javiana_control.gff:salmonella:salmonella:II912lz2815orJaviana
  Salmonella_enterica_Johannesburg_control.gff:salmonella:salmonella:Johannesburg
  Salmonella_enterica_Kentucky_control.gff:salmonella:salmonella:Kentucky
  Salmonella_enterica_Kiambu_contigs_control.gff:salmonella:salmonella:KiambuorII141227z15
  Salmonella_enterica_Kintambo_contigs_control.gff:salmonella:salmonella:KintamboorWashington
  Salmonella_enterica_Kottbus_contigs_control.gff:salmonella:salmonella:Kottbus
  Salmonella_enterica_Litchfield_contigs_control.gff:salmonella:salmonella:PakistanorLitchfield
  Salmonella_enterica_Livingstone_contigs_control.gff:salmonella:salmonella:Livingstone
  Salmonella_enterica_Lomalinda_contigs_control.gff:salmonella:salmonella:LomalindaorII1912aenx
  Salmonella_enterica_London_contigs_control.gff:salmonella:salmonella:London
  Salmonella_enterica_Manhattan_control.gff:salmonella:salmonella:YovokomeorManhattan
  Salmonella_enterica_Mbandaka_control.gff:salmonella:salmonella:Mbandaka
  Salmonella_enterica_Meleagridis_contigs_control.gff:salmonella:salmonella:Meleagridis
  Salmonella_enterica_Miami_contigs_control.gff:salmonella:salmonella:Miami
  Salmonella_enterica_Michigan_contigs_control.gff:salmonella:salmonella:Michigan
  Salmonella_enterica_Minnesota_control.gff:salmonella:salmonella:Minnesota
  Salmonella_enterica_Mississippi_contigs_control.gff:salmonella:salmonella:IbadanorMississippiorII11323b15
  Salmonella_enterica_Monschaui_contigs_control.gff:salmonella:salmonella:Monschaui
  Salmonella_enterica_Montevideo_control.gff:salmonella:salmonella:Montevideo
  Salmonella_enterica_Muenchen_contigs_control.gff:salmonella:salmonella:VirginiaorMuenchen
  Salmonella_enterica_Muenster_control.gff:salmonella:salmonella:Muenster
  Salmonella_enterica_Napoli_contigs_control.gff:salmonella:salmonella:Napoli
  Salmonella_enterica_Newport_control.gff:salmonella:salmonella:Newport
  Salmonella_enterica_Norwich_contigs_control.gff:salmonella:salmonella:Norwich
  Salmonella_enterica_Ohio_FA_control.gff:salmonella:salmonella:Ohio
  Salmonella_enterica_Onderstepoort_contigs_control.gff:salmonella:salmonella:BahrenfeldorOnderstepoort
  Salmonella_enterica_Oranienburg_control.gff:salmonella:salmonella:OranienburgorII67mt-
  Salmonella_enterica_Oslo_contigs_control.gff:salmonella:salmonella:Oslo
	Salmonella_enterica_Panama_control.gff:salmonella:salmonella:PanamaorHouston
	Salmonella_enterica_ParatyphiA_control.gff:salmonella:salmonella:ParatyphiA
	Salmonella_enterica_ParatyphiB_control.gff:salmonella:salmonella:ParatyphiB
	Salmonella_enterica_ParatyphiC_control.gff:salmonella:salmonella:ParatyphiC
	Salmonella_enterica_Pomona_control.gff:salmonella:salmonella:Pomona
	Salmonella_enterica_Poona_control.gff:salmonella:salmonella:FarmsenorPoona
	Salmonella_enterica_Potsdam_contigs_control.gff:salmonella:salmonella:Potsdam
	Salmonella_enterica_Reading_contigs_control.gff:salmonella:salmonella:Reading
	Salmonella_enterica_Rissen_contigs_control.gff:salmonella:salmonella:Rissen
	Salmonella_enterica_Rubislaw_control.gff:salmonella:salmonella:Rubislaw
	Salmonella_enterica_Saintpaul_control.gff:salmonella:salmonella:Saintpaul
	Salmonella_enterica_Sandiego_contigs_control.gff:salmonella:salmonella:Sandiego
	Salmonella_enterica_Schwarzengrund_control.gff:salmonella:salmonella:Schwarzengrund
	Salmonella_enterica_Senftenberg_control.gff:salmonella:salmonella:SenftenbergorDessau
	Salmonella_enterica_Singapore_contigs_control.gff:salmonella:salmonella:Singapore
	Salmonella_enterica_Stanley_control.gff:salmonella:salmonella:Stanley
	Salmonella_enterica_Stanleyville_control.gff:salmonella:salmonella:Stanleyville
	Salmonella_enterica_Sundsvall_contigs_control.gff:salmonella:salmonella:Sundsvall
	Salmonella_enterica_Tallahassee_contigs_control.gff:salmonella:salmonella:Tallahassee
	Salmonella_enterica_Telelkebir_contigs_control.gff:salmonella:salmonella:Telelkebir
	Salmonella_enterica_Tennessee_control.gff:salmonella:salmonella:II67z29z42orTennessee
	Salmonella_enterica_Thompson_control.gff:salmonella:salmonella:Thompson
	Salmonella_enterica_Typhi_FA_control.gff:salmonella:salmonella:Typhi
	Salmonella_enterica_Typhimurium_FA_control.gff:salmonella:salmonella:Typhimurium
	Salmonella_enterica_Uganda_contigs_control.gff:salmonella:salmonella:Uganda
	Salmonella_enterica_Urbana_contigs_control.gff:salmonella:salmonella:Urbana
	Salmonella_enterica_Virchow_contigs_control.gff:salmonella:salmonella:Virchow
	Salmonella_enterica_Wandsworth_control.gff:salmonella:salmonella:Wandsworth
	Salmonella_enterica_Weltevreden_control.gff:salmonella:salmonella:Weltevreden
	Salmonella_enterica_Worthington_contigs_control.gff:salmonella:salmonella:Worthington
	Serratia_marcescens_control.gff:other:Serratia:marcescens
	Shigella_flexneri_control.gff:other:Shigella:flexneri
	Shigella_sonnei_FA_control.gff:other:Shigella:sonnei
	Variovorax_paradoxus_control.gff:other:Variovorax:paradoxus
	Vibrio_navarrensis_control.gff:other:Vibrio:navarrensis
	Vibrio_parahaemolyticus_control.gff:other:Vibrio:parahaemolyticus
  )

  for control in ${control_list[@]}
  do
    control_parts=($(echo "$control_path/$control" | awk -F ":" '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' ))
    if [ -d "$out/$d/${control_parts[1]}/all/all" ] ; then ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/all/all/. ; fi

    if [ -d "$out/$d/${control_parts[1]}/${control_parts[2]}/all" ]
    then
      echo "Adding ${control_parts[0]}"
      ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/${control_parts[2]}/all/.

      if [ -d "$out/$d/other/Salmonella/enterica" ] && [ "${control_parts[2]}" == "salmonella" ]
      then
        ln -s ${control_parts[0]} $out/$d/other/Salmonella/enterica/.
        ln -s ${control_parts[0]} $out/$d/other/Salmonella/all/.
      elif [ -d "$out/$d/other/Escherichia/coli" ] && [ "${control_parts[2]}" == "ecoli" ]
      then
        ln -s ${control_parts[0]} $out/$d/other/Escherichia/coli/.
        ln -s ${control_parts[0]} $out/$d/other/Escherichia/all/.
      fi

      if [ -d "$out/$d/${control_parts[1]}/${control_parts[2]}/${control_parts[3]}" ]
      then
        ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/${control_parts[2]}/${control_parts[3]}/.
      fi
    fi
  done

  date
  echo "Adding controls has completed."
}
run_organize.abricate ()
{
  gff_file=$1
  sample=$2
  out=$3
  d=$4

  mkdir -p $out/$d/ecoli
  mkdir -p $out/$d/ecoli/all
  mkdir -p $out/$d/ecoli/all/all
  if [ -n "$(find $out/$d/serotyping_results/abricate -iname serotypefinder*$sample*out.tab | head -n 1 )" ]
  then
    abricate_file=$(ls $out/$d/serotyping_results/abricate/serotypefinder*$sample*out*tab | head -n 1 )
    O_group=($(cat $abricate_file | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "O" | sed 's/\///g' ))
    H_group=($(cat $abricate_file | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "H" | sed 's/\///g' ))
    if [ -n "${O_group[0]}" ] && [ -n "${H_group[0]}" ] ; then ln -s $gff_file $out/$d/ecoli/all/all/. ; fi
    if [ -z "${O_group[0]}" ] ; then O_group=("none") ; fi
    if [ -z "${H_group[0]}" ] ; then H_group=("none") ; fi
    echo "Abricate results for $sample: O group: ${O_group[@]} H group: ${H_group[@]}"

    for O in ${O_group[@]}
    do
      mkdir -p $out/$d/ecoli/$O
      mkdir -p $out/$d/ecoli/$O/all
      ln -s $gff_file $out/$d/ecoli/$O/all/.
      for H in ${H_group[@]}
      do
        mkdir -p $out/$d/ecoli/$O/$H
        ln -s $gff_file $out/$d/ecoli/$O/$H/.
      done
    done
  else
    echo -e "$sample\tno_abricate_serotypefinder" >> $out/$d/logs/samples.rm
  fi
}
run_organize.mash ()
{
  gff_file=$1
  sample=$2
  out=$3
  d=$4
  species=$5

  if [ -z "$species" ] ; then species="." ; fi
  mkdir -p $out/$d/other
  if [ -n "$(find $out/$d/serotyping_results/mash -iname $sample*sorted.txt | head -n 1 )" ]
  then
    mash_file=$(ls $out/$d/serotyping_results/mash/$sample*sorted.txt | head -n 1 )
    other_serotype=($(head -n 1 $mash_file | grep "$species" | cut -f 1 | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "\t" $2 }'))
    echo "Mash results for $sample: ${other_serotype[@]}"

    if [ -n "${other_serotype[0]}" ] && [ -n "${other_serotype[1]}" ]
    then
      mkdir -p $out/$d/other/${other_serotype[0]}
      mkdir -p $out/$d/other/${other_serotype[0]}/all
      mkdir -p $out/$d/other/${other_serotype[0]}/${other_serotype[1]}
      ln -s $gff_file $out/$d/other/${other_serotype[0]}/all/.
      ln -s $gff_file $out/$d/other/${other_serotype[0]}/${other_serotype[1]}/.
    else
      echo -e "$sample\tno_mash" >> $out/$d/logs/samples.rm
    fi
  else
    echo -e "$sample\tno_mash" >> $out/$d/logs/samples.rm
  fi
}
run_organize.seqsero ()
{
  gff_file=$1
  sample=$2
  out=$3
  d=$4

  mkdir -p $out/$d/salmonella
  mkdir -p $out/$d/salmonella/salmonella
  mkdir -p $out/$d/salmonella/salmonella/all
  if [ -n "$(find $out/$d/serotyping_results/seqsero -iname $sample*Seqsero_result.txt | head -n 1 )" ]
  then
    seqsero_file=$(ls $out/$d/serotyping_results/seqsero/$sample*Seqsero_result.txt | head -n 1 )
    seqsero_serotype=$(grep "Predicted serotype(s):" $seqsero_file | awk '{ $1=$2="" ; print $0 }' | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/O5-//g')
    echo "Seqsero results for $sample: $seqsero_serotype"

    mkdir -p $out/$d/salmonella/salmonella/$seqsero_serotype
    ln -s $gff_file $out/$d/salmonella/salmonella/$seqsero_serotype/.
    ln -s $gff_file $out/$d/salmonella/salmonella/all/.
  else
    echo -e "$sample\tno_seqsero" >> $out/$d/logs/samples.rm
  fi
}

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
while getopts 'ac:d:e:f:hi:m:o:p:s:t:P:' OPTION
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
    d)
    echo "Looking for files that have been modified $OPTARG days prior"
    age=$OPTARG
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
    echo "Input file is located at $OPTARG"
    input_files=("${input_files[@]}" "$OPTARG")
    if [ ! -f "$OPTARG" ]
    then
      echo "$OPTARG could not be located"
      echo $USAGE
      exit
    fi
    ;;
    h)
    echo "$HELP"
    exit 1
    ;;
    i)
    echo "Creating an input file from the following path: $OPTARG"
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
    serotyping=("mash")
    echo "Will look through mash results for samples that are $species"
    ;;
    s)
    echo "Ending date is now $OPTARG instead of $d"
    d=$OPTARG
    ;;
    t)
    if [[ "$OPTARG" == "mash" || "$OPTARG" == "abricate" || "$OPTARG" == "seqsero" ]]
    then
      echo "Results are restricted to those found in $OPTARG"
      serotyping=("$OPTARG")
    else
      echo "Invalid option: $OPTARG can only accept mash/seqsero/abricate"
      echo "$USAGE"
      exit 1
    fi
    ;;
    P)
    echo "Adding $OPTARG to PATH"
    PATH=$PATH:$OPTARG
    echo "PATH is now $PATH"
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

if [ -z "$out" ] ; then echo "Output directory is required" ; echo "$USAGE" ; exit ; fi
mkdir -p $out
mkdir -p $out/$d
if [ ! -d "$out/$d" ] ; then echo "Could not create directory for results!" ; exit ; fi
mkdir -p $out/$d/files
mkdir -p $out/$d/logs
pd=$(date -d "$d - $age days" +%Y-%m-%d)
rd=$(date -d "$d - 10 days" +%Y-%m-%d)
if [ ! -f "$out/$d/logs/samples.rm" ] ; then echo -e "sample\treason" > $out/$d/logs/samples.rm ; fi

if [ -d "$search_path" ]
then
  mkdir -p $out/$d/serotyping_results
  mkdir -p $out/$d/serotyping_results/abricate
  mkdir -p $out/$d/serotyping_results/seqsero
  mkdir -p $out/$d/serotyping_results/mash

  list_of_files="$out/$d/files/all_gff_files.txt"
  list_of_recnt="$out/$d/files/recent_gff_files.txt"

  find $search_path/*/ALL_gff -iname *gff -newermt $pd ! -newermt $d -type f > $list_of_files
  find $search_path/*/ALL_gff -iname *gff -newermt $rd ! -newermt $d -type f > $list_of_recnt
  if [ -f "$exclude_file" ]
  then
    grep -v -f $exclude_file $out/$d/files/all_gff_files.txt > $out/$d/files/filtered_gff_files.txt
    grep -v -f $exclude_file $out/$d/files/recent_gff_files.txt > $out/$d/files/filtered_recent_gff_files.txt
    list_of_files="$out/$d/files/filtered_gff_files.txt"
    list_of_recnt="$out/$d/files/filtered_recent_gff_files.txt"
  fi

  cut $list_of_recnt -f 7 -d "/" > $out/$d/files/recent_files.txt

  date
#  for database in argannot card ecoh ecoli ecoli_vf ncbi plasmidfinder resfinder serotypefinder vfdb cpd
  for database in ncbi serotypefinder
  do
    echo "Gathering abricate results for $database"
    ln -s $search_path/*/abricate_results/$database/*.out.tab $out/$d/serotyping_results/abricate/.
    ln -s $search_path/*/ALL_assembled_plasmids/$database/*.out.tab $out/$d/serotyping_results/abricate/.
  done
  if [ -n "$(echo ${serotyping[@]} | grep seqsero)" ]
  then
    date
    echo "Gathering seqsero results"
    ln -s $search_path/*/SeqSero/*Seqsero_result.txt $out/$d/serotyping_results/seqsero/.
  fi
  if [ -n "$(echo ${serotyping[@]} | grep mash)" ]
  then
    date
    echo "Gathering mash results"
    ln -s $search_path/*/mash/*sorted.txt $out/$d/serotyping_results/mash/.
  fi

  while read gff_file
  do
    date
    sample=$(echo $gff_file | cut -f 7 -d '/' | sed 's/.gff//g' | sort | uniq )
    echo "Finding serotying files for $sample"
    for kind in ${serotyping[@]}
    do
      run_organize.$kind "$gff_file" "$sample" "$out" "$d" "$species"
    done
  done < $list_of_files
fi

date
echo "Removing samples with files that are too small"
run_clean.size $out $d

echo "recent flag is $recent_flag"
date
echo "Listing samples in directory"
run_clean.recent $out $d $recent_flag

date
echo "Removing samples with abnormal serotypes"
run_clean.prune $out $d

date
echo "Adding controls to samples"
run_control $out $d $control_path

date
echo "Removing directories with too few or too many samples"
run_clean.files $out $d $max_samples

date
echo "log files can be found at $out/$d/logs"
echo "OUTBREAK_120 has completed!"
