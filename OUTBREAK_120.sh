#!/bin/bash

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

USAGE="
This is a script that takes the output from prokka, abricate, mash, and/or
SeqSero through roary, and iqtree in order to create a phylogenetic tree of
organisms with files created less than 120 days ago.

This current script is dated September 24, 2018 by Erin Young.

Usage: ./OUTBREAK_120.sh -f /input/file -o /output/directory/

A full list of options is available with
./OUTBREAK_120.sh -h
"
#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
HELP="
OUTBREAK_120 looks through a directory (or directories) for contig files.

REQUIRED OPTIONS:

-f <path to input file(s)>

Full path to input file(s). Input file is tab deliminated file with columns
\"date \t sample_name \t path/to/prokka/file \t type_of_serotyping \t /path/to/serotyping\"
Where type_of_serotyping is either MASH, SeqSero (Salmonella only) or ABRICATE
(E. coli only) and the date is yyyy-mm-dd format.

Samples that must be included for the analysis should have columns
\"include \t path/to/prokka/file \t type_of_serotyping \t /path/to/serotyping\"

-o <path for outfolder>

This path is where the results and subfolders for this script will be created.

OPTIONAL OPTIONS:

-i <path to look through for input files>

Instead of providing an input file with \"-f\", create one based on file
modification time signatures. OUTBREAK_120 will search this directory for
files matching the following patterns:
Prokka files:     <path>/*/ALL_gff/NAME.gff
Mash files:       <path>/*/mash/NAME.msh.distance.txt.sorted.txt
SeqSero files:    <path>/*/SeqSero/NAME.Seqsero_result.txt
Abricate results: <path>/*/abricate_results/Ecoli_serotypefinder_out.txt

As of right now, there is no way to update this. If your files do not meet these
patterns, you will need to create your own script to create an input file.

-e <names samples to exclude, even if younger than 120 days>

List of sample names to be removed from analysis if found.

-t mash/seqsero/abricate

Choose between results from mash, SeqSero (Salmonella only), or Abricate
(E. coli only) for the final phylogenetic trees.

-p <string specifiying what species to look for>

Looks for samples matching the string in mash results. No spaces.

-s <yyyy-mm-dd>

Specify start date. Format must be in yyyy-mm-dd.

-d <number of days>

The default is 120 days prior to the start date.

-m <number of samples>

Adjust maximum number of samples. The default is 100 samples.

-x

Runs raxmlHPC-SSE3 instead of IQTREE.

-q

Runs both IQTREE and raxmlHPC-SSE3.

-a

Evaluates all subfolders, not just those samples less than 10 days old.

-v

Display versions and exit.

-I

Does not follow full pipeline. Just searches through path and creates an input
file.

-O

Does not follow full pipeline. Just uses input file(s) to create subfolders and
organizes samples for future Roary and IQTREE analyses.

-E

Does not follow full pipeline. Removes samples identified via an exclude file
specified with \"-x\"

-R

Does not follow full pipeline. Runs roary on gff files in subfolders specified in
\"-o\" option. Warning: Requires 5> number of samples>100

-Q

Does not follow full pipeline. Runs IQTREE on roary core alignment files in
subfolders specified in \"-o\" option.

-X

Does not follow full pipeline. Runs RAXML on roary core alignment files in
subfolders specified in \"-o\" option.

-A

Does not follow full pipeline. Runs organizes abricate results for GGTREE heatmaps.
ABRICATE results should be in the <outfolder>/<date>/all_abricate_results.txt

-G

Does not follow full pipeline. Runs GGTREE Rscript on iqtree treefile as well as the
core alignment fasta (core_gene_alignment.aln) and corresponding gene Rtab file.

-V

Check path for required programs.

-P

Add path to PATH. Can be used to specify multiple paths.

"

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

input_files=()
exclude_file=""
d=$(date +%Y-%m-%d)
run_date=$(date +%Y-%m-%d)
age=120
pd=$(date -d "$d - $age days" +%Y-%m-%d)
control_path=/home/Bioinformatics/Data/NCBI_refs/OUTBREAK_120_controls/ASSEMBLED_GENOMES
kraken_path=/home/IDGenomics_NAS/kraken_mini_db/minikraken_20141208
max_samples=100
search_flag=""
organize_flag=""
clean_flag=""
control_flag=""
roary_flag=""
iqtree_flag=""
abricate_flag=""
ggtree_flag=""
recent_flag=""
both_flag=""
raxml_flag=""

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

run_abricate ()
{
    out=$1
    d=$2

    date
    if [ -z "$(find $out/$d/ -iname "*_abricate_results.txt") " ]
    then
	echo "Coult not find $out/$d/*_abricate_results.txt"
	exit
    fi

    echo "Parsing ABRICATE results:"
    echo "abricate --summary directory/ABRICATE/abricate_database/*tab > directory/ABRICATE/abricate_database/summary.txt"

    directories=($(ls -d $out/$d/*/*/*/ ))
    for directory in ${directories[@]}
    do
	abricate_databases=($(ls $out/$d/*.abricate_results.txt | awk -F "/" '{ print $NF }' | sed 's/.abricate_results.txt//g' ))
	for abricate_database in ${abricate_databases[@]}
	do
	    if [ -f "$directory/IQTREE/$d.iqtree.iqtree" ] || [ -f "$directory/RAXML/RAxML_bipartitionsBranchLabels.raxml" ]
	    then
		date
		mkdir -p $directory/ABRICATE
		mkdir -p $directory/ABRICATE/$abricate_database

	    	ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | sort | uniq | parallel --jobs 10 "echo -e '#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT' > $directory/ABRICATE/$abricate_database/{}.tab"
		ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | sort | uniq | parallel --jobs 10 "grep {} $out/$d/$abricate_database.abricate_results.txt >> $directory/ABRICATE/$abricate_database/{}.tab"
		abricate --summary $directory/ABRICATE/$abricate_database/*tab > $directory/ABRICATE/$abricate_database/$abricate_database.summary.txt
		cat $directory/ABRICATE/$abricate_database/$abricate_database.summary.txt | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }' | sed 's/\t/,/g' | sed 's/ /,/g' | sed 's/[.],/0,/g' | sed 's/,[.]/,0/g' | sed 's/,,/,/g' > $out/ABRICATE/$abricate_database/$abricate_database.summary.csv
	    else
		echo "Not pursuing abricate results for $directory"
	    fi
	done
    done

    date
    echo "Combined ABRICATE gene presence and absence results:"
    ls $out/$d/*/*/*/ABRICATE/*/*summary.csv
    echo "ABRICATE is complete!"
}
run_clean.exclude ()
{
    out=$1
    d=$2
    exclude_file=$3

    date
    if [ -f "$exclude_file" ]
    then
	echo "Ensuring that samples are not listed in $exclude_file"
	while read exclusion
	do
	    if [ -n "$(find $out/$d/*/*/* -name *$exclusion*gff | head -n 1 )" ]
	    then
		ls $out/$d/*/*/*/$exclusion*gff | awk '{ print $0 "\texclude_file" }' >> $out/$d/$d.samples.rm
		rm $out/$d/*/*/*/$exclusion*gff
	    else
		echo "$exclusion was not found among gff files"
	    fi
	done < $exclude_file
    else
	echo "$exclude_file could not be found!"
    fi
    echo "All samples listed in $exclude_file were removed"
}
run_clean.prune ()
{
    out=$1
    d=$2

    date
    echo "Removing empty directories"
    directories=($(ls -d $out/$d/*/*/*))
    for directory in ${directories[@]}
    do
	if [ -z "$(find $directory -iname '*gff' | head -n 1 )" ]
	then
	    rm -R $directory
	fi
    done

    echo "Removing E coli samples missing O or H groups"
    ls $out/$d/ecoli/none/*/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tmissing_O_group" }' >> $out/$d/$d.samples.rm
    ls $out/$d/ecoli/*/none/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tmissing_H_group" }' >> $out/$d/$d.samples.rm
    rm -R $out/$d/ecoli/none
    ls $out/$d/ecoli/*/none | parallel "rm -R {}"

    if [ -d "$out/$d/other/Salmonella/enterica" ] && [ -d "$out/$d/salmonella/salmonella/all" ]
    then
	date
	echo "To avoid duplicates, Salmonella samples in $out/$d/other will be removed."
	ls $out/$d/other/Salmonella/enterica/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/$d.samples.rm
	ls $out/$d/other/Salmonella/all/*gff      | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/$d.samples.rm
	rm -R $out/$d/other/Salmonella/enterica
	rm -R $out/$d/other/Salmonella/all
    fi

    if [ -d "$out/$d/other/Escherichia/coli" ] && [ -d "$out/$d/ecoli/all/all" ]
    then
	date
	echo "To avoid duplicates, E. coli samples in $out/$d/other will be removed."
	ls $out/$d/other/Escherichia/coli/*gff | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/$d.samples.rm
	ls $out/$d/other/Escherichia/all/*gff  | sed 's!.*/!!' | cut -d "." -f 1 | awk '{ print $0 "\tduplicate" }' >> $out/$d/$d.samples.rm
	rm -R $out/$d/other/Escherichia/coli
	rm -R $out/$d/other/Escherichia/all
    fi

    date
    echo "A list of the samples removed can be found at $out/$d/$d.samples.rm"
}
run_clean.size ()
{
    out=$1
    d=$2

    date
    if [ -z "$(find $out/$d/*/*/* -name *gff | head -n 1 )" ]
    then
	echo "Cannot find necessary files."
	exit
    fi

    echo "Going through samples and ensuring adequate size"
    wc -l $out/$d/*/*/*/*gff | awk '{ if ( $1 < 10000) print $2 "\tsize_too_small"}' >> $out/$d/$d.samples.rm
    wc -l $out/$d/*/*/*/*gff | awk '{ if ( $1 < 10000) print $2 }' | xargs rm

    date
    echo "Removing samples based on size is complete"
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
	echo "Adding ${control_parts[0]}"
	if [ -d "$out/$d/${control_parts[1]}/all/all" ]
	then
	    ln -s ${control_parts[0]} $out/$d/${control_parts[1]}/all/all/.
	fi

	if [ -d "$out/$d/${control_parts[1]}/${control_parts[2]}/all" ]
	then
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
run_dependencies ()
{
    for program in roary iqtree parallel awk R abricate raxmlHPC-SSE3
    do
	if [ -z "$(which $program)" ]
	then
	    depend_check=1
	    echo "FAILURE: $program was not found"
	else
	    which $program
	fi
    done

    #ggtree

    if [ -n "$depend_check" ]
    then
	echo "OUTBREAK_120.sh cannot run. Please put these in your path!"
	exit
    fi
}
run_ggtree ()
{
    out=$1
    d=$2
    date
    if [ -z "$(find $out/$d/*/*/*/Roary_out -iname core_gene_alignment.aln | head -n 1 )" ] || [ -z "$(find $out/$d/*/*/*/Roary_out gene_presence_absence.Rtab | head -n 1 )" ]
    then
	if [ -z "$(find $out/$d/*/*/*/IQTREE -iname $d.iqtree.iqtree | head -n 1 )" ] && [ -z "$(find $out/$d/*/*/*/RAXML -iname RAxML_bipartitionsBranchLabels.raxml | head -n 1 )" ]
	then
	    echo "Cannot find required files"
	    exit
	fi
    fi

    echo "Running GGTREE on treefiles:"
    echo "Rscript PLOTS_OUTBREAK_120.R"

    cd $out/$d
    file_directories=($(ls -d $out/$d/*/*/* ))
    for directory in ${file_directories[@]}
    do
	if [ -f "$directory/Roary_out/gene_presence_absence.Rtab" ] && [ -f "$directory/Roary_out/core_gene_alignment.aln" ] && [ -f "$directory/ABRICATE/abricate_presence_absence.tab.R" ]
	then
	    cp /home/eyoung/sandbox/OUTBREAK_120/PLOTS_OUTBREAK_120.R $directory/.
	    split_paths=()
	    split_paths=($(echo "$directory" | awk -F "/" '{print $(NF-3) "\t" $(NF-2) "\t" $(NF-1) "\t" $(NF) }' ))
	    PRESABSC=$(echo "$directory/Roary_out/gene_presence_absence.Rtab" | sed 's/\//\\\//g' )
	    ROARYALN=$(echo "$directory/Roary_out/core_gene_alignment.aln" | sed 's/\//\\\//g' )
	    ABRICALN=$(echo "$directory/ABRICATE/abricate_presence_absence.tab.R" | sed 's/\//\\\//g' )
	    ABSENIMG=$(echo "$directory/${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}.roary_gene_presence" | sed 's/\//\\\//g' )
	    DISTAIMG=$(echo "$directory/${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}.nucleotide_distance" | sed 's/\//\\\//g' )
	    ABRICIMG=$(echo "$directory/${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}.abricate_resistance" | sed 's/\//\\\//g' )
	    TREEFILE="TREEFILE"
	    RXMLFILE="RXMLFILE"
	    if [ -f "$directory/IQTREE/$d.iqtree.treefile" ]
	    then
		TREEFILE=$(echo "$directory/IQTREE/$d.iqtree.treefile" | sed 's/\//\\\//g' )
		perl -pi -e "s/##IQTREE//g" $directory/PLOTS_OUTBREAK_120.R
	    fi
	    if [ -f "$directory/RAXML/RAxML_bipartitionsBranchLabels.raxml" ]
	    then
		RXMLFILE=$(echo "$directory/RAXML/RAxML_bipartitionsBranchLabels.raxml" | sed 's/\//\\\//g' )
		perl -pi -e "s/##RAXML//g" $directory/PLOTS_OUTBREAK_120.R
	    fi
	    perl -pi -e "s/FULLPATHTOIQTREETREEFILE/$TREEFILE/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTORAXML_TREEFILE/$RXMLFILE/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTOROARYGENETABLE/$PRESABSC/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTOROARYALIGNMENT/$ROARYALN/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTOABRICATE_TABLE/$ABRICALN/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTODISTANCEIMAGE/$DISTAIMG/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTOGENETABLIMAGE/$ABSENIMG/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/FULLPATHTORESISTENIMAGE/$ABRICIMG/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/TITLEOFTREE/${split_paths[1]}.${split_paths[2]}.${split_paths[3]}/g" $directory/PLOTS_OUTBREAK_120.R
	    perl -pi -e "s/DATEOFRUN/$d/g" $directory/PLOTS_OUTBREAK_120.R
	else
	    echo "Could not find required files! for $directory"
	fi
    done

    ls $out/$d/*/*/*/PLOTS_OUTBREAK_120.R | parallel --jobs 25 "Rscript {}"
    mkdir -p $out/$d/results
    cp $out/$d/*/*/*/*tiff $out/$d/results/.
    cp $out/$d/*/*/*/*pdf  $out/$d/results/.

    date
    echo "GGTREE has created these tiff files:"
    ls $out/$d/results/*
    echo "Which are also located in $out/$d"
    echo "GGTREE is complete!"
}
run_iqtree ()
{
    out=$1
    d=$2
    date
    if [ -z "$(find $out/$d/*/*/*/Roary_out -iname core_gene_alignment.aln -type f | head -n 1 )" ]
    then
	echo "Cannot find required files"
	exit
    fi

    echo "Running IQtree on core alignment files:"
    echo "iqtree -s directory/Roary_out/core_gene_alignment.aln -t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -pre directory/IQTREE/$d.iqtree -nt AUTO -o control"

    iqtree_array=()

    ls -d $out/$d/*/*/*/Roary_out | sed 's/Roary_out//g' | parallel " mkdir -p {}/IQTREE "
    roary_directories=($(ls -d $out/$d/*/*/*))
    for directory in ${roary_directories[@]}
    do
	if [ -f "$directory/Roary_out/core_gene_alignment.aln" ]
	then
	    if [ -n "$(find $directory -iname '*_control.gff' | head -n 1 )" ]
	    then
		outgroup=$(ls $directory/*_control.gff | head -n 1 )
		outgroup_basename=$(basename "$outgroup" .gff)
		iqtree_array+=("$(echo "iqtree___-s___$directory/Roary_out/core_gene_alignment.aln___-t___RANDOM___-m___GTR+F+I___-bb___1000___-alrt___1000___-pre___$directory/IQTREE/$d.iqtree___-nt___AUTO___-o___$outgroup_basename " )")
	    else
		iqtree_array+=("$(echo "iqtree___-s___$directory/Roary_out/core_gene_alignment.aln___-t___RANDOM___-m___GTR+F+I___-bb___1000___-alrt___1000___-pre___$directory/IQTREE/$d.iqtree___-nt___AUTO" )")
	    fi
	fi
    done

    history -p ${iqtree_array[@]} | sed 's/___/ /g' | parallel --jobs 10 '{}'

    date
    echo "IQTREE has created these iqtree files:"
    ls $out/$d/*/*/*/IQTREE/*iqtree
    echo "IQTREE is complete!"
}
run_organize.abricate ()
{
    out=$1
    d=$2
    pd=$3
    input_file=$4

    date
    echo "Organizing E. coli gff files based off of abricate O and H groups"
    if [ ! -d "$out/$d" ] || [ ! -f "$input_file" ]
    then
	echo "Required directories and files were not found."
	exit
    fi

    mkdir -p $out/$d/ecoli
    mkdir -p $out/$d/ecoli/all
    mkdir -p $out/$d/ecoli/all/all
    if [ ! -f "$out/$d/ecoli/ecoli_samples.txt" ]
    then
	echo -e "date\tsample\tprokka_file\tserotyping_file_type\tserotyping_file" > $out/$d/ecoli/ecoli_samples.txt
    fi

    while read ecoli
    do
	ecoli_files=($(echo $ecoli | awk '{ print $2 "\t" $3 "\t" $5 }' ))
	ecoli_O_group=($(grep ${ecoli_files[0]} ${ecoli_files[2]} | cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "O" | sed 's/\///g' ))
	ecoli_H_group=($(grep ${ecoli_files[0]} ${ecoli_files[2]} | cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "H" | sed 's/\///g' ))

	if [ -z "${ecoli_O_group[0]}" ] && [ -z "${ecoli_H_group[0]}" ]
	then
	    ecoli_O_group=("none")
	    ecoli_H_group=("none")
	elif [ -z "${ecoli_O_group[0]}" ] && [ -n "${ecoli_H_group[0]}" ]
	then
	    ecoli_O_group=("none")
	elif [ -n "${ecoli_O_group[0]}" ] && [ -z "${ecoli_H_group[0]}" ]
	then
	    ecoli_H_group=("none")
	else
	    ln -s ${ecoli_files[1]} $out/$d/ecoli/all/all/.
	fi

	echo -e "$ecoli\tO:${ecoli_O_group[@]};H:${ecoli_H_group[@]}" >> $out/$d/ecoli/ecoli_samples.txt

	for O_group in ${ecoli_O_group[@]}
	do
	    mkdir -p $out/$d/ecoli/$O_group
	    mkdir -p $out/$d/ecoli/$O_group/all
	    ln -s ${ecoli_files[1]} $out/$d/ecoli/$O_group/all/.
	    for H_group in ${ecoli_H_group[@]}
	    do
		mkdir -p $out/$d/ecoli/$O_group/$H_group
		ln -s ${ecoli_files[1]} $out/$d/ecoli/$O_group/$H_group/.
	    done
	done
    done < <(grep -w "abricate" $input_file | awk -v pd=$pd '{ if ( $1 > pd || $1=="include" ) print $0 }' )

    date
    if [ -f "$out/$d/ecoli/ecoli_samples.txt" ]
    then
	echo "$out/$d/ecoli/ecoli_samples.txt has been created:"
	echo " "
	echo "File organization is for samples with abricate results is complete!"
    else
	echo "Something went wrong and file were not moved appropriately. Sorry!"
    fi
}
run_organize.mash ()
{
    out=$1
    d=$2
    pd=$3
    input_file=$4

    date
    echo "Organizing gff files on species identified as the top hit in MASH"
    if [ ! -d "$out/$d" ] || [ ! -f "$input_file" ]
    then
	echo "Required directories and files were not found."
	exit
    fi

    mkdir -p $out/$d/other
    mkdir -p $out/$d/other/all
    mkdir -p $out/$d/other/all/all
    if [ ! -f "$out/$d/other/mash_samples.txt" ]
    then
	echo -e "date\tsample\tprokka_file\tserotyping_file_type\tserotyping_file" > $out/$d/other/mash_samples.txt
    fi

    while read other
    do
	other_serotype=($(echo $other | awk '{ print $5 }' | parallel " head -n 1 {} " | cut -f 1 | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "\t" $2 }'))
	mkdir -p $out/$d/other/${other_serotype[0]}
	mkdir -p $out/$d/other/${other_serotype[0]}/all
	mkdir -p $out/$d/other/${other_serotype[0]}/${other_serotype[1]}
	echo -e "$other\t${other_serotype[0]}""_""${other_serotype[1]}" >> $out/$d/other/mash_samples.txt
	other_file=$(echo $other | awk '{ print $3 }')
	ln -s $other_file $out/$d/other/${other_serotype[0]}/all/.
	ln -s $other_file $out/$d/other/${other_serotype[0]}/${other_serotype[1]}/.
    done < <(grep -w "mash" $input_file | awk -v pd=$pd '{ if ( $1 > pd || $1=="include" ) print $0 }' )

    date
    if [ -f "$out/$d/other/mash_samples.txt" ]
    then
	echo "$out/$d/other/mash_samples.txt has been created"
	echo " "
	echo "File organization is for samples with mash results is complete!"
    else
	echo "Something went wrong and file were not moved appropriately. Sorry!"
    fi
}
run_organize.mash.subtype ()
{
    out=$1
    d=$2
    pd=$3
    input_file=$4
    species=$5

    date
    echo "Organizing gff files on $species identified as the top hit in MASH"
    if [ ! -d "$out/$d" ] || [ ! -f "$input_file" ] || [ -z "$species" ]
    then
	echo "Required directories and files were not found."
	exit
    fi

    mkdir -p $out/$d/other
    if [ ! -f "$out/$d/other/mash_samples.$species.txt" ]
    then
	echo -e "date\tsample\tprokka_file\tserotyping_file_type\tserotyping_file" > $out/$d/other/mash_samples.$species.txt
    fi

    while read other
    do
	other_serotype=($(echo $other | awk '{ print $5 }' | xargs head -n 1 | grep "$species" | cut -f 1 | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "\t" $2 }'))
	if [ -n "${other_serotype[0]}" ]
	then
	    mkdir -p $out/$d/other/${other_serotype[0]}
	    mkdir -p $out/$d/other/${other_serotype[0]}/all
	    mkdir -p $out/$d/other/${other_serotype[0]}/${other_serotype[1]}
	    echo -e "$other\t${other_serotype[0]}""_""${other_serotype[1]}" >> $out/$d/other/mash_samples.txt
	    other_file=$(echo $other | awk '{ print $3 }')
	    ln -s $other_file $out/$d/other/${other_serotype[0]}/all/.
	    ln -s $other_file $out/$d/other/${other_serotype[0]}/${other_serotype[1]}/.
	fi
    done < <(grep -w "mash" $input_file | awk -v pd=$pd '{ if ( $1 > pd || $1=="include" ) print $0 }' )

    date
    if [ -f "$out/$d/other/mash_samples.$species.txt" ]
    then
	echo "$out/$d/other/mash_samples.$species.txt has been created"
	echo " "
	echo "File organization is for $species with mash results is complete!"
    else
	echo "Something went wrong and files were not moved appropriately. Sorry!"
    fi
}
run_organize.seqsero ()
{
    out=$1
    d=$2
    pd=$3
    input_file=$4
    date
    echo "Organizing salmonella samples based off of serotype groups identied in seqsero"

    if [ ! -d "$out/$d" ] || [ ! -f "$input_file" ]
    then
	echo "Required directories and files were not found."
	exit
    fi
    mkdir -p $out/$d/salmonella
    mkdir -p $out/$d/salmonella/salmonella
    mkdir -p $out/$d/salmonella/salmonella/all

    if [ ! -f "$out/$d/salmonella/salmonella_samples.txt" ]
    then
	echo -e "date\tsample\tprokka_file\tserotyping_file_type\tserotyping_file" > $out/$d/salmonella/salmonella_samples.txt
    fi

    while read salmonella
    do
	seqsero_serotype=$(echo $salmonella | awk '{ print $5 }' | xargs cat | grep "Predicted serotype(s):" | awk '{ $1=$2="" ; print $0 }' | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/O5-//g')
	echo -e "$salmonella\t$seqsero_serotype" >> $out/$d/salmonella/salmonella_samples.txt
	mkdir -p $out/$d/salmonella/salmonella/$seqsero_serotype
	salmonella_contig=$(echo $salmonella | awk '{ print $3 }')
	ln -s $salmonella_contig $out/$d/salmonella/salmonella/$seqsero_serotype/.
	ln -s $salmonella_contig $out/$d/salmonella/salmonella/all/.
    done < <( grep -w "seqsero" $input_file | awk -v pd=$pd '{ if ( $1 > pd || $1=="include" ) print $0 }' )

    date
    if [ -f "$out/$d/salmonella/salmonella_samples.txt" ]
    then
	echo "$out/$d/salmonella/salmonella_samples.txt has been created"
	echo " "
	echo "File organization for samples with seqsero results is complete!"
    else
	echo "Something went wrong and file were not moved appropriately. Sorry!"
    fi
}
run_raxml ()
{
    out=$1
    d=$2

    random_number1=$(echo "$RANDOM % 10000 + 1" | bc )
    random_number2=$(echo "$RANDOM % 10000 + 1" | bc )

    date
    if [ -z "$(find $out/$d/*/*/*/Roary_out -iname core_gene_alignment.aln -type f | head -n 1 )" ]
    then
	echo "Cannot find required files"
	exit
    fi

    echo "Running raxmlHPC-SSE3 on core alignment files:"
    echo "raxmlHPC-SSE3 -s directory/Roary_out/core_gene_alignment.aln -d -m GTRGAMMA -w directory/RAXML -n raxml -p $random_number1 -N 100 -f a -x $random_number2 -o control "

    raxml_array=()

    ls -d $out/$d/*/*/*/Roary_out | sed 's/Roary_out//g' | parallel " mkdir -p {}/RAXML "
    roary_directories=($(ls -d $out/$d/*/*/*))
    for directory in ${roary_directories[@]}
    do
	if [ -f "$directory/Roary_out/core_gene_alignment.aln" ]
	then
	    if [ -n "$(find $directory -iname '*_control.gff' | head -n 1 )" ]
	    then
		outgroup=$(ls $directory/*_control.gff | head -n 1 )
		outgroup_basename=$(basename "$outgroup" .gff)
		raxml_array+=("$(echo "raxmlHPC-SSE3___-s___$directory/Roary_out/core_gene_alignment.aln___-d___-m___GTRGAMMA___-w___$directory/RAXML___-n___raxml___-p___$random_number1""___-N___100___-f___a___-x___$random_number2""___-o___$outgroup_basename " )")
	    else
		raxml_array+=("$(echo "raxmlHPC-SSE3___-s___$directory/Roary_out/core_gene_alignment.aln___-d___-m___GTRGAMMA___-w___$directory/RAXML___-n___raxml___-p___$random_number1""___-N___100___-f___a___-x___$random_number2" )")
	    fi
	fi
    done

    history -p ${raxml_array[@]} | sort | uniq | sed 's/___/ /g' | parallel --jobs 10 '{}'

    date
    echo "RAXML has created these raxmlHPC-SSE3 files:"
    ls $out/$d/*/*/*/RAXML/*
    echo "RAXML is complete!"
}
run_recent ()
{
    out=$1
    d=$2
    recent_date=$(date -d "$d - 10 days" +%Y-%m-%d)

    date
    echo "Identifying recent samples"

    if [ -f "$out/$d/ecoli/ecoli_samples.txt" ]
    then
	while read ecoli
	do
	    sample=$(echo "$ecoli" | awk '{print $2}' )
	    if [ -d "$out/$d/ecoli/all/all" ]
	    then
		echo "$sample" >> $out/$d/ecoli/all/all/recent_samples.txt
	    fi
	    O_group=($(echo "$ecoli" | awk '{print $6}' | sed 's/O://g' | sed 's/H://g' | awk -F ";" '{ print $1}' ))
	    H_group=($(echo "$ecoli" | awk '{print $6}' | sed 's/O://g' | sed 's/H://g' | awk -F ";" '{ print $2}' ))

	    for o in ${O_group[@]}
	    do
		if [ -d "$out/$d/ecoli/$o/all" ]
		then
		    echo "$sample" >> $out/$d/ecoli/$o/all/recent_samples.txt
		fi
		for h in ${H_group[@]}
		do
		    if [ -d "$out/$d/ecoli/$o/$h" ]
		    then
			echo "$sample" >> $out/$d/ecoli/$o/$h/recent_samples.txt
		    fi
		done
	    done
	done < <( grep "abricate" $out/$d/ecoli/ecoli_samples.txt | awk -v pd=$recent_date '{ if ( $1 > pd ) print $0 }' )
    fi

    if [ -f "$out/$d/other/mash_samples.txt" ]
    then
	while read other
	do
	    sample=$(echo "$other" | awk '{print $2}')
	    other_serotype=($(echo "$other" | awk '{print $6}' | awk -F "_" '{ print $1 "\t" $2 }' ))
	    if [ -d "$out/$d/other/${other_serotype[0]}/${other_serotype[1]}" ]
	    then
		echo "$sample" >> $out/$d/other/${other_serotype[0]}/${other_serotype[1]}/recent_samples.txt
	    fi
	    if [ -d "$out/$d/other/${other_serotype[0]}/all" ]
	    then
		echo "$sample" >> $out/$d/other/${other_serotype[0]}/all/recent_samples.txt
	    fi
	done < <( grep "mash" $out/$d/other/mash_samples.txt | awk -v pd=$recent_date '{ if ( $1 > pd ) print $0 }' )
    fi

    if [ -f "$out/$d/salmonella/salmonella_samples.txt" ]
    then
	while read salmonella
	do
	    sample=$(echo "$salmonella" | awk '{print $2}')
	    serotype=$(echo "$salmonella" | awk '{print $6}')
	    if [ -d "$out/$d/salmonella/salmonella/all" ]
	    then
		echo "$sample" >> $out/$d/salmonella/salmonella/all/recent_samples.txt
	    fi
	    if [ -d "$out/$d/salmonella/salmonella/$serotype" ]
	    then
		echo "$sample" >> $out/$d/salmonella/salmonella/$serotype/recent_samples.txt
	    fi
	done < <( grep "seqsero" $out/$d/salmonella/salmonella_samples.txt | awk -v pd=$recent_date '{ if ( $1 > pd ) print $0 }' )
    fi

    cat $out/$d/*/*/*/recent_samples.txt > $out/$d/recent_samples.txt
    echo "Identification of recent samples is complete"
}
run_roary ()
{
    out=$1
    d=$2
    max_samples=$3
    kraken_path=$4
    recent_flag=$5

    date
    if [ ! -d "$out/$d" ] || [ -z "$(find $out/$d/*/*/* -iname '*gff' | head -n 1 )" ]
    then
	echo "Required directories and files were not found."
	exit
    fi

    echo "running Roary:"
    echo "roary -p 48 -f directory/Roary_out -e -n -qc -k $kraken_path directory/*.gff"

    roary_array=()
    directories=""
    if [ -z "$recent_flag" ]
    then
	directories=($(ls -d $out/$d/*/*/*/recent_samples.txt | sed 's/\/recent_samples.txt//g' ))
    else
	directories=($(ls -d $out/$d/*/*/*))
    fi

    for directory in ${directories[@]}
    do
	if [ -n "$(find $directory -iname '*gff' | head -n 1 )" ]
	then
	    number_of_gff=$(ls $directory/*gff | wc -l )
	    if (( $number_of_gff > 4 )) && (( $number_of_gff < $max_samples ))
	    then
		roary_array+="$(echo " roary___-p___48___-f___$directory/Roary_out___-e___-n___-qc___-k___""$kraken_path""___$directory/*.gff " )"
	    elif (( $number_of_gff > 4 )) && (( $number_of_gff >= $max_samples ))
	    then
		echo "$directory will be skipped because it has $number_of_gff gff files. This will take too long for iqtree to run. Sorry"
		ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | parallel 'echo -e "{}\ttoo_many_samples" >> $out/$d/$d.samples.rm '
	    else
		echo "$directory does not have enough gff files for roary and iqtree"
		ls $directory/*gff | sed 's!.*/!!' | cut -d "." -f 1 | parallel 'echo -e "{}\ttoo_few_samples" >> $out/$d/$d.samples.rm '
	    fi
	else
	    echo "$directory does not have *gff files"
	fi
    done

    date
    history -p ${roary_array[@]} | sort | uniq | sed 's/___/ /g' | parallel --jobs 10 '{}'

    date
    echo "Sample,Genus,Species" > $out/$d/Roary_qc_report.csv
    cat $out/$d/*/*/*/Roary_out/qc_report.csv | grep -v "Sample,Genus,Species" | sort | uniq >> $out/$d/Roary_qc_results.csv
    organisms=($(cat $out/$d/*/*/*/Roary_out/qc_report.csv | grep -v "Sample,Genus,Species" | cut -f 3 -d "," | sort | uniq | tr ' ' '_' ))
    header="Sample"
    for organism in ${organisms[@]}
    do
	header=$(echo "$header,$organism" )
    done
    echo -e "$header" > $out/$d/Roary_qc_report.csv

    date
    echo -e "Set\tCore_genes\tSoft_core_genes\tShell_genes\tCloud_genes" > $out/$d/Roary_summary_statistics.txt
    for directory in ${directories[@]}
    do
	split_paths=($(echo "$directory" | awk -F "/" '{print $(NF-3) "\t" $(NF-2) "\t" $(NF-1) "\t" $(NF) }' ))
	#ls $directory/Roary_out/gene_presence_absence.Rtab | parallel "cp {} {}.${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}.summary_mqc.txt"

	core_genes=$(grep "Core genes"      $directory/Roary_out/summary_statistics.txt | cut -f 3 )
	soft_genes=$(grep "Soft core genes" $directory/Roary_out/summary_statistics.txt | cut -f 3 )
	shel_genes=$(grep "Shell genes"     $directory/Roary_out/summary_statistics.txt | cut -f 3 )
	clod_genes=$(grep "Cloud genes"     $directory/Roary_out/summary_statistics.txt | cut -f 3 )
	totl_genes=$(grep "Total genes"     $directory/Roary_out/summary_statistics.txt | cut -f 3 )
	echo -e "${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}\t$core_genes\t$soft_genes\t$shel_genes\t$clod_genes" >> $out/$d/Roary_summary_statistics.txt

	organism_counts="${split_paths[0]}.${split_paths[1]}.${split_paths[2]}.${split_paths[3]}"
	for organism in ${organisms[@]}
	do
	    number=$(cat $directory/*/*/*/Roary_out/qc_report.csv | grep -v "Sample,Genus,Species" | cut -f 3 -d "," | tr ' ' '_' | sort | uniq -c | grep "$organism" | awk '{ print $1 }' )
	    if [ -z "$number" ]
	    then
		number="0"
	    fi
	    organism_counts="$organism_counts,$number"
	done
	echo $organism_counts >> $out/$d/Roary_qc_report.csv

    done

    date
    if [ -n "$(find $out/$d/*/*/*/Roary_out -iname core_gene_alignment.aln | head -n 1 )" ]
    then
	echo "core gene alignment files are located in"
	ls $out/$d/*/*/*/Roary_out/core_gene_alignment.aln
	echo "Roary has been successful"
    else
	echo "Something went wrong with Roary. Sorry!"
    fi
}
run_search ()
{
    out=$1
    d=$2
    pd=$3
    search_path=$4

    date
    echo "Searching through $search_path for prokka gff files"

    if [ -z "$(find $search_path/*/ALL_gff -iname '*gff' | head -n 1 )" ]
    then
	echo "No gff files were found."
	exit
    fi

    gff_files=($(find $search_path/*/ALL_gff -iname *gff -newermt $pd ! -newermt $d -type f ))

    date
    echo "Searching through $search_path for seqsero, abricate, and mash results for samples"

    if [ ! -f "$out/$d/input_file.txt" ]
    then
	echo -e "date\tsample\tprokka_file\tserotyping_file_type\tserotyping_file" > $out/$d/input_file.txt
    fi

    for prokka in ${gff_files[@]}
    do
	modification_date=$(stat -c %y "$prokka" | awk '{ print $1}' )
	sample_name=$(echo $prokka | sed 's!.*/!!' | cut -d "." -f 1 )

	seqsero_result=()
	abricate_result=()
	mash_result=()

	if [ -n "$(find $search_path/*/mash -iname $sample_name*.msh.distance.txt.sorted.txt | head -n 1 )" ]
	then
	    mash_result=($(ls $search_path/*/mash/$sample_name*.msh.distance.txt.sorted.txt))
	    echo -e "$modification_date\t$sample_name\t$prokka\tmash\t${mash_result[0]}" >> $out/$d/input_file.txt
	else
	    echo "Could not find mash result for $sample_name"
	fi

	if [ -n "$(find $search_path/*/SeqSero -iname $sample_name*.Seqsero_result.txt | head -n 1 )" ]
	then
	    seqsero_result=($(ls $search_path/*/SeqSero/$sample_name*.Seqsero_result.txt ))
	    echo -e "$modification_date\t$sample_name\t$prokka\tseqsero\t${seqsero_result[0]}" >> $out/$d/input_file.txt
	else
	    echo "Could not find seqsero result for $sample_name" | grep "PNUSAS"
	fi

	if [ -n "$(grep $sample_name $search_path/*/abricate_results/*serotypefinder*out.txt | head -n 1 )" ]
	then
	    abricate_result=($(grep $sample_name $search_path/*/abricate_results/*serotypefinder*out.txt | cut -f 1 -d ":" ))
	    echo -e "$modification_date\t$sample_name\t$prokka\tabricate\t${abricate_result[0]}" >> $out/$d/input_file.txt
	else
	    echo "Could not find abricate result for $sample_name " | grep "PNUSAE"
	fi

	if [ -z "${mash_result[0]}" ] && [ -z "${seqsero_result[0]}" ] && [ -z "$abricate_result[0]}" ]
	then
	    echo -e "$sample_name\tno_serotype" >> $out/$d/$d.samples.rm
	    echo "WARN:Could not find serotyping file for $prokka !"
	    echo -e "$modification_date\t$sample_name\t$prokka\tunknown\tunknown" >> $out/$d/input_file.txt
	fi
    done

    date
    echo "Combining Abricate results"
    ls $search_path/*/abricate_results/*out.txt | awk -F "/" '{ print $NF }' | sed 's/.out.txt//g' | sort | uniq | parallel "cat $search_path/*/abricate_results/{}*out.txt > $out/$d/{}_abricate_results.txt "

    date
    if [ -f "$out/$d/input_file.txt" ]
    then
	echo "The input_file has been created and is located at $out/$d/input_file.txt"
	echo "Input file creation is complete!"
    else
	echo "Something went wrong and no input file was created. Sorry!"
    fi
}
#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
while getopts 'ac:d:e:f:hi:k:m:o:p:qs:t:vxAC:EGIOP:QRVX' OPTION
do
    case "$OPTION" in
	a)
	    echo "Running roary, iqtree, and ggtree on all directories, not just those with samples less than 10 days old"
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
	k)
	    kraken_path=$OPTARG
	    if [ ! -f "$kraken_path" ]
	    then
		echo "Path for kraken could not be found"
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
	    echo "Will look through mash results for samples that are $species"
	    ;;
	q)
	    both_flag=1
	    echo "Will run both IQTREE and raxmlHPC-SSE3"
	    ;;
	s)
	    echo "Ending date is now $OPTARG instead of $d"
	    d=$OPTARG
	    ;;
	t)
	    if [[ "$OPTARG" == "mash" || "$OPTARG" == "abricate" || "$OPTARG" == "seqsero" ]]
	    then
		echo "Results are restricted to those found in $OPTARG"
		serotyping=("${serotyping[@]}" "$OPTARG")
	    else
		echo "Invalid option: $OPTARG can only accept mash/seqsero/abricate"
		echo "$USAGE"
		exit 1
	    fi
	    ;;
	v)
	    date
	    echo "OUTBREAK_120.sh v0.2"
	    roary_version="$(roary -w)"
	    echo "Roary $roary_version"
	    iqtree -version
	    raxml_version=$(raxmlHPC-SSE3 -v | grep "version" )
	    echo "raxmlHPC-SSE3 $raxml_version"
	    # ape
	    # ggtree
	    exit 1
	    ;;
	x)
	    echo "Running raxml instead of iqtree"
	    raxml_flag=1
	    ;;
	A)
	    flag=1
	    echo "Not using full pipeline. Organizing abricate results for resistence information"
	    if [ -f "$out/$d/all_abricate_results.txt" ]
	    then
		abricate_flag=1
	    else
		echo "Cannot find abricate results at $out/$d/all_abricate_results.txt"
	    fi
	    ;;
	C)
	    control_file=$OPTARG
	    if [ ! -f "$control_file" ]
	    then
		echo "Cannot find file listing controls and serotypes"
		echo $USAGE
		exit
	    fi
	    ;;
	E)
	    flag=1
	    echo "Not using full pipeline. Will clean up files and add controls"
	    if [ -d "$out/$d" ]
	    then
		clean_flag=1
	    else
		echo "Cannot find directories"
		echo $USAGE
		exit
	    fi
	    ;;
	G)
	    flag=1
	    echo "Not using full pipeline. Running GGTREE"
	    if [ -d "$out/$d" ]
	    then
		ggtree_flag=1
	    else
		"Cannot find directories with iqtree and roary results"
		exit
	    fi
	    ;;
	I)
	    flag=1
	    echo "Not using full pipeline. Creating input file."
	    if [ -d "$search_path" ] && [ -n "$out" ]
	    then
		search_flag=1
	    else
		echo "Cannot create input file from $search_path"
	    fi
	    ;;
	O)
	    flag=1
	    echo "Not using full pipeline. Organizing files from input file(s)"
	    organize_flag=1
	    ;;
	P)
	    echo "Adding $OPTARG to PATH"
	    PATH=$PATH:$OPTARG
	    echo "PATH is now $PATH"
	    ;;
	Q)
	    flag=1
	    echo "Not using full pipeline. Running IQTREE on core genomes in subfolders in $out/$d"
	    if [ -d "$out/$d" ]
	    then
		iqtree_flag=1
		iqtree_flag_force=1
	    else
		echo "Cannot find directories with core alignment files for IQtree"
		echo $USAGE
		exit 1
	    fi
	    ;;
	R)
	    flag=1
	    echo "Not using full pipeline. Running Roary on subfolders under $out/$d"
	    if [ -d "$out/$d" ]
	    then
		roary_flag=1
	    else
		echo "Cannot find directory for Roary"
		echo $USAGE
		exit 1
	    fi
	    ;;
	V)
	    flag=1
	    echo "Checking dependencies"
	    run_dependencies
	    exit 1
	    ;;
	X)
	    flag=1
	    echo "Not using full pipeline. Running raxmlHPC-SSE3 on core genomes in subfolders $out/$d"
	    raxml_flag=1
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

if [ -z "$flag" ]
then
    search_flag=1
    organize_flag=1
    clean_flag=1
    control_flag=1
    roary_flag=1
    iqtree_flag=1
    abricate_flag=1
    ggtree_flag=1
fi

if [ -n "$raxml_flag" ]
then
    iqtree_flag=""
fi

if [ -n "$both_flag" ]
then
    iqtree_flag=1
    raxml_flag=1
fi

if [ -n "$iqtree_flag_force" ]
then
    iqtree_flag=1
fi

run_dependencies

if [ ! -f "${input_files[0]}" ]
then
    echo "WARNING: no input file was specified!"
fi

if [ -z "$out" ]
then
    echo "Output directory is required"
    echo "$USAGE"
    exit
fi

mkdir -p $out
mkdir -p $out/$d
if [ ! -d "$out/$d" ]
then
    echo "Could not create directory for results!"
    exit
fi

mkdir -p $out/$d/logs
pd=$(date -d "$d - $age days" +%Y-%m-%d)

if [ -n "$search_flag" ]
then
    if [ ! -d "$search_path" ]
    then
	echo "$search_path is not a directory that can be searched through!"
	echo "$USAGE"
	exit
    fi

    date
    echo "Searching in $search_path for Prokka and serotyping files"
    run_search $out $d $pd $search_path > $out/$d/logs/$run_date.search.log 2> $out/$d/logs/$run_date.search.err
    input_files=("${input_files[@]}" "$out/$d/input_file.txt")
fi

if [ -n "$organize_flag" ]
then
    if [ ! -f "${input_files[0]}" ]
    then
	echo "Could not find input files:  ${input_files[@]}"
	echo "$USAGE"
	exit
    fi
    for input in ${input_files[@]}
    do
	echo "Organizing samples from $input"
	if [ -n "${serotyping[0]}" ] && [ -z "$species" ] && [ -f "$input" ]
	then
	    for kind in ${serotyping[@]}
	    do
		date
		echo "Restricting results to those with $kind results"
		run_organize.$kind $out $d $pd $input > $out/$d/logs/$run_date.organize.$kind.log 2>$out/$d/logs/$run_date.organize.$kind.err
	    done
	elif [ -n "$species" ] && [ -f "$input" ]
	then
	    date
	    echo "Restricting results to $species with mash results"
	    run_organize.mash.subtype $out $d $pd $input $species > $out/$d/logs/$run_date.organize.$species.log 2>$out/$d/logs/$run_date.organize.$species.err
	else
	    for kind in mash abricate seqsero
	    do
		date
		echo "Starting organization of samples with $kind results"
		run_organize.$kind $out $d $pd $input > $out/$d/logs/$run_date.organize.$kind.log 2>$out/$d/logs/$run_date.organize.$kind.err
	    done
	fi
    done

    date
    echo "Specifying folders containing samples less than 10 days old"
    run_recent $out $d > $out/$d/logs/$run_date.recent.log 2>$out/$d/logs/$run_date.recent.err
fi

if [ -n "$clean_flag" ]
then
    if [ ! -f "$out/$d/$d.samples.rm" ]
    then
    	echo -e "sample\treason" > $out/$d/$run_date.samples.rm
    fi

    date
    echo "Removing samples with files that are too small"
    run_clean.size $out $d > $out/$d/logs/$run_date.clean.log 2> $out/$d/logs/$run_date.clean.err

    if [ -f "$exclude_file" ]
    then
	date
	echo "Removing samples specied in $exclude_file"
	run_clean.exclude $out $d $exclude_file >> $out/$d/logs/$run_date.clean.log 2>> $out/$d/logs/$run_date.clean.err
    fi

    date
    echo "Removing samples with abnormal serotypes (like missing O groups in E. coli)"
    run_clean.prune $out $d >> $out/$d/logs/$run_date.clean.log 2>> $out/$d/logs/$run_date.clean.err

    date
    echo "Adding controls to samples"
    run_control $out $d $control_path >> $out/$d/logs/$run_date.clean.log 2>> $out/$d/logs/$run_date.clean.err
fi

if [ -n "$roary_flag" ]
then
    date
    echo "Running roary"
    run_roary $out $d $max_samples $kraken_path $recent_flag > $out/$d/logs/$run_date.roary.log 2> $out/$d/logs/$run_date.roary.err
fi

if [ -n "$iqtree_flag" ]
then
    date
    echo "Running IQTREE"
    run_iqtree $out $d > $out/$d/logs/$run_date.iqtree.log 2> $out/$d/logs/$run_date.iqtree.err
fi

if [ -n "$raxml_flag" ]
then
    date
    echo "Running raxmlHPC-SSE3"
    run_raxml $out $d > $out/$d/$run_date.raxml.log 2> $out/$d/$run_date.raxml.err
fi

if [ -n "$abricate_flag" ]
then
    date
    echo "Creating Abricate resfinder presence/absence tables"
    run_abricate $out $d > $out/$d/logs/$run_date.resistence.log 2> $out/$d/logs/$run_date.resistence.err
fi

if [ -n "$ggtree_flag" ]
then
    date
    echo "Running GGTREE"
    run_ggtree $out $d > $out/$d/logs/$run_date.ggtree.log 2> $out/$d/logs/$run_date.ggtree.err
fi

date
#multiqc $out/$d -d --outdir $out/$d/logs

date
echo "log files can be found at $out/$d/logs"
echo "OUTBREAK_120 has completed!"


# TO DO:
# add in SNPs
# get better controls
