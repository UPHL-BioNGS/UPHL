#!/bin/bash

#----#----#----#----#----#----#----#----#----#----

USAGE="
This BASH script is meant to take fastq files from
bacterial sequencing to an annotated genome.

USAGE: ./WGS_NGS_pipeline.sh -i /path/to/fastq -o /out/directory

A full list of options is listed with 
./UPHL_reference_free.sh -h
"

#----#----#----#----#----#----#----#----#----#----

HELP="
WGS_NGS_pipeline.sh moves the fastq files to the
output directory and runs multiple algorithms.

For QC:
1. FASTQC (unpublished)
2. QUAST (doi: 10.1093/bioinformatics/btt086)
3. SEQYCLEAN (doi: 10.1145/3107411.3107446)
4. MULTIQC (doi: 10.1093/bioinformatics/btw354)

For SEROTYPING or RESISTENCE GENE DETECTION
1. MASH (doi: 10.1186/s13059-016-0997-x)
2. SEQSERO (doi: 10.1128/JCM.00323-15)
3. ABRICATE (unpublished)

For ALIGNMENT and GENOME ANNOTATION
1. SHOVILL (unpublished)
2. PROKKA (doi: 10.1093/bioinformatics/btu153)

REQUIRED OPTIONS:

-i /path/to/fastq/files          Full path to fastq files. All fastq files must be in the same folder.
-o /output/directory             The full path of the directory where all files will be moved or created.
-r <BaseSpace Sequencing Run>    If BaseSpace has been mounted at /home/BaseSpace, the run identifier can be used instead of using \"-o\" or \"-i\"

OPTIONAL FLAGS:

-C /path/to/seqyclean/adaptors   Path to SEQYCLEAN adaptors
-M /path/to/refseq/mash/sketches Path to RefSeqSketchesDefaults.msh
-v                               Display versions and paths of dependencies and exit
-b                               Display script commands and exit
-P <PATH>                        Add path(s) to PATH variable

PARTIAL PIPELINE FLAGS (Does not run full pipeline):
-c                               Run SEQYCLEAN
-l                               Run SHOVILL  (Requires SEQYCLEAN files)
-p                               Run PROKKA   (Requires SHOVILL contigs)
-f                               Run FASTQC   (Requires fastq files, but will also work on SEQYCLEAN files)
-q                               Run Quast    (Requires SHOVILL contigs) 
-m                               Run Mash     (Requires SEQYCLEAN files)
-g                               Locate Salmonella and E. coli samples from MASH results 
-s                               Run SEQSERO  (Requires SEQYCLEAN files); NOTE: Salmonella samples can be listed in out/directory/salmonella.txt
-a                               Run ABRICATE (Requires SHOVILL contigs); NOTE: E. coli samples can be listed in out/directory/ecoli.txt
-x                               Check if all files have been generated

examples:
./UPHL_reference_free.sh -r UT-M78932-170504 -l -p
./UPHL_reference_free.sh -o /home/NGS -m -g -s -a

"
run_date=$(date +%Y%m%d)
seqyclean_adaptors=/home/Bioinformatics/Data/SeqyClean_data/PhiX_174_plus_Adapters.fasta
mash_sketches=/home/Bioinformatics/Data/RefSeqSketchesDefaults.msh

#----#----#----#----#--FUNCTIONS--#----#----#----#----#----

run_abricate ()
{
    out=$1
    date
    if [ -z "$out" ] || [ -z "$(find $out/ALL_assembled -name '*_contigs.fa' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi
    
    echo "Now running ABRicate on the samples for resistance, virulence, and E. Coli serotyping"
    echo "COMMAND: abricate -db resfinder      out/ALL_assembled/*_contigs.fa > out/abricate_results/resfinder_out.txt "
    echo "COMMAND: abricate -db vfdb           out/ALL_assembled/*_contigs.fa > out/abricate_results/vfdb_out.txt "
    echo "COMMAND: abricate -db serotypefinder out/ALL_assembled/*_contigs.fa > out/abricate_results/Ecoli_serotypefinder_out.txt "
    echo "COMMAND: abricate -db ecoh           out/ALL_assembled/*_contigs.fa > out/abricate_results/ecoh_out.txt"

    mkdir -p $out/abricate_results

    abricate -db resfinder $out/ALL_assembled/*_contigs.fa > $out/abricate_results/resfinder_out.txt
    abricate -db vfdb      $out/ALL_assembled/*_contigs.fa > $out/abricate_results/vfdb_out.txt
    
    ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "echo -e '#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT' | tee $out/abricate_results/resfinder.{}.tab $out/abricate_results/vfdb.{}.tab " 
    ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "grep {} $out/abricate_results/resfinder_out.txt >> $out/abricate_results/resfinder.{}.tab"
    ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "grep {} $out/abricate_results/vfdb_out.txt >> $out/abricate_results/vfdb.{}.tab"
    abricate --summary $out/abricate_results/resfinder*tab > $out/abricate_results/resfinder_summary.tab
    abricate --summary $out/abricate_results/vfdb*tab      > $out/abricate_results/vfdb_summary.tab
    cat $out/abricate_results/resfinder_summary.tab | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }'> $out/abricate_results/resfinder_summary.Rtab
    cat $out/abricate_results/vfdb_summary.tab      | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }' > $out/abricate_results/vfdb_summary.Rtab

    if  [ -f "$out/ecoli.txt" ] || [ -n "$(find $out/Sequencing_reads/Raw/ -name 'PNUSAE*' | head -n 1 )" ]
    then
	abricate -db serotypefinder $out/ALL_assembled/*_contigs.fa > $out/abricate_results/Ecoli_serotypefinder_out.txt
	abricate -db ecoh           $out/ALL_assembled/*_contigs.fa > $out/abricate_results/ecoh_out.txt
	ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "echo -e '#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT' | tee $out/abricate_results/serotypefinder.{}.tab $out/abricate_results/ecoh.{}.tab " 
	ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "grep {} $out/abricate_results/ecoh_out.txt >> $out/abricate_results/ecoh.{}.tab"
	ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 "grep {} $out/abricate_results/Ecoli_serotypefinder_out.txt >> $out/abricate_results/serotypefinder.{}.tab"
	abricate --summary $out/abricate_results/ecoh*tab > $out/abricate_results/ecoh_summary.tab
	abricate --summary $out/abricate_results/serotypefinder*tab  > $out/abricate_results/serotypefinder_summary.tab
	cat $out/abricate_results/ecoh_summary.tab           | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }' > $out/abricate_results/ecoh_summary.Rtab
	cat $out/abricate_results/serotypefinder_summary.tab | sed 's/#//g' | sed 's/.tab//g' | awk '{ sub("^.*/", "", $1); print}' | awk '{ for (i=1;i<=NF;i++) if ($i ~ ";" )gsub(";.*$","",$i)g ; else continue}{print $0}' | awk '{ $2="" ; print $0 }' > $out/abricate_results/serotypefinder_summary.Rtab

    else
	echo "No ecoli samples were indicated, so Ecoli serotyping was skipped."
    fi
    
    date
    if [ -n "$(find $out/abricate_results/ -iname '*out.txt' | head -n 1 )" ]
    then 
	echo "ABRicate is finished."
	echo "The results are located in $out/abricate_results"
    else
	echo "Something went wrong with Abricate. There are no results!"
    fi
}
run_basespace ()
{
    run=$1
    out=$2

    if [ -z "$out" ] || [ -z "$run" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    mkdir -p $out/Sequencing_reads
    mkdir -p $out/Sequencing_reads/Raw

    date
    echo "Moving FASTQ files from BaseSpace to the desired output file."
    echo "The fastq files from run $run were located in the following locations on BaseSpace:"
    
    while read line
    do
	SAMPLE=($(echo $line | awk -F "," '{print $1 "\t" $2 "\t" $9 }' ))
	if [ -n "$(find /home/BaseSpace/Projects/"${SAMPLE[2]}"/Samples/${SAMPLE[0]}/Files/ -iname ${SAMPLE[1]}*fastq* | head -n 1 )" ]
	then
	    ls /home/BaseSpace/Projects/"${SAMPLE[2]}"/Samples/${SAMPLE[0]}/Files/${SAMPLE[1]}*fastq* | tee -a $out/Logs/BaseSpace_locations.txt
	    cp /home/BaseSpace/Projects/"${SAMPLE[2]}"/Samples/${SAMPLE[0]}/Files/${SAMPLE[1]}*fastq* $out/Sequencing_reads/Raw/.
	else
	    echo "Could not find the fastq files in /home/BaseSpace/Projects/${SAMPLE[2]}/Files/"
	    echo "Expanding search for sample ${SAMPLE[0]} to other projects."
	    ls /home/BaseSpace/Projects/*/Samples/${SAMPLE[0]}*/Files/${SAMPLE[1]}*fastq* | tee -a $out/Logs/BaseSpace_locations.txt
	    cp /home/BaseSpace/Projects/*/Samples/${SAMPLE[0]}*/Files/${SAMPLE[1]}*fastq* $out/Sequencing_reads/Raw/.
	fi
    done < <(grep "," /home/BaseSpace/Runs/$run/Files/SampleSheet.csv | grep -v "Sample_ID" | grep -v "Chemistry" | grep -v "Workflow" | grep -v "Module" | grep -v "Date" | grep -v "Experiment Name" | grep -v "adapter" | grep -v "Library Prep" | grep -v "Local Run Manager" | grep -v "IEMFileVersion" | grep -v "Investigator" | grep -v "Application" | grep -v "Assay" | grep -v "Description" | grep -v "ReverseComplement" | grep -v "Adapter" )

    date
    if [ -n "$(find $out/Sequencing_reads/Raw/ -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The BaseSpace FASTQ files from $run have now been moved to"
	ls $out/Sequencing_reads/Raw/*
	echo "Moving files from BaseSpace is complete!"
    else
	echo "Something went wrong! No files were moved!!!"
    fi
}
run_complete ()
{
    out=$1
    date

    if [ -z "$(find $out/Sequencing_reads/Raw/ -iname '*fastq*' | head -n 1 )" ]
    then
	echo "No fastq files were found in $out/Sequencing_reads/Raw"
	exit
    fi

    echo -e "sample_id\trun_id\tout_directory\tR1_fastq\tR2_fastq\tClean_R1_fastq\tClean_R2_fastq\tfastqc_R1_raw\tfastqc_R1_raw_reads\tfastqc_R2_raw\tfastqc_R2_raw_reads\tfastqc_R1_clean\tfastqc_R1_clean_reads\tfastqc_R2_clean\tfastqc_R2_clean_reads\tshovill_file\tprokka_file\tquast_file\tmash_file\tmash_result\tsimple_mash_result\tseqsero_file\tseqsero_serotype\tseqsero_antigenic_profile\tsimple_seqsero_result\tecoh_file\tecoh_result_O\tecoh_result_H\tabricate_serotypefinder_file\tabricate_serotypefinder_result_O\tabricate_serotypefinder_H\tresfinder_file\tresfinder_result\tvfdb_file\tvvfdb_result" > $out/Run_summary.txt
    samples=($(ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq ))
    for pnusa in ${samples[@]}
    do
	sample=$(echo $pnusa | sed 's/-UT.*//g' )
	# the raw fastq files
	if [ -n "$(find $out/Sequencing_reads/Raw -iname $pnusa* )" ]
	then
	    raw_reads=($(ls $out/Sequencing_reads/Raw/$pnusa* )) 
	else
	    raw_reads=("not_found" "not_found" )
	    echo "WARNING: Raw reads were not found for $sample"
	fi
	
	# cleaned fastq files
	if [ -n "$(find $out/Sequencing_reads/QCed -iname $pnusa* )" ]
	then
	    clean_reads=($(ls $out/Sequencing_reads/QCed/$pnusa* ))
	else
	    clean_reads=("not_found" "not_found" )
	    echo "WARNING: Clean reads were not found for $sample"
	fi
	
	# mash results
	if [ -n "$(find $out/mash -iname $pnusa* )" ]
	then
	    mash_file=$(ls $out/mash/$pnusa*msh.distance.txt.sorted.txt | head -n 1 )
	    mash_results=($(head -n 1 $mash_file | cut -f 1 | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' ))
	    mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}""_""${mash_results[2]}""_""${mash_results[3]}")
	    simple_mash_result=$(echo "${mash_results[0]}""_""${mash_results[1]}")
	    if [ -z "$mash_result" ]
	    then
		mash_result="no_result"
		simple_mash_result="no_result"
		echo "WARNING: MASH results were not found for $sample"
	    fi
	    ecoli_flag=""
	    salm_flag=""
	    if [ -n "$(echo $mash_result | grep "Escherichia_coli" )" ]
	    then
		ecoli_flag=1
	    elif [ -n "$(echo $mash_result | grep "Salmonella_enterica" )" ]
	    then
		salm_flag=1
	    fi
	else
	    echo "WARNING: MASH results were not found for $sample"
	    mash_file="not_found"
	    mash_result="no_result"
	    simple_mash_result="no_result"
	fi
	
	# abricate results
	# ecoh
	if [ -f "$out/abricate_results/ecoh_out.txt" ] && [ -n "$ecoli_flag" ]
	then
	    ecoh_file=$out/abricate_results/ecoh_out.txt
	    ecoh_results=$(grep $pnusa $out/abricate_results/ecoh_out.txt | head -n 1 ) 
	    if [ -n "$ecoh_results" ]
	    then
		ecoh_O_group=($(grep $pnusa $out/abricate_results/ecoh_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | cut -f 2 -d "-" | sort | uniq | grep "O" | sed 's/\///g' ))
		ecoh_H_group=($(grep $pnusa $out/abricate_results/ecoh_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | cut -f 2 -d "-" | sort | uniq | grep "H" | sed 's/\///g' ))
		ecoh_result_O=$(echo ${ecoh_O_group[@]} | tr ' ' '_' )
		ecoh_result_H=$(echo ${ecoh_H_group[@]} | tr ' ' '_' )
		if [ -z "$ecoh_result_O" ]
		then
		    ecoh_result_O="none"
		fi
		if [ -z "$ecoh_result_H" ]
		then
		    ecoh_result_H="none"
		fi
	    else
		if [ -n "$ecoli_flag" ]
		then
		    echo "WARNING: ABRICATE ECOH results were not found for $sample"
		fi
		ecoh_result_O="not_in_file"
		ecoh_result_H="not_in_file"
	    fi
	else
	    if [ -n "$ecoli_flag" ]
	    then
		echo "WARNING: ABRICATE ECOH results were not found for $sample"
	    fi
	    ecoh_file="not_found"
	    ecoh_result_O="no_result"
	    ecoh_result_H="no_result"
	fi
	# serotype_finder
	if [ -f "$out/abricate_results/Ecoli_serotypefinder_out.txt" ] && [ -n "$ecoli_flag" ]
	then
	    abricate_serotypefinder_file=$out/abricate_results/Ecoli_serotypefinder_out.txt
	    abricate_serotypefinder_results=$(grep $pnusa $out/abricate_results/Ecoli_serotypefinder_out.txt | head -n 1 )
	    if [ -n "$abricate_serotypefinder_results" ]
	    then
		ecoli_O_group=($(grep $pnusa $out/abricate_results/Ecoli_serotypefinder_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }'| cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "O" | sed 's/\///g' ))
		ecoli_H_group=($(grep $pnusa $out/abricate_results/Ecoli_serotypefinder_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }'| cut -f 5 | cut -f 4 -d "_" | sort | uniq | grep "H" | sed 's/\///g' ))
		abricate_serotypefinder_result_O=$(echo ${ecoli_O_group[@]} | tr ' ' '_' )
		abricate_serotypefinder_result_H=$(echo ${ecoli_H_group[@]} | tr ' ' '_' )
		if [ -z "$abricate_serotypefinder_result_O" ]
		then
		    abricate_serotypefinder_result_O="none"
		fi
		if [ -z "$abricate_serotypefinder_result_H" ]
		then
		    abricate_serotypefinder_result_H="none"
		fi
	    else
		if [ -n "$ecoli_flag" ]
		then
		    echo "WARNING: ABRICATE SEROTYPE results were not found for $sample"
		fi
		abricate_serotypefinder_result_O="not_in_file"
		abricate_serotypefinder_result_H="not_in_file"
	    fi
	else
	    if [ -n "$ecoli_flag" ]
	    then
		echo "WARNING: ABRICATE SEROTYPE results were not found for $sample"
	    fi
	    abricate_serotypefinder_file="not_found"
	    abricate_serotypefinder_result_O="no_result"
	    abricate_serotypefinder_result_H="no_result"
	fi
	# resfinder
	if [ -f "$out/abricate_results/resfinder_out.txt" ]
	then
	    resfinder_file=$out/abricate_results/resfinder_out.txt
	    resfinder_results=($(grep $pnusa $out/abricate_results/resfinder_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq ))
	    resfinder_result=$(echo ${resfinder_results[@]} | tr ' ' '_' )
	    if [ -z "$resfinder_result" ]
	    then
		resfinder_result="not_in_file"
	    fi
	else
	    echo "WARNING: ABRICATE RESFINDER results were not found for $sample"
	    resfinder_file="not_found"
	    resfinder_result="not_in_file"
	fi
	# vfdb
	if [ -f "$out/abricate_results/vfdb_out.txt" ]
	then
	    vfdb_file=$out/abricate_results/vfdb_out.txt
	    vfdb_results=($(grep $pnusa $out/abricate_results/vfdb_out.txt | awk '{ if ($10 > 80) print $0 }' | awk '{ if ($9 > 80) print $0 }' | cut -f 5 | sort | uniq ))
	    vfdb_result=$(echo ${vfdb_results[@]} | tr ' ' '_' )
	    if [ -z "$vfdb_result" ]
	    then
		vfdb_result="not_in_file"
	    fi
	else
	    echo "WARNING: ABRICATE VFDB results were not found for $sample"
	    vfdb_file="not_found"
	    vfdb_result="no_result"
	fi
	
	# shovill results
	if [ -n "$(find $out/ALL_assembled -iname $pnusa* )" ]
	then
	    shovill_result=$(ls $out/ALL_assembled/$pnusa* | head -n 1 )
	else
	    echo "WARNING: SHOVILL results were not found for $sample"
	    shovill_result="not_found"
	fi
	
	# prokka results
	if [ -n "$(find $out/ALL_gff -iname $pnusa* )" ]
	then
	    prokka_result=$(ls $out/ALL_gff/$pnusa* | head -n 1 )
	else
	    echo "WARNING: PROKKA results were not found for $sample"
	    prokka_result="not_found"
	fi
	
	# fastqc files
	if [ -n "$(find $out/fastqc -iname $pnusa* )" ]
	then
	    raw_files=($(ls $out/fastqc/$pnusa*R*zip | grep -v "clean_SE_fastqc" | grep -v "clean_PE" | head -n 2 ))
	    clean_files=($(ls $out/fastqc/$pnusa*clean_PE*zip | grep -v "clean_SE_fastqc" | head -n 2 ))

	    fastqc_files[0]=${raw_files[0]}
	    fastqc_files[1]=${raw_files[1]}
	    fastqc_files[2]=${clean_files[0]}
	    fastqc_files[3]=${clean_files[1]}

	    if [ -n "${fastqc_files[0]}" ]
	    then
		fastq_summary_R1=$(unzip -l ${fastqc_files[0]} | grep fastqc_data.txt | awk '{ print $4 }' )
		fastqc_result_R1=$(unzip -p ${fastqc_files[0]} $fastq_summary_R1 | grep "Total Sequences" | awk '{ print $3 }' )
	    else
		fastqc_files[0]="not_found"
		fastqc_result_R1="no_result"
	    fi
	    if [ -n "${fastqc_files[1]}" ]
	    then
		fastq_summary_R2=$(unzip -l ${fastqc_files[1]} | grep fastqc_data.txt | awk '{ print $4 }' )
		fastqc_result_R2=$(unzip -p ${fastqc_files[1]} $fastq_summary_R2 | grep "Total Sequences" | awk '{ print $3 }' )
	    else
		fastqc_files[1]="not_found"
		fastqc_result_R2="no_result"
	    fi
	    if [ -n "${fastqc_files[2]}" ]
	    then
		fastq_summary_P1=$(unzip -l ${fastqc_files[2]} | grep fastqc_data.txt | awk '{ print $4 }' )
		fastqc_result_P1=$(unzip -p ${fastqc_files[2]} $fastq_summary_P1 | grep "Total Sequences" | awk '{ print $3 }' )
	    else
		fastqc_files[2]="not_found"
		fastqc_result_P1="no_result"
	    fi
	    if [ -n "${fastqc_files[3]}" ]
	    then
		fastq_summary_P2=$(unzip -l ${fastqc_files[3]} | grep fastqc_data.txt | awk '{ print $4 }' )
		fastqc_result_P2=$(unzip -p ${fastqc_files[3]} $fastq_summary_P2 | grep "Total Sequences" | awk '{ print $3 }' )
	    else
		fastqc_files[3]="not_found"
		fastqc_result_P2="no_result"
	    fi
	    if [ -z "$fastqc_result_R1" ]
	    then
		fastqc_result_R1="no_result"
	    fi
	    if [ -z "$fastqc_result_R2" ]
	    then
		fastqc_result_R2="no_result"
	    fi
	    if [ -z "$fastqc_result_P1" ]
	    then
		fastqc_result_P1="no_result"
	    fi
	    if [ -z "$fastqc_result_P2" ]
	    then
		fastqc_result_P2="no_result"
	    fi
	else
	    echo "WARNING: FASTQC results were not found for $sample"
	    fastqc_files=("not_found" "not_found" "not_found" "not_found" )
	    fastqc_result_R1="no_result"
	    fastqc_result_R2="no_result"
	    fastqc_result_P1="no_result"
	    fastqc_result_P2="no_result"
	fi
	
	# quast results
	if [ -n "$(find $out/quast -iname $pnusa* )" ]
	then
	    quast_result=$(ls $out/quast/$pnusa*/report.html | head -n 1 ) 
	else
	    echo "WARNING: QUAST results were not found for $sample"
	    quast_result="not_found"
	fi
	
	# seqsero results
	if [ -n "$(find $out/SeqSero -iname $pnusa* )" ] && [ -n "$salm_flag" ]
	then
	    seqsero_file=$(ls $out/SeqSero/$pnusa*.Seqsero_result.txt | head -n 1 )
	    seqsero_serotype=$(grep "Predicted serotype" $seqsero_file | cut -f 2 | tr ' ' '_' )
	    seqsero_profile=$(grep "Predicted antigenic profile" $seqsero_file | cut -f 2 | tr ' ' '_' )
	    simple_seqsero_result=$(echo $seqsero_serotype | perl -pe 's/[^\w.-]+//g' | sed 's/potentialmonophasicvariantof//g' | sed 's/potential_monophasic_variant_of_//g' | sed 's/O5-//g' )
	else
	    if [ -n "$salm_flag" ]
	    then
		echo "WARNING: SEQSERO results were not found for $sample"
	    fi
	    seqsero_file="not_found"
	    seqsero_result="no_result"
	    simple_seqsero_result="no_result"
	fi
	
	if [ -z "$salm_flag" ]
	then
	    seqsero_file="not_salmonella"
	    seqsero_serotype="not_salmonella"
	    seqsero_profile="not_salmonella"
	    simple_seqsero_result="not_salmonella"
	fi
	if [ -z "$ecoli_flag" ]
	then
	    abricate_serotypefinder_file="not_ecoli"
	    abricate_serotypefinder_result_O="not_ecoli"
	    abricate_serotypefinder_result_H="not_ecoli"
	    ecoh_file="not_ecoli"
	    ecoh_result_O="not_ecoli"
	    ecoh_result_H="not_ecoli"
	fi
	echo -e "$sample\t$pnusa\t$out\t${raw_reads[0]}\t${raw_reads[1]}\t${clean_reads[0]}\t${clean_reads[1]}\t${fastqc_files[0]}\t$fastqc_result_R1\t${fastqc_files[1]}\t$fastqc_result_R2\t${fastqc_files[2]}\t$fastqc_result_P1\t${fastqc_files[3]}\t$fastqc_result_P2\t$shovill_result\t$prokka_result\t$quast_result\t$mash_file\t$mash_result\t$simple_mash_result\t$seqsero_file\t$seqsero_serotype\t$seqsero_profile\t$simple_seqsero_result\t$ecoh_file\t$ecoh_result_O\t$ecoh_result_H\t$abricate_serotypefinder_file\t$abricate_serotypefinder_result_O\t$abricate_serotypefinder_result_H\t$resfinder_file\t$resfinder_result\t$vfdb_file\t$vfdb_result" >> $out/Run_summary.txt
    done
    
    echo "Mash results count"
    cut -f 21 $out/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" 

    echo "Seqsero results count"
    cut -f 25 $out/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" | grep -v "not_salmonella"
    
    echo "Abricate serotype results count"
    cut -f 30,31 $out/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"
    
    echo "Abricate ecoh results count"
    cut -f 27,28 $out/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"

    date
    echo "Finding each file for each sample complete!"
}
run_dependencies ()
{
    date
    for program in shovill prokka parallel sed awk seqyclean
    do
	if [ -z "$(which $program)" ]
	then
	    echo "FAILURE: $program was not found"
	    path_flag=1
	else
	    which $program
	fi
    done

    for program in fastqc quast mash abricate basemount-cmd multiqc SeqSero.py
    do
	if [ -z "$(which $program)" ]
	then
	    echo "WARN: $program was not found"
	else
	    which $program
	fi
    done

    date
    if [ -n "$path_flag" ]
    then
	echo "UPHL_reference_free.sh cannot run. Please put these in your path!!!"
	echo "Paths can be added to PATH via \"-P\""
	exit
    fi
}
run_fastqc ()
{
    date
    out=$1
    if [ -z "$out" ] || [ -z "$(find $out/Sequencing_reads/Raw -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    echo "Now beginning fastqc."
    echo "COMMAND: fastqc --outdir out/fastqc --threads 48 Sequencing_reads/*/sample.fastq "
    mkdir -p $out/fastqc
    
    ls $out/Sequencing_reads/*/*fastq* | parallel " fastqc --outdir $out/fastqc --threads 48 {} "
    
    date
    if [ -n "$(find $out/fastqc/ -iname '*.html' | head -n 1 )" ]
    then 
	echo "FastQC is complete. The results are located in $out/fastqc."
    else
	echo "Something went wrong! There are no fastq results!!!"
    fi
}
run_mash ()
{
    date
    out=$1
    mash_sketches=$2

    if [ -z "$out" ] || [ -z "$(find $out/Sequencing_reads/Raw -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    if [ ! -f "$mash_sketches" ]
    then
	echo "Cannot find refseq sketches file for mash!"
	exit
    fi

    echo "Beginning MASH analysis. First the files will be concantinated."
    echo "COMMAND: cat out/Sequencing_reads/QCed/sample_clean_PE1.fastq out/Sequencing_reads/QCed/sample_clean_PE2.fastq out/Sequencing_reads/QCed/sample_clean_SE.fastq > out/mash/sample.clean_all.fastq " 
    echo "COMMAND: mash sketch -m 2 -p 48 sample.clean_all.fastq "
    echo "COMMAND: mash dist -p 48 $mash_sketches sample.clean_all.msh > sample.distance.txt "          
    echo "COMMAND: sort -gk3 sample.distance.txt > sample.distance.txt.sorted.txt "                                                                                                             
    mkdir -p $out/mash
    ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 " cat $out/Sequencing_reads/QCed/{}_clean_PE1.fastq $out/Sequencing_reads/QCed/{}_clean_PE2.fastq $out/Sequencing_reads/QCed/{}_clean_SE.fastq > $out/mash/{}.clean_all.fastq "

    ls $out/mash/*clean_all.fastq | parallel --jobs 10 " mash sketch -m 2 -p 48 {} "
    ls $out/mash/*clean_all.*.msh | parallel --jobs 10 " mash dist -p 48 /home/Bioinformatics/Data/RefSeqSketchesDefaults.msh {} > {}.distance.txt "         
    ls $out/mash/*sh.distance.txt | parallel --jobs 10 " sort -gk3 {} > {}.sorted.txt "                                                                                                             
    
    date
    if [ -n "$(find $out/mash/ -iname '*sorted.txt' | head -n 1 )" ]
    then
	echo "MASH is complete."
	echo "The MASH results are located in $out/mash"
    else
	echo "Something went wrong! There are no MASH results!!!"
    fi
}
run_prokka ()
{
    date
    out=$1
    if [ -z "$out" ] || [ -z "$(find $out/ALL_assembled -iname '*_contigs.fa' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi
    
    echo "Running PROKKA"
    echo "COMMAND: prokka --cpu 48 --compliant --centre --UPHL --mincontiglen 500 --outdir out/Prokka/sample --locustag locus_tag --prefix sample --genus mash_serotype --species mash_serotype --force out/ALL_assembled/sample_contigs.fa" 
    sample_list=($(ls $out/ALL_assembled/*_contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 ))
    for sample in ${sample_list[@]}
    do
	mash_file=$(ls $out/mash/$sample*sorted.txt | head -n 1 )
	mash_serotype=($(head -n 1 $mash_file | cut -f 1 | awk -F "-.-" '{ print $NF }' | sed 's/.fna//g' | awk -F "_" '{ print $1 "\t" $2 }' ))
	prokka_command+=("$(echo "prokka___--cpu___48___--compliant___--centre___--UPHL___--mincontiglen___500___--outdir___$out/Prokka/$sample""___--locustag___locus_tag___--prefix___$sample""___--genus___${mash_serotype[0]}""___--species___${mash_serotype[1]}""___--force___$out/ALL_assembled/$sample""_contigs.fa" )")
    done
    history -p ${prokka_command[@]} | sort | uniq | sed 's/___/ /g' | parallel --jobs 10 '{}'
	
    mkdir -p $out/ALL_gff
    cp $out/Prokka/*/*gff $out/ALL_gff/.

    if [ -n "$(find $out/ALL_gff -iname '*.gff' | head -n 1 )" ]
    then
	date
	echo "PROKKA is complete."
	echo "Results can be found at $out/Prokka and $out/ALL_gff"
    else
	date
	echo "Something went wrong and Prokka did not run properly"
    fi
}
run_quast ()
{
    date
    out=$1
    if [ -z "$out" ] || [ -z "$(find $out/ALL_assembled -iname '*contigs.fa' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi
    echo "Now beginning QUAST."
    echo "COMMAND: quast out/ALL_assembled/sample_name_contigs.fa --output-dir out/quast/sample_name "
    mkdir -p $out/quast
    
    ls $out/ALL_assembled/*contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | parallel --jobs 10 " quast $out/ALL_assembled/{}_contigs.fa --output-dir $out/quast/{} "
    ls $out/ALL_assembled/*contigs.fa | sed 's!.*/!!' | cut -d "_" -f 1 | parallel --jobs 10 " cp $out/quast/{}/report.txt $out/quast/{}.report.txt "

    date
    if [ -n "$(find $out/quast/ -iname '*report.txt' | head -n 1 )" ]
    then 
	echo "QUAST is complete. The results are located in $out/quast."
    else
	echo "Something went wrong! There are no quast results!!!"
    fi
}
run_reorganize ()
{
    date
    out=$1
    if [ -z "$out" ] || [ -z "$(find $out/mash -iname '*sh.distance.txt.sorted.txt' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    echo "Beginning to look through mash results that may benefit from SeqSero or Abricate analysis."
    head -n 1 $out/mash/*sh.distance.txt.sorted.txt | grep -v "PNUSAE" | grep "Escherichia" | grep "coli" | awk '{ print $2 }' | sed 's!.*/!!' | cut -d "." -f 1 | sort >> $out/ecoli.txt
    head -n 1 $out/mash/*sh.distance.txt.sorted.txt | grep -v "PNUSAS" | grep "Salmonella" | awk '{ print $2 }' | sed 's!.*/!!' | cut -d "." -f 1 | sort >> $out/salmonella.txt

    if [ -s "$out/ecoli.txt" ]
    then
	echo "$out/ecoli.txt should contain non-PNUSAE samples that look like E. coli"
    else
	rm $out/ecoli.txt
    fi
    if [ -s "$out/salmonella.txt" ]
    then
	echo "$out/salmonella.txt should contain non-PNUSAS samples that look like Salmonella."
    else
	rm $out/salmonella.txt
    fi

    date
    echo "Reorganization is complete."
}
run_seqsero ()
{
    out=$1

    date
    if [ -z "$out" ] || [ -z "$(find $out/Sequencing_reads/QCed -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    if [ -z "$(which SeqSero.py)" ]
    then
	echo "SEQSERO is not found. Please add to your path with \"-P\""
	exit
    fi

    mkdir -p $out/SeqSero
    cd $out/SeqSero

    echo "Starting SeqSero to determine Salmonella subtype for PNUSAS* samples"
    echo "COMMAND: SeqSero.py -m 2 -i out/Sequencing_reads/QCed/clean_PE1.fastq out/Sequencing_reads/QCed/clean_PE2.fastq "

    ls $out/Sequencing_reads/Raw/PNUSAS*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 10 " SeqSero.py -m 2 -i $out/Sequencing_reads/QCed/{}*clean_PE1.fastq $out/Sequencing_reads/QCed/{}*clean_PE2.fastq "
    
    if [ -f "$out/salmonella.txt" ]
    then
	echo "Starting SeqSero to determine Salmonella subtype for samples in $out/salmonella.txt"
	cat $out/salmonella.txt | parallel --jobs 5 " SeqSero.py -m 2 -i $out/Sequencing_reads/QCed/{}*clean_PE1.fastq $out/Sequencing_reads/QCed/{}*clean_PE2.fastq "
    fi
    
    mv SeqSero_result* $out/SeqSero
    
    RESULTS=$(ls $out/SeqSero/*/Seqsero_result.txt)

    for result in ${RESULTS[@]}
    do
	SAMPLE=$(head -n 1 $result | awk '{ print $3 }' | awk -F "." '{ print $1 }' )
	cp $result $out/SeqSero/$SAMPLE.Seqsero_result.txt
    done

    date
    if [ -n "$(find $out/SeqSero -iname '*Seqsero_result.txt' | head -n 1 )" ]
    then
	echo "SeqSero is complete."
	echo "Results can be found at $out/SeqSero"
    else
	echo "Something went wrong and SeqSero did not complete properly"
    fi
}
run_seqyclean ()
{
    out=$1
    seqyclean_adaptors=$2

    date
    if [ -z "$out" ] || [ -z "$(find $out/Sequencing_reads/Raw -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    if [ -z "$(which seqyclean)" ] || [ ! -f "$seqyclean_adaptors" ]
    then
	echo "A path to seqyclean and adaptors is required. This can be set with \"-C\""
	exit
    fi

    echo "Now beginning SeqyClean"
    echo "COMMAND: seqyclean -minlen 25 -qual -c $seqyclean_adaptors -1 R1.fastq -2 R2.fastq -o out/Sequencing_reads/QCed/file_name_clean"
    mkdir -p $out/Sequencing_reads/Logs
    mkdir -p $out/Sequencing_reads/QCed

    ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq
    FASTQ=$(ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq )
    echo ${FASTQ[@]}
    sample_reads=()
    
    for file_name in ${FASTQ[@]}                                                                                                          
    do
	READS=($(ls $out/Sequencing_reads/Raw/$file_name*.fastq*))
	sample_reads+=("$(echo "seqyclean___-minlen___25___-qual___-c___$seqyclean_adaptors""___-1___${READS[0]}___-2___${READS[1]}___-o___$out/Sequencing_reads/QCed/$file_name""_clean" )")
    done

    history -p ${sample_reads[@]} | sed 's/___/ /g' | parallel --jobs 10 '{}'

    mv $out/Sequencing_reads/QCed/*SummaryStatistics* $out/Sequencing_reads/Logs/.

    date    
    if [ -n "$(find $out/Sequencing_reads/Logs/ -iname '*SummaryStatistics*' | head -n 1 )" ]
    then
	echo "SeqyClean has completed."
	echo "Clean fastq files are located in $out/Sequencing_reads/QCed"
	echo "The log files are located in $out/Sequencing_reads/Logs"
    else
	echo "Something went wrong. There are no clean files!!!"
    fi
}
run_shovill ()
{
    out=$1
    date
    if [ -z "$out" ] || [ -z "$(find $out/Sequencing_reads/QCed -iname '*fastq*' | head -n 1 )" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    echo "Now starting alignment with Shovill"
    echo "COMMAND: shovill --cpu 48 --ram 200 --outdir out/shovill_result/sample_name --R1 out/Sequencing_reads/QCed/sample_name_PE1.fastq --R2 out/Sequencing_reads/QCed/sample_name_PE2.fastq"
    mkdir -p $out/shovill_result
    mkdir -p $out/ALL_assembled

    ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel --jobs 2 " shovill --cpu 48 --ram 200 --outdir $out/shovill_result/{} --R1 $out/Sequencing_reads/QCed/{}_clean_PE1.fastq --R2 $out/Sequencing_reads/QCed/{}_clean_PE2.fastq "

    sample_names=($(ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq ))
    for sample_name in ${sample_names[@]}
    do
	if [ ! -f "$out/shovill_result/$sample_name/contigs.fa" ]
	then
	    sample_name_1=$(echo "$sample_name""_clean_PE1.fastq" )
	    sample_name_2=$(echo "$sample_name""_clean_PE2.fastq" )
	    shovill --cpu 48 --ram 200 --outdir $out/shovill_result/$sample_name --R1 $out/Sequencing_reads/QCed/$sample_name_1 --R2 $out/Sequencing_reads/QCed/$sample_name_2
	fi
    done

    ls $out/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | parallel " cp $out/shovill_result/{}/contigs.fa $out/ALL_assembled/{}_contigs.fa"                                   

    date
    if [ -n "$(find $out/ALL_assembled/ -iname '*contigs.fa' | head -n 1 )" ]
    then
	echo "Alignment is completed."
	echo "Aligned reads (contigs.fa) are located in $out/shovill_result and $out/ALL_assembled"
    else
	echo "Something went wrong and Shovill did not create the *.contig files. Sorry!"
    fi
}
run_wait ()
{
    run=$1
    if [ -z "$run" ] || [ ! -d "/home/BaseSpace" ]
    then
	echo "The necessary arguments (like the output directory) and/or steps are missing!"
	exit 1
    fi

    if [ -z "$(which basemount-cmd)" ]
    then
	echo "Basemount is not mounted correctly at /home/BaseSpace"
	exit
    fi

    echo "Waiting for $run in /home/BaseSpace/Runs/$run"

    cd /home/BaseSpace/Runs
    basemount-cmd refresh

    while [ ! -f "/home/BaseSpace/Runs/$run/Files/SampleSheet.csv" ]
    do
	date
	echo "Did not find SampleSheet.csv for run $run. Sleeping for 20 min."
	sleep 20m
	basemount-cmd refresh
    done
    
    date
    echo "Sample file was located at /home/BaseSpace/Runs/$run/Files/SampleSheet.csv "

    echo "Looking for a test sample in the BaseSpace directory. The test sample from this line of information:"
    LAST_SAMPLE=($(tail -n 1 /home/BaseSpace/Runs/$run/Files/SampleSheet.csv | awk -F "," '{ print $1 "\t" $2 }' ))
    TEST_SAMPLE=$(ls /home/BaseSpace/Projects/*/Samples/${LAST_SAMPLE[0]}*/Files/${LAST_SAMPLE[1]}*R1_001.fastq* | head -n 1 )
    
    cd /home/BaseSpace/Projects
    
    while [ -z "$TEST_SAMPLE" ]
    do
	date
	echo "looking for /home/BaseSpace/Projects/*/Samples/${LAST_SAMPLE[0]}*/Files/${LAST_SAMPLE[1]}*fastq*"
	echo "Did not find fastq files, yet. Sleeping for 20 min."
	sleep 20m
	
	basemount-cmd refresh
	TEST_SAMPLE=$(ls /home/BaseSpace/Projects/*/Samples/${LAST_SAMPLE[0]}*/Files/${LAST_SAMPLE[1]}*R1_001.fastq* | head -n 1 )
    done
    
    date
    echo "Run $run has completed!"
    echo "Now starting UPHL_reference_free pipeline"
}

flag=""
fastqc_flag=""
quast_flag=""
shovill_flag=""
abricate_flag=""
mash_flag=""
reorganize_flag=""
seqsero_flag=""
seqyclean_flag=""
abricate_flag=""
prokka_flag=""
multiqc_flag=""

while getopts 'i:r:o:fcmgaslxpqvbP:hC:M:' OPTION
do
    case "$OPTION" in
	h)
	    echo "$HELP"
	    exit 1
	    ;;
	b)
	    echo "The following commands will be run:"
	    echo "     SEQYCLEAN: seqyclean -minlen 25 -qual -c seqyclean_adaptors -1 R1.fastq -2 R2.fastq -o out/Sequencing_reads/QCed/file_name_clean"
	    echo "     SHOVILL:   shovill --cpu 48 --ram 200 --outdir out/shovill_result/sample_name --R1 out/Sequencing_reads/QCed/sample_name_PE1.fastq --R2 out/Sequencing_reads/QCed/sample_name_PE2.fastq"
	    echo "     PROKKA:    prokka --cpu 48 --compliant --centre --UPHL --mincontiglen 500 --outdir out/Prokka/sample --locustag locus_tag --prefix sample --genus mash_serotype --species mash_serotype --force out/ALL_assembled/sample_contigs.fa" 	    
	    echo "     MASH:      cat out/Sequencing_reads/QCed/sample_clean_PE1.fastq out/Sequencing_reads/QCed/sample_clean_PE2.fastq out/Sequencing_reads/QCed/sample_clean_SE.fastq > out/mash/sample.clean_all.fastq " 
	    echo "     MASH:      mash sketch -m 2 -p 48 sample.clean_all.fastq "
	    echo "     MASH:      mash dist -p 48 refseq_sketches sample.clean_all.msh > sample.distance.txt "          
	    echo "     MASH:      sort -gk3 sample.distance.txt                        > sample.distance.txt.sorted.txt "
	    echo "     ABRICATE:  abricate -db resfinder      out/ALL_assembled/*_contigs.fa > out/abricate_results/resfinder_out.txt "
	    echo "     ABRICATE:  abricate -db vfdb           out/ALL_assembled/*_contigs.fa > out/abricate_results/vfdb_out.txt "
	    echo "     ABRICATE:  abricate -db serotypefinder out/ALL_assembled/*_contigs.fa > out/abricate_results/Ecoli_serotypefinder_out.txt "
	    echo "     ABRICATE:  abricate -db ecoh           out/ALL_assembled/*_contigs.fa > out/abricate_results/ecoh_out.txt"
	    echo "     SEQSERO:   SeqSero.py -m 2 -i out/Sequencing_reads/QCed/clean_PE1.fastq out/Sequencing_reads/QCed/clean_PE2.fastq "
	    echo "     FASTQC:    fastqc --outdir out/fastqc --threads 48 Sequencing_reads/*/sample.fastq "
            echo "     QUAST:     quast out/ALL_assembled/sample_name_contigs.fa --output-dir out/quast/sample_name "
	    echo "     MULTIQC:   multiqc --outdir out/Logs --cl_config \"prokka_fn_snames: True\" out"
	    exit 1
	    ;;
	i)
	    input_directory=$OPTARG
	    echo "Will use the fastq files in $OPTARG : "
	    ls $OPTARG/*fastq*
	    ;;
	o)
	    output_directory=$OPTARG
	    echo "Will use the following directory for outfiles: $OPTARG"
	    mkdir -p $output_directory
	    mkdir -p $output_directory/Logs
	    ;;
	r)
	    run=$OPTARG
	    echo "Will use the files from run $OPTARG"
	    output_directory=/home/IDGenomics_NAS/WGS_Serotyping/$run
	    ;;
	f)
	    flag=1
	    fastqc_flag=1
	    echo "Not running full pipeline, but will run FASTQC"
	    ;;
	c)
	    flag=1
	    seqyclean_flag=1
	    echo "Not running full pipeline, but will run SEQYCLEAN"
	    ;;
	C)
	    seqyclean_adaptors=$OPTARG
	    if [ -f "$seqyclean_adaptors" ]
	    then
		echo "SEQYCLEAN adaptors are located at $seqyclean_adaptors"
	    else
		echo "Could not find SEQYCLEAN adaptors"
		echo "$USAGE"
		exit
	    fi
	    ;;
	m)
	    flag=1
	    mash_flag=1
	    echo "Not running full pipeline, but will run MASH"
	    ;;
	M)
	    mash_sketches=$OPTARG
	    echo "refseq sketches for mash is located at $mash_sketches"
	    ;;
	g)
	    flag=1
	    reorganize_flag=1
	    echo "Looking through mash results to find Ecoli samples that are not PNUSAE*"
	    echo "and/or Salmonella samples that are not PNUSAS*"
	    ;;
	s)
	    flag=1
	    seqsero_flag=1
	    echo "Not running full pipeline, but will run SEQSERO"
	    ;;
	l)
	    flag=1
	    shovill_flag=1
	    echo "Not running full pipeline, but will run SHOVILL"
	    ;;
	a)
	    flag=1
	    abricate_flag=1
	    echo "Not running full pipeline, but will run ABRICATE"
	    ;;
	p)
	    flag=1
	    prokka_flag=1
	    echo "Not running full pipeline, but will run PROKKA"
	    ;;
	q)
	    flag=1
	    quast_flag=1
	    echo "Not running full pipeline, but will run QUAST"	    
	    ;;
	x)
	    flag=1
	    multiqc_flag=1
	    echo "Checking files for pipeline completeness"
	    ;;
	P)
	    echo "Adding $OPTARG to PATH"
	    PATH=$PATH:$OPTARG
	    echo "Path is now"
	    echo "$PATH"
	    ;;
	v)
	    echo "WGS_NGS_pipeline v0.14"
	    fastqc --version
	    quast --version
	    shovill --version
	    mash_version=$(mash --version)
	    echo "mash $mash_version"
	    seqy_version=$(seqyclean | grep "Version")
	    echo "seqyclean $seqy_version"
	    abricate --version
	    prokka --version
	    echo "basemount-cmd no_version"
	    echo "SeqSero.py no_version"
	    multiqc --version

	    run_dependencies
	    exit 1
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

#----#----#----#----#--PIPELINE--#----#----#----#----#----

if [ -z "$flag" ]
then
    wait_flag=1
    fastqc_flag=1
    quast_flag=1
    shovill_flag=1
    abricate_flag=1
    mash_flag=1
    reorganize_flag=1
    seqsero_flag=1
    seqyclean_flag=1
    abricate_flag=1
    prokka_flag=1
    multiqc_flag=1
    date
    echo "Beginning UPHL reference-free pipeline"
fi

run_dependencies

if [ -z "$output_directory" ]
then
    echo "$USAGE"
    exit 1
fi

mkdir -p $output_directory/Sequencing_reads
mkdir -p $output_directory/Sequencing_reads/Raw
mkdir -p $output_directory/Logs

if [ ! -d "$output_directory" ]
then
    echo "Directory for results could not be created"
    exit
fi

if [ -n "$run" ] && [ -n "$wait_flag" ]
then
    run_wait "$run"
    run_basespace "$run" "$output_directory"
fi

if [ -z "$run" ] && [ -z "$flag" ]
then
    echo "The fastq files are"
    ls $input_directory/*fastq*
    cp $input_directory/*fastq* $output_directory/Sequencing_reads/Raw/. 
    echo "The fastq files have been moved to $output_directory/Sequencing_reads/Raw"
fi

if [ -z "$flag" ]
then
    date
    echo "The sample names in this pipeline are"
    ls $output_directory/Sequencing_reads/Raw/*fastq* | sed 's!.*/!!' | cut -d "_" -f 1 | sort | uniq | tee $output_directory/samples.txt 
fi

if [ -n "$seqyclean_flag" ]
then
    date
    echo "Running SEQYCLEAN"
    run_seqyclean "$output_directory" "$seqyclean_adaptors" > $output_directory/Logs/$run_date.seqyclean.log 2> $output_directory/Logs/$run_date.seqyclean.err
fi

if [ -n "$fastqc_flag" ]
then
    date
    echo "Running FASTQC"
    run_fastqc "$output_directory" > $output_directory/Logs/$run_date.fastqc.log 2> $output_directory/Logs/$run_date.fastqc.err
fi

if [ -n "$mash_flag" ]
then
    date
    echo "Running MASH"
    run_mash "$output_directory" "mash_sketches" > $output_directory/Logs/$run_date.mash.log 2> $output_directory/Logs/$run_date.mash.err
fi

if [ -n "$reorganize_flag" ]
then
    date
    echo "Looking for salmonella and E. coli samples that do not have their respsective PNUSAS* or PNUSAE* label"
    run_reorganize "$output_directory" > $output_directory/Logs/$run_date.reorganize.log 2> $output_directory/Logs/$run_date.reorganize.err
fi

if [ -n "$seqsero_flag" ]
then
    date
    echo "If there are any salmonella samples, SEQSERO will be used to serotype them"
    if [ -n "$(find $output_directory/Sequencing_reads/Raw -iname 'PNUSAS*' | head -n 1 )" ] || [ -f "$output_directory/salmonella.txt" ]
    then
	cd $output_directory
	run_seqsero "$output_directory" > $output_directory/Logs/$run_date.seqsero.log 2> $output_directory/Logs/$run_date.seqsero.err
	date
	echo "Seqsero is complete!"
    fi
fi

if [ -n "$shovill_flag" ]
then
    date
    echo "Running SHOVILL"
    run_shovill "$output_directory" > $output_directory/Logs/$run_date.shovill.log 2> $output_directory/Logs/$run_date.shovill.err
fi

if [ -n "$abricate_flag" ]
then
    date
    echo "Running ABRICATE. E coli samples will also be serotyped."
    run_abricate "$output_directory" > $output_directory/Logs/$run_date.abricate.log 2> $output_directory/Logs/$run_date.abricate.err
fi

if [ -n "$prokka_flag" ]
then
    date
    echo "Running PROKKA"
    run_prokka "$output_directory" > $output_directory/Logs/$run_date.prokka.log 2> $output_directory/Logs/$run_date.prokka.err
fi

if [ -n "$quast_flag" ]
then
    date
    echo "Running QUAST"
    run_quast "$output_directory" > $output_directory/Logs/$run_date.quast.log 2> $output_directory/Logs/$run_date.quast.err
fi

if [ -n "$multiqc_flag" ]
then
    date
    echo "Checking files for pipeline completeness"
    multiqc --outdir $output_directory/Logs --cl_config "prokka_fn_snames: True" $output_directory > $output_directory/Logs/$run_date.multiqc.log 2> $output_directory/Logs/$run_date.multiqc.err
    run_complete "$output_directory" 2> $output_directory/Logs/$run_date.summary.err | tee $output_directory/Logs/$run_date.summary.log | grep "WARNING"
    
    echo "Mash results count"
    cut -f 21    $output_directory/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" 

    echo "Seqsero results count"
    cut -f 25    $output_directory/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c | sort -k 1 -n | grep -v "no_result" | grep -v "not_salmonella"
    
    echo "Abricate serotype results count"
    cut -f 30,31 $output_directory/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"
    
    echo "Abricate ecoh results count"
    cut -f 27,28 $output_directory/Run_summary.txt | awk '{if(NR>1)print}' | sort | uniq -c |  sort -k 1 -n | grep -v "no_result" | grep -v "not_ecoli"

    date
    echo "A summary is located at $output_directory/Logs/$run_date.summary.log"
    echo "File location and results can be found at $output_directory/Run_summary.txt"
fi

date
echo "UPHL reference free pipeline complete!"
echo "Log and error files are located in $output_directory/Logs"

# TO DO: 
# MULTIQC modules
