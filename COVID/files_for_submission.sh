#/bin/bash
# ~/sandbox/UPHL/COVID/files_for_submission.sh $(pwd)

out=$1

covid_samples=$out/covid_samples.txt
covid_summary=$out/covid/summary.txt

if [ -n "$out" ] && [ -f "$covid_samples" ] && [ -f "$covid_summary"]
then
  mkdir -p $out/covid/submission_files/
else
  echo "Usage: files_for_submission.sh /home/IDGenomics_NAS/WGS_Serotyping/UT-M70330-200313"
  echo "Needs $covid_samples with tab delimiated columns of lab_accession  submission_id collection_date"
  echo "Needs $covid_summary from Illumina_covid_V3.nf pipeline"
  exit
fi

while read line
do
  sample_id=$(echo $line       | awk '{print $1}' )
  submission_id=$(echo $line   | awk '{print $2}' )
  collection_date=$(echo $line | awk '{print $3}' )
  sample_run_id=$(grep $sample_id $out/covid_samples.txt)

  # getting the consensus fasta file
  # changing the fasta header
  cat $out/covid/consensus/$sample_run_id.consensus.fa | sed "s/>.*/>$submission_id/g" > $out/covid/submission_files/$submission_id.consensus.fa
##  # replacing all degenerate bases with N, removing leading Ns, and folding sequence to 75 bp wide
  # removing leading Ns, and folding sequencing to 75 bp wide
  echo ">$submission_id" > $out/covid/submission_files/$submission_id.nondegenerate.fa
#  grep -v ">" $out/covid/submission_files/$submission_id.consensus.fa | sed 's/[BDEFHIJKLMOPQRSUVWXYZ]/N/g' | sed 's/^N*N//g' | fold -w 75 >> $out/covid/submission_files/$submission_id.nondegenerate.fa
  grep -v ">" $out/covid/submission_files/$submission_id.consensus.fa | sed 's/^N*N//g' | fold -w 75 >> $out/covid/submission_files/$submission_id.nondegenerate.fa
  cat $out/covid/submission_files/$submission_id.nondegenerate.fa | awk -v id=$submission_id -v cldt=$collection_date '{ if ($0 ~ ">") print $0 " [organism=Severe acute respiratory syndrome coronavirus 2][isolate=SARS-CoV-2/human/USA/" id "/2020][host=Human][country=USA][collection_date=" cldt "]" ; else print $0 }' > $out/covid/submission_files/$submission_id.genbank.fa

  # copying fastq files and changing the file name
  cp $out/Sequencing_reads/Raw/$sample_run_id*R1_001.fastq.gz $out/covid/submission_files/$submission_id.R1.fastq.gz
  cp $out/Sequencing_reads/Raw/$sample_run_id*R2_001.fastq.gz $out/covid/submission_files/$submission_id.R2.fastq.gz

  num_n=$(grep $covid_summary | cut -f 7 -d ',')
  run_id=${sample_run_id: -6}

  # preparing fasta for gisaid submission
  if [ $num_n < 14952 ]
  then
    cat $out/covid/submission_files/$submission_id.consensus.fa >> $out/covid/submission_files/$run_id.gisaid_submission.fasta
    # preparing fasta for genbank submission
    if [ $num_n < 4903 ]
    then
      cat $out/covid/submission_files/$submission_id.genbank.fa >> $out/covid/submission_files/$run_id.genbank_submission.fasta
    fi
  fi
done < $covid_samples
