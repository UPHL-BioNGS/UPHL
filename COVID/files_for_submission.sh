#/bin/bash

# ~/sandbox/UPHL/COVID/files_for_submission.sh $(pwd)
# out/covid/submission_files/sample_submission_id.txt should be lab_accession\tSubmission_ID\tCollection_Date

out=$1
if [ -n "$out" ] && [ -f "$out/covid_samples.txt" ] && [ -f "$out/covid/submission_files/sample_submission_id.txt" ]
then
  mkdir -p $out/covid/submission_files/
else
  echo "Usage: files_for_submission.sh /home/IDGenomics_NAS/WGS_Serotyping/UT-M70330-200313"
  echo "Needs out/covid_samples.txt and out/covid/submission_files/sample_submission_id.txt"
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
  # replacing all degenerate bases with N, removing leading Ns, and folding sequence to 75 bp wide
  echo ">$submission_id" > $out/covid/submission_files/$submission_id.nondegenerate.fa
  grep -v ">" $out/covid/submission_files/$submission_id.consensus.fa | sed 's/[BDEFHIJKLMOPQRSUVWXYZ]/N/g' | sed 's/^N*N//g' | fold -w 75 >> $out/covid/submission_files/$submission_id.nondegenerate.fa
  cat $out/covid/submission_files/$submission_id.nondegenerate.fa | awk -v id=$submission_id -v cldt=$collection_date '{ if ($0 ~ ">") print $0 " [organism=Severe acute respiratory syndrome coronavirus 2][isolate=SARS-CoV-2/human/USA/" id "/2020][host=Human][country=USA][collection_date=" cldt "]" ; else print $0 }' > $out/covid/submission_files/$submission_id.genbank.fa
  # copying fastq files and changing the file name
  cp $out/Sequencing_reads/Raw/$sample_run_id*R1_001.fastq.gz $out/covid/submission_files/$submission_id.R1.fastq.gz
  cp $out/Sequencing_reads/Raw/$sample_run_id*R2_001.fastq.gz $out/covid/submission_files/$submission_id.R2.fastq.gz
  # extracting human reads from fastq file that need them and changing the name
  if [ -n "$(grep $sample_run_id $out/covid/summary.txt | cut -f 2 -d ',' | grep -v none)" ]
  then
     cp $out/covid/bwa/$sample_run_id.sorted.bam $out/covid/submission_files/$submission_id.bam
     samtools fastq -1 $out/covid/submission_files/$submission_id.filtered.R1.fastq.gz -2 $out/covid/submission_files/$submission_id.filtered.R2.fastq.gz $out/covid/submission_files/$submission_id.bam
  fi
done < $out/covid/submission_files/sample_submission_id.txt

# creating mulifasta
cat $out/covid/submission_files/*.consensus.fa > $out/covid/submission_files/all.consensus.fasta
