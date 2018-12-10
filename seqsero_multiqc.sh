#!/bin/bash

out=$1
echo -e "Sample\tInput_files\tO_antigen_prediction\tH1_antigen_prediction(fliC)\tH2_antigen_prediction(fljB)\tPredicted_antigenic_profile\tPredicted_serotype(s)" > $out/SeqSero/Seqsero_serotype_results.txt
RESULTS=$(ls $out/SeqSero/*/Seqsero_result.txt)
for result in ${RESULTS[@]}
do
  SAMPLE=$(head -n 1 $result | awk '{ print $3 }' | awk -F "_" '{ print $1 }' )
  seqsero_inputfile=$(grep "Input files"                 $result | cut -f 2 | tr ' ' '_' )
  seqsero_Oantipred=$(grep "O antigen prediction"        $result | cut -f 2 | tr ' ' '_' )
  seqsero_H1antpred=$(grep "H1 antigen prediction(fliC)" $result | cut -f 2 | tr ' ' '_' )
  seqsero_H2antpred=$(grep "H2 antigen prediction(fljB)" $result | cut -f 2 | tr ' ' '_' )
  seqsero_antigenic=$(grep "Predicted antigenic profile" $result | cut -f 2 | tr ' ' '_' )
  seqsero_serotypes=$(grep "Predicted serotype(s)"       $result | cut -f 2 | tr ' ' '_' )
  echo -e "$SAMPLE\t$seqsero_inputfile\t$seqsero_Oantipred\t$seqsero_H1antpred\t$seqsero_H2antpred\t$seqsero_antigenic\t$seqsero_serotypes" >> $out/SeqSero/Seqsero_serotype_results_all.txt
  grep -v "O--" $out/SeqSero/Seqsero_serotype_results_all.txt >> $out/SeqSero/Seqsero_serotype_results.txt
done
