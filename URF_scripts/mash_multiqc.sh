#!/bin/bash

out=$1
organisms=($(cat $out/mash/*sorted.txt | awk '{ if( $4 == 0 ) print $1 }' | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d "_" | cut -f 1 -d "." | sort | uniq -c | sort -rhk 1,1 | awk '{ print $2 }' ))
header="Sample"

for organism in ${organisms[@]}
do
  header=$(echo "$header\t$organism" )
done
echo -e "$header" > $out/mash/mash_results.txt

mash_results=($(ls $out/mash/*sorted.txt | sed 's!.*/!!' | cut -d "_" -f 1 | cut -d '.' -f 1 | sort | uniq ))
for mash_result in ${mash_results[@]}
do
  mash_line=$mash_result
  for organism in ${organisms[@]}
  do
    number=$(cat $out/mash/$mash_result*sorted.txt | awk '{ if( $4 == 0 ) print $1 }' | cut -f 8 -d "-" | sed 's/^_\(.*\)/\1/' | cut -f 1,2 -d "_" | cut -f 1 -d "." | sort | uniq -c | grep $organism | awk '{ print $1 }' )
    if [ -z "$number" ]
    then
      number="0"
    fi
    mash_line=$(echo "$mash_line\t$number" )
  done
echo -e "$mash_line" >> $out/mash/mash_results.txt
done
